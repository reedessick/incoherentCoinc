description = "a module that houses methods that represent our observations"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

import numpy as np

#-------------------------------------------------
### general utility methods

def poisson_dt( r ):
    """
    randomly draws time-between-events for poisson series with rate="r"
    """
    return -np.log(1-np.random.rand())/r

def poisson_series( T, r ):
    """
    returns a list of Trigger objects drawn from poisson series of rate="r" and duration="T"
    """
    t = 0
    ans = []
    while t < T:
        t += poisson_dt( r )
        if t > T: 
            break
        ans.append( Trigger( t ) )
    return ans

#-------------------------------------------------
### classes to wrap around data

class Trigger(object):
    """
    wrapper around data for a single trigger
    """

    def __init__(self, time, amp=None):
        self.time = time
        self.amp = amp

class Realization(object):
    """
    one realization of the observed triggers
    can either be fed lists of triggers or it will generate them itself.
    """

    def __init__(self, duration, Ndet=2):
        assert isinstance(Ndet, int), "Ndet must be an integer"
        assert Ndet > 1, "must have at least 2 detectors"

        self.Ndet = Ndet
        self.dur = duration
        self.trgs = [[]]*Ndet ### instantiate empty lists initially

    def simulate(self, rateC, rates):
        """
        simulate a realization with the associated parameters
          rateC is the rate of coincident events
          rates = list with length = Ndet with associated rates for each detector
          dur = duration of experiment
        """
        assert len(rates) == self.Ndet, "len(rates) must be equal to self.Ndet"

        s = poisson_series( self.dur, rateS ) ### coincident stream 
        for d, r in enumerate(rates):
            strm = poisson_series( self.dur, r ) + s ### noise in stream "d" and coinc events
            strm.sort(key=lambda l: l.time)
            self.trgs[d] = strm

    def get_dn(self, n):
        """
        number of events seen in detector "n"
        """
        return len(self.trgs[n])

    def get_nC(self, tau):
        """
        number of coincidences with window="tau"
        """
        Ns = [self.get_dn(n) for n in xrange(self.Ndet)] ### number of triggers in each stream
        inds = [[0,0] for _ in xrange(self.Ndet-1)] ### the indecies for each stream (besides stream 0, which we treat separately)

        count = 0
        for trg in self.trgs[0]: ### iterate over triggers in "stream 0", counting the total number of coincs for each
            t = trg.time ### time of interest

            for n, (s, e) in enumerate(inds): ### iterate through the rest of the streams and find the idecies that bracket "t"
                N = Ns[n] ### number of events in "stream n"

                while (s<N) and (self.trgs[n][s].time < t-tau): ### start is too early, so we bump it
                    s += 1

                if e < s: ### end must be at least as big as the start
                    e = s 
                inds[n] = [s, e] ### update start and end points

                ### ensure that the start is still a good coinc
                if self.trgs[n][s].time > t+tau:
                    break ### no good coincs here! so we stop working and move to the next "steam 0" trigger

                while (e+1<N) and (self.trgs[n][e+1].time < t+tau): ### end time is too early
                    e += 1
                inds[n][1] = e ### update end point

            else: ### there must be at least one participating trigger in all streams, otherwise we would have broken.
                c = 1
                for s, e in inds:
                    c *= (e-s+1) ### total number of coincs is the product of the number participating in each stream
                count += c

        return count

    def get_nP(self, tau):
        """
        number of coincidences from all possible time slides keeping zero-lag coincidences with window="tau"
        """
        raise NotImplementedError

    def get_nM(self, tau):
        """
        number of coincidences from all possible time slides removing zero-lag coincidences with window="tau"
        """
        raise NotImplementedError
