description = "a module that houses methods that represent our observations"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

import numpy as np
import copy

#-------------------------------------------------
### related to Poisson Distribution
#-------------------------------------------------

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
    ans = TrigList()
    while t < T:
        t += poisson_dt( r )
        if t > T: 
            break
        ans.append( Trigger( t ) )
    return ans

#-------------------------------------------------
### coincidence counting and such
#-------------------------------------------------

def coinc(trgs, tau):
    """
    find coincs between all triggers in trgs ( list of lists of Triggers )
    assumes triggers are already chronologically ordered

    returns the number of coincidences, does not modify trgs
    """
    Ndet = len(trgs)
    Ns = [len(_) for _ in trgs] ### number of triggers in each stream
    inds = [[0,0] for _ in xrange(Ndet-1)] ### the indecies for each stream (besides stream 0, which we treat separately)

    count = 0
    for trg in trgs[0]: ### iterate over triggers in "stream 0", counting the total number of coincs for each
        t = trg.time ### time of interest

        for n, (s, e) in enumerate(inds): ### iterate through the rest of the streams and find the idecies that bracket "t"
            N = Ns[n] ### number of events in "stream n"

            while (s<N) and (trgs[n][s].time < t-tau): ### start is too early, so we bump it
                s += 1

            if e < s: ### end must be at least as big as the start
                e = s
            inds[n][0] = s ### update start points

            ### ensure that the start is still a good coinc
            if trgs[n][s].time > t+tau:
                inds[n][1] = e ### update end point
                break ### no good coincs here! so we stop working and move to the next "steam 0" trigger

            while (e+1<N) and (trgs[n][e+1].time < t+tau): ### end time is too early
                e += 1
            inds[n][1] = e ### update end point

        else: ### there must be at least one participating trigger in all streams, otherwise we would have broken.
            c = 1
            for s, e in inds:
                c *= (e-s+1) ### total number of coincs is the product of the number participating in each stream
            count += c

    return count

def removeCoinc(trgs, tau):
    """
    rinds all coincs between all triggers in trgs and removes them from the lists
    assumes chronological ordering

    returns trgs including modifications (for symantic ease)
    """
    Ndet = len(trgs)
    Ns = [len(_) for _ in trgs]
    inds = [[0,0] for _ in xrange(Ndet)]

    remove = [set() for _ in xrange(Ndet)] ### list of indicies to be removed

    for ind, trg in enumerate(trgs[0]):
        t = trg.time ### time of interest

        for n, (s, e) in enumerate(inds): ### iterate through the rest of the streams and find the idecies that bracket "t"
            N = Ns[n] ### number of events in "stream n"

            while (s<N) and (trgs[n][s].time < t-tau): ### start is too early, so we bump it
                s += 1

            if e < s: ### end must be at least as big as the start
                e = s
            inds[n][0] = s ### update start points

            ### ensure that the start is still a good coinc
            if trgs[n][s].time > t+tau:
                inds[n][1] = e ### update end point
                break ### no good coincs here! so we stop working and move to the next "steam 0" trigger

            while (e+1<N) and (trgs[n][e+1].time < t+tau): ### end time is too early
                e += 1
            inds[n][1] = e ### update end point

        else: ### there must be at least one participating trigger in all streams, otherwise we would have broken.

            rem[0].add(ind)
            for rem, (s, e) in zip(remove[1:], inds): ### add indecies to be removed to sets
                for i in xrange(s,e+1):
                    rem.add(i)
            
    for n, rem in enumerate(remove): ### iterate through lists and remove indecies
        for i in sorted(rem, reverse=True): ### start from the back and work forward
            trgs[n].pop(i)

    return trgs

def slideCoinc(trgs, tau, dt, dur, ind=0):
    """
    perform a cyclic time shift and computes coincs for a list of TrigLists
    assumes time ordered trgs and maintains this

    recurse and pass the index of which TrigList in trgs should be slid.
    we then slide that by one and then recurse to a lower level to figure out the coincs.
    At the bottom level, we slide and count coincs without another recursive call.
    If we sum coincs back up the recursive chain, we should get the total number of coincs

    returns the number of coincidences
    """
    if ind < len(trgs)-1: ### not the last one, so slide list and then recurse
        count = 0
        cpTrgs = copy.deepcopy(trgs[ind:]) ### make a copy of everything "below" ind so we don't mess with original trgs throughout the recursion
        for _ in xrange(int(dur/dt)): ### iterate over slides
            cpTrgs[0].slide(dt, dur) ### slide the first element
            count += slideCoinc(trgs[:ind]+cpTrgs, tau, dt, dur, ind=ind+1) ### recurse
        return count

    else: ### this is the last one, so we just count coincs
        count = 0
        cpTrgs = copy.deepcopy(trgs[ind])
        for _ in xrange(int(dur/dt)): ### iterate over slides
            cpTrgs.slide(dt, dur) ### slide 
            count += coinc(trgs[:ind]+[cpTrgs], tau) ### count coincs
        return count

#-------------------------------------------------
### classes to wrap around data
#-------------------------------------------------

class Trigger(object):
    """
    wrapper around data for a single trigger
    """
    def __init__(self, time, amp=None):
        self.time = time
        self.amp = amp

#------------------------

class TrigList(list):
    """
    represents a time-ordered series of triggers

    NOTE: be careful with duration because you don't check for sensible input anywhere...
    duration should be bigger than the largest time of any Trigger in the list.

    you also don't check whether the items are actually Trigger objects...
    """
    def slide(self, dt, dur):
        """
        perform a cyclic time shift
        """
        N = len(self)
        i = 0
        while (i<N) and (trg[i]+dt < dur): ### don't wrap around yet so just modify in place
            self.trgs[i].time += dt
            i += 1
        j = i
        while (i<N): ### does wrap around so we manipulate the placement in the list
            trigger = self.trgs.pop(i)
            trigger.time += dt - dur
            self.trgs.insert(i-j, trigger)
            i += 1

#------------------------

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
        self.trgs = [None]*Ndet ### instantiate empty lists initially

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
        return coinc(self.trgs, tau)

    def get_nP(self, tau, dt):
        """
        number of coincidences from all possible time slides keeping zero-lag coincidences with window="tau"
        """
        assert dt >= tau, "time slide window must be larger than coinc window"
        assert dt < dur,  "slide must be smaller than experiment duration"

        return slideCoinc( copy.deepcopy(self.trgs), tau, dt, self.dur )
        
    def get_nM(self, tau):
        """
        number of coincidences from all possible time slides removing zero-lag coincidences with window="tau"
        """
        assert dt >= tau, "time slide window must be larger than coinc window"
        assert dt < dur,  "slide must be smaller than experiment duration"

        return slideCoinc( removeCoinc(copy.deepcopy(self.trgs), tau), tau, dt, self.dur )

#------------------------

class Data(object):
    """
    a representation of statistics calculated from a single realization
    """
    def __init__(self, dur):
        self.dur = dur
