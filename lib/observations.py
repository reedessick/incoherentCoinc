description = "a module that houses methods that represent our observations"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

import numpy as np
import copy

import sympy
import itertools

#-------------------------------------------------
### combinatorics
#-------------------------------------------------

def kPartitions(M, k):
    """Returns a list of all unique k-partitions of range(M).

    Each partition is a list of parts, and each part is a tuple.

    taken with slight modifications from:
        http://stackoverflow.com/questions/39192777/python-split-a-list-into-n-groups-in-all-possible-combinations-of-group-length-a
    """
    groups = []  # a list of lists, currently empty
    seq = range(M) ### indecies which we group

    def generate_partitions(i):
        if i >= M: ### we're out of elements, so just map them to tuples
            yield list(map(tuple, groups))

        else: ### still elements that are unassigned
            if M - i > k - len(groups): ### there are enough elements left to still get enough groups without them all be new groups
                for group in groups:
                    group.append(seq[i])
                    for x in generate_partitions(i + 1): ### recurse
                        yield x
                    group.pop() ### not exactly sure what this is doing, but it works.

            if len(groups) < k: ### we should add another group because we don't have enough
                groups.append([seq[i]])
                for x in generate_partitions(i + 1): ### recurse
                    yield x
                groups.pop() ### not exactly sure what this is doing, but it works

    return generate_partitions(0)

#-------------------------------------------------
### related to Poisson Distribution
#-------------------------------------------------

def poisson_dt( r ):
    """
    randomly draws time-between-events for poisson series with rate="r"
    """
    if r == 0:
        return np.infty
    else:
        return -np.log(1-np.random.rand())/r

#------------------------

def poisson_series( T, r, amp=None, label=None ):
    """
    returns a list of Trigger objects drawn from poisson series of rate="r" and duration="T"
    """
    t = 0
    ans = TrigList()
    while t < T:
        t += poisson_dt( r )
        if t > T: 
            break
        ans.append( Trigger( t, amp=amp, label=label ) )
    return ans

#------------------------

def poisson_moment( k, n):
    """
    returns the moment of x**n with expectation value k

    CURRENTLY A SET OF HARD CODED EXPRESSIONS! VERY FRAGILE!
      --> would be *much* better if we could do this algorithmically
    """
    if n==0:
        return 1

    elif n==1:
        return k

    elif n==2:
        return k**2 + k

    elif n==3:
        return k**3 + 3*k**2 + k

    elif n==4:
        return k**4 + 6*k**3 + 7*k**2 + k

    elif n==5:
        return k**5 + 10*k**4 + 25*k**3 + 15*k**2 + k

    elif n==6:
        return k**6 + 15*k**5 + 65*k**4 + 90*k**3 + 31*k**2 + k

    elif n==7:
        return k**7 + 21*k**6 + 140*k**5 + 350*k**4 + 301*k**3 + 63*k**2 + k

    elif n==8:
        return k**8 + 28*k**7 + 266*k**6 + 1050*k*85 + 1701*k**4 + 966*k**3 + 127*k**2 + k

    else:
        raise NotImplementedError('currently only support n<=8')

#------------------------

def poisson_expect( expr, symbs, cnts ):
    """
    returns the numerical expectation value of a combination of poisson distributed counts

    assumes expression is polynomial in the rates (really, rates=rate*tau~dimensionless)
    """
    assert len(symbs)==len(cnts), "symbs and rates must have the same dimension"

    if len(symbs): ### we still have at least one symbol to substitute
        coeffs = sympy.poly_from_expr(expr, symbs[0])[0].all_coeffs() ### make into a polynomial and extract coefficients
        dim = len(coeffs)-1 ### dimension+1 of the polynomial

        ans = 0.0 ### holder for the answer
        for n, coeff in enumerate(coeffs):
            ### add to this the expectation value of the coefficient * the associated moment of this variable
            ans += poisson_expect( coeff, symbs[1:], cnts[1:] ) * poissson_moment( cnts[0], dim-n ) 

    else: ### we're at the bottom, so just return the expression cast as a float
        ans = float(expr)

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

        for n, (s, e) in enumerate(inds, 1): ### iterate through the rest of the streams and find the idecies that bracket "t"
            N = Ns[n] ### number of events in "stream n"

            while (s<N) and (trgs[n][s].time < t-tau): ### start is too early, so we bump it
                s += 1

            if e < s: ### end must be at least as big as the start
                e = s
            inds[n-1][0] = s ### update start points

            ### ensure that the start is still a good coinc
            # too large      too far away
            if (s>=N) or (trgs[n][s].time > t+tau):
                inds[n-1][1] = e ### update end point
                break ### no good coincs here! so we stop working and move to the next "steam 0" trigger

            while (e+1<N-1) and (trgs[n][e+1].time < t+tau): ### end time is too early
                e += 1
            inds[n-1][1] = e ### update end point

        else: ### there must be at least one participating trigger in all streams, otherwise we would have broken.
            c = 1
            for s, e in inds:
                c *= (e-s+1) ### total number of coincs is the product of the number participating in each stream
            count += c

    return count

#------------------------

def removeCoinc(trgs, tau):
    """
    rinds all coincs between all triggers in trgs and removes them from the lists
    assumes chronological ordering

    returns trgs including modifications (for symantic ease)
    """
    Ndet = len(trgs)
    Ns = [len(_) for _ in trgs]
    inds = [[0,0] for _ in xrange(Ndet-1)]

    remove = [set() for _ in xrange(Ndet)] ### list of indicies to be removed

    for ind, trg in enumerate(trgs[0]):
        t = trg.time ### time of interest

        for n, (s, e) in enumerate(inds, 1): ### iterate through the rest of the streams and find the idecies that bracket "t"
            N = Ns[n] ### number of events in "stream n"

            while (s<N) and (trgs[n][s].time < t-tau): ### start is too early, so we bump it
                s += 1

            if e < s: ### end must be at least as big as the start
                e = s
            inds[n-1][0] = s ### update start points

            ### ensure that the start is still a good coinc
            # too big        too far away
            if (s>=N) or (trgs[n][s].time > t+tau):
                inds[n-1][1] = e ### update end point
                break ### no good coincs here! so we stop working and move to the next "steam 0" trigger

            while (e+1<N) and (trgs[n][e+1].time < t+tau): ### end time is too early
                e += 1
            inds[n-1][1] = e ### update end point

        else: ### there must be at least one participating trigger in all streams, otherwise we would have broken.

            remove[0].add(ind)
            for rem, (s, e) in zip(remove[1:], inds): ### add indecies to be removed to sets
                for i in xrange(s,e+1):
                    rem.add(i)
            
    for n, rem in enumerate(remove): ### iterate through lists and remove indecies
        for i in sorted(rem, reverse=True): ### start from the back and work forward
            trgs[n].pop(i)

    return trgs

#------------------------

def slideCoinc(trgs, tau, dt, dur, ind=1):
    """
    perform a cyclic time shift and computes coincs for a list of TrigLists
    assumes time ordered trgs and maintains this

    recurse and pass the index of which TrigList in trgs should be slid.
    we then slide that by one and then recurse to a lower level to figure out the coincs.
    At the bottom level, we slide and count coincs without another recursive call.
    If we sum coincs back up the recursive chain, we should get the total number of coincs

    NOTE: ind defaults to starting at 1 because we don't need to slide the first row

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
    def __init__(self, time, amp=None, label=None):
        self.time = time
        self.amp = amp
        self.label = label

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
        while (i<N) and (self[i].time+dt < dur): ### don't wrap around yet so just modify in place
            self[i].time += dt
            i += 1
        j = i
        while (i<N): ### does wrap around so we manipulate the placement in the list
            trigger = self.pop(i)
            trigger.time += dt - dur
            self.insert(i-j, trigger)
            i += 1

    def __add__(self, other):
        '''
        add in triggers so they're in sequence
        returns a *new* object
        '''
        new = TrigList()
        i = 0
        j = 0
        while (i<len(self)) and (j<len(other)):
            if self[i].time < self[j].time:
                new.append(self[i])
                i += 1

            else:
                new.append(self[j])
                j += 1

        while i<len(self):
            new.append(self[i])
            i += 1

        while j<len(other):
            new.append(other[j])
            j += 1

        return new
            
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

        s = poisson_series( self.dur, rateC, label='coinc' ) ### coincident stream 
        for d, r in enumerate(rates):
            strm = poisson_series( self.dur, r, label='det-%d'%d ) + s ### noise in stream "d" and coinc events
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
        assert dt < self.dur,  "slide must be smaller than experiment duration"

        return slideCoinc( copy.deepcopy(self.trgs), tau, dt, self.dur )
        
    def get_nM(self, tau, dt):
        """
        number of coincidences from all possible time slides removing zero-lag coincidences with window="tau"
        """
        assert dt >= tau, "time slide window must be larger than coinc window"
        assert dt < self.dur,  "slide must be smaller than experiment duration"

        return slideCoinc( removeCoinc(copy.deepcopy(self.trgs), tau), tau, dt, self.dur )

#------------------------

class Likelihood(object):
    """
    a representation of statistics calculated from a Realization object
    callable, which will return the likelihood 
    """
    def __init__(self, realization, taus=[]):
        self.realization = realization

        ### ensure all taus are integer multiples of eachother
        taus = sorted(taus)
        for ind, tau1 in enumerate(taus):
            for tau2 in taus[ind+1:]:
                assert tau2%tau1 <= 1e-10*tau1
        self.taus = taus

        self.Ndim = realization.Ndet + 3*len(taus) ### dimensionality of the likelihood

        self.vector = None

    def vectorize(self):
        """
        returns a vectorized list of inputs with standard format.
        *must* be called before logLikelihood. Otherwise, will be called within logLikelihood
        """
        self.vector = np.empty((self.Ndim,), dtype=float)

        ind = self.realization.Ndet
        self.vector[:ind] = [self.realization.get_dn(n) for n in xrange(ind)] ### get the number of events in each stream
        for tau in self.taus: ### iterate over taus and compute for each

            dt = tau ### FIXME: this may be a fragile/bad choice but it's consistent with what my theory work predicts

            self.vector[ind] = self.realization.get_nC(tau)
            ind += 1

            self.vector[ind] = self.realization.get_nP(tau, dt)
            ind += 1

            self.vector[ind] = self.realization.get_nM(tau, dt)
            ind += 1

        self.vector
        return self.vector

    def means(self, rateC, rates):
        """
        returns the means of the statistics
        """
        m = np.zeros((self.Ndim,), dtype=float) ### array to hold numerical results
        Ndet = self.realization.Ndet

        ### set up variables
        allRates = np.array([rateC]+rates) ### array of numerical rates

        symbs = sympy.symbols('rC '+' '.join('r%d'%n for n in xrange(Ndet))) ### sympy variables
        rC = symbs[0]
        rs = symbs[1:]

        ### compute expected number in each detector separately
        m[:Ndet] = [self.realization.dur*(rates[n]+rateC) for n in xrange(Ndet)]

        ### compute expected number of coincidences
        ind = Ndet ### counter used within next loop
        for tau in self.taus:
            N = self.realization.dur/tau ### FIXME: this may be fragile/bad choice!
            cnts = allRates*tau

            #------------
            ### nC
            #------------
            ### construct equation in terms of symbolic variables
            nC = N
            for r in rs:
                nC *= (rC+r)
            m[ind] = poisson_expect( nC, symbs, cnts )
            ind += 1

            #------------
            ### nP
            #------------
            prefact = 1.0
            terms = [rC+r for r in rs] ### the possible terms that are grouped together by the sum
            for m in xrange(Ndet):
                prefact *= (N-m) ### increment prefactor

                ans = 0
                for combo in kPartitions( Ndet, m+1 ): ### all possible partitions of Ndet terms into m groups
                    term = 1
                    for termInds in combo: ### compute expectation value for each term separately
                        expr = 1
                        for termInd in termInds:
                            expr *= terms[termInd]
                        term *= poisson_expect( expr, symbs, cnts ) ### multiply them together
                    ans += term ### add to total

                m[ind] += prefact * ans ### add to mean
            ind += 1

            #------------
            ### nM
            #------------
            raise NotImplementedError('need to compute mean{nM} for an arbitrary number of detectors')
            ### fact that we require some counts to be zero means we'll have to rethink this a bit...
            ind += 1

        return m

    def covar(self, rateC, rates):
        """
        returns the covariance matrix
        """
        c = np.zeros((self.Ndim,self.Ndim), dtype=float)

        raise NotImplementedError

        return c

    def logLikelihood(self, rateC, rates):
        """
        returns ln(likelihood) of this data given rates
        """
        if self.vector==None: ### has not been computed yet
            v = self.vectorize()
        else: ### has been computed, so we reference stored values (faster and ~cheap)
            v = self.vector

        ### compute model-dependent parameters
        m = self.means(rateC, rates)
        c = self.covar(rateC, rates)

        ### actually compute the logLikelihood
        raise NotImplementedError('actually compute the logLikelihood via matrix inversion and multiplication')
