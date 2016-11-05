#!/usr/bin/python
usage       = "simulate.py [--options]"
description = "simulate streams of poisson series and regression"
author      = "reed.essick@ligo.org"

#-------------------------------------------------

import time

import observations

from optparse import OptionParser

#-------------------------------------------------

parser = OptionParser(usage=usage, description=description)

### verbosity and output

parser.add_option('-v', '--verbose', default=False, action='store_true')

### parameters about the number of simulations

parser.add_option('-N', '--Ntrials', default=1, type='int',
  help='the number of repetitions of the experiment')

### parameters of each simulation

parser.add_option('-r', '--rate', default=[], type='float', action='append', 
  help='specify a noise rate for a detector. The number of detectors simulated \
  will be taken as the total number of "--rate" options provided.')

parser.add_option('-R', '--rateC', default=0, type='float', 
  help='the rate of the coincident stream')

parser.add_option('-T', '--duration', default=100, type='float',
  help='the duration of the experiment. DEFAULT=100')

parser.add_option('-t', '--tau', default=[], type='float', action='append',
  help='the coincidence window to use. Can be repeated. If not supplied, will not compute coincidences in likelihood')

opts, args = parser.parse_args()

#-------------------------------------------------

### instantiate likelihood and realization
realization = observations.Realization( opts.duration, Ndet=len(opts.rate) )
likelihood = observations.Likelihood( realization, taus=opts.tau ) 
for m in xrange(1, opts.Ntrials+1):

    print "trial : %d"%m

    to = time.time()
    realization.simulate( opts.rateC, opts.rate ) ### draw from simulation
    print "    realization.simulate : %.6f sec"%(time.time()-to)

    for n, rate in enumerate(opts.rate):
        to = time.time()
        print "    d(%d)\t\t"%n, realization.get_dn(n)
        print "        realization.get_dn : %.6f sec"%(time.time()-to)

    for tau in opts.tau:

        dt = tau

        to = time.time()
        print "    nC(%.3f)\t\t"%tau, realization.get_nC(tau)
        print "        realization.get_nC : %.6f sec"%(time.time()-to)

        to = time.time()
        print "    nP(%.3f, %.3f)\t"%(tau, tau), realization.get_nP(tau, dt)
        print "        realization.get_nP : %.6f sec"%(time.time()-to)

        to = time.time()
        print "    nM(%.3f, %.3f)\t"%(tau, tau), realization.get_nM(tau, dt)
        print "        realization.get_nM : %.6f sec"%(time.time()-to)

    to = time.time()
    print "    vect\t\t %s"%(" ".join("%d"%_ for _ in likelihood.vectorize()))
    print "        likelihood.vectorize : %.6f sec"%(time.time()-to)
