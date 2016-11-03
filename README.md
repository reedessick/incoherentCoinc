# incoherentCoinc

the idea here is to use the "Realization" class to simulate experiments repeatedly. We can then call that class's methods to extract statistics for each realization. We feed these into a "Data" object which represents all the observations from a single realization. Our likelihood will then take in a "Data" object and the associated rates.

We should write some wrapper function that takes in a "Realization" and spits out a "Data" object, ensuring the "Data" object is built correctly and all that jazz. "Data" objects might need a reference to a "Realization", and then can compute the statistics they need on the fly...

Currently in the middle of figuring out how to perform time-slides and count coincs efficiently for an arbitrary unspecified number of detectors... (~lib/observations.slideCoinc(trgs, tau, dt, dur))
