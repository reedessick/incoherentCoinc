# incoherentCoinc

the idea here is to use the "Realization" class to simulate experiments repeatedly. We can then call that class's methods to extract statistics for each realization. We feed these into a "Data" object which represents all the observations from a single realization. Our likelihood will then take in a "Data" object and the associated rates.

Currently in the middle of implementing the "Realization" class...
