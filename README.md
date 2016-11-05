# incoherentCoinc

the idea here is to use the "Realization" class to simulate experiments repeatedly. We can then call that class's methods to extract statistics for each realization. We feed these into a "Likelihood" object which will then delegate to Realization to compute statistics in a predictable vectorize representation. The Likelihood is callable and returns logLikelihood using the data and specified parameters

----------------------------------------------------------------------------------------------------

Need to finish implementing "Likelihood" means, covar, and logLikelihood

specifically, need to implement:
  - mean{nM}
  - all of the covariance matrix
