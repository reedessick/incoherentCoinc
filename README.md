# incoherentCoinc

the idea here is to use the "Realization" class to simulate experiments repeatedly. We can then call that class's methods to extract statistics for each realization. We feed these into a "Data" object which represents all the observations from a single realization. Our likelihood will then take in a "Data" object and the associated rates.

We should write some wrapper function that takes in a "Realization" and spits out a "Data" object, ensuring the "Data" object is built correctly and all that jazz. "Data" objects might need a reference to a "Realization", and then can compute the statistics they need on the fly...

Need to do set up "Data" objects with references back to "Realization"

We then need to set up a likelihood that takes in a list of "Data" and computes things automatically using the information contained only in the Data objects (ie, the Data objects must be sufficient) and a set of model parameters: rateS, rates (list with length==numIFOs)
