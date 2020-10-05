%%% Author Sergios Agapiou

These files contain the code corresponding to Compartmental Model 1 in Section 2.4.1
of "Modeling of Covid-19 Pandemic in Cyprus" by Agapiou et al.

The code is based on the code of "Substantial undocumented infection
facilitates the rapid dissemination of novel coronavirus (sars-cov-2)",
in Science 368, by Li et al.
This code is simplified by not considering the metapopulation
structure. An independence sampler is used instead of Kalman filter
techniques.


To run the code one needs to run the file MCMC.m (for all incidents) or MCMC_local.m (for incidents due to local transmission only).

To plot the results use plot_figures.m. 

To produce comparable figures to the ones in the article, set the number of iterations of the independence sampler to iter=10000 (currently iter is set to 100, so that the code runs in a couple of minutes).

