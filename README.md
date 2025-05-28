# Files
1) swine_model_asf_v9_anonymized.nlogo - NetLogo simulation model of a potential ASF spread in the Netherlands with several biosecurity and control options that can be activated.
2) swine_model_asf_v9_anonymized.html - NetLogo Web version of the ABM-ASF model. This can be run in a web browser without installing the NetLogo desktop version.
3) swine_asf_model_lhs.R - R code for generating the prior distribution to be used in the parameter estimation using the Latin Hypercube Sampling (LHS).
4) sim_asf_model_nlrx.R - R code for simulating the NetLogo model headless in R without activating biosecurity and control measures. The user has the option to indicate the number of simulations and seeds.
5) sim_asf_model_nlrx_biosecurity.R - R code for simulating the NetLogo model headless in R with biosecurity and control measures. The user has the option to indicate the number of simulations and seeds. The user can also indicate the combinations of parameters to use for scenario analysis.
6) param_est_ABC.Rmd - R notebook that performs the Approximate Bayesian Computation methods to estimate patch-to-patch transmission parameters using rejection, loclinear, and neuralnet methods.
