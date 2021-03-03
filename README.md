# infant_lookaway_JAGS
Bayesian hierarchical model for the detection of individual differences in infants' looking behavior

The model works on the data that was collected in Poli et al. (2020): https://advances.sciencemag.org/content/6/39/eabb5053

Data is available in the repository http://hdl.handle.net/11633/aadgbtz7 or can be made available by me upon request.

The model can be extended to work on any dataset in which participants' free-looking behavior is recorded.

- lookAway.txt contains the JAGS model
- lookAway.R is the main file where you can upload the data and run the model
- samples contains an example of the output of the model
- analysis.R works on samples to return plots of the results
- analysisFunctions.R and DBDA2E-few-utilities.R are used in analysis.R
