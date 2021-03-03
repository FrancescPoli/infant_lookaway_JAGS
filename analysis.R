### Set up workspace
## Clear graphics and memory
graphics.off() # Close all of R's graphics windows
rm(list=ls())  # Clear all of R's memory
memory.limit(size = 100000000000) # Extend R's memory

## Set Libraries and sources
library(R2jags)               # JAGS
library(ggmcmc)               # mcmc to ggplot2
library(ggthemes)             # ggplot2 themes
library(logspline)            # logarithmic spline fitting
library(scales)                   # logarithmic scaling
source("analysisFunctions.R") # wrapper functions

### Load samples
load("samples")
samples.mcmc <- as.mcmc(samples)
N=length(samples$BUGSoutput$mean$alpha)
S=samples$BUGSoutput$n.sims

### Preprocessing
## Retrieve posteriors
a.posterior = samples$BUGSoutput$sims.list$alpha
b0.posterior = samples$BUGSoutput$sims.list$b0
bI.posterior = samples$BUGSoutput$sims.list$bI
bD.posterior = samples$BUGSoutput$sims.list$bD

## Calculate prior points
a.prior.point = dgamma(1,1,2)
b0.prior.point = dnorm(0, dnorm(0,0.0001), 1/dgamma(1,1))
bI.prior.point = dnorm(0, dnorm(0,0.0001), 1/dgamma(1,1))
bD.prior.point = dnorm(0, dnorm(0,0.0001), 1/dgamma(1,1))

## Make tidy objects using ggs()
a_ggs <- ggs(samples.mcmc, family="^alpha")
b0_ggs <- ggs(samples.mcmc, family="^b0")
bI_ggs <- ggs(samples.mcmc, family="^bI")
bD_ggs <- ggs(samples.mcmc, family="^bD")

### Perform calculations
## Calculate Bayes factors
a.BF10 = bayes_factors(a.posterior, a.prior.point, 0, rows=1:1200, lbound=0)
b0.BF10 = bayes_factors(b0.posterior, b0.prior.point, 0, rows=1:1200)
bI.BF10 = bayes_factors(bI.posterior, bI.prior.point, 0, rows=1:1200)
bD.BF10 = bayes_factors(bD.posterior, bD.prior.point, 0, rows=1:1200)

## Calculate HDIs
# on the individual level
a.HDIs = individual_HDIs(a.posterior, names=c("a.min", "a.max"))
b0.HDIs = individual_HDIs(b0.posterior, names=c("b0.min", "b0.max"))
bI.HDIs = individual_HDIs(bI.posterior, names=c("bI.min", "bI.max"))
bD.HDIs = individual_HDIs(bD.posterior, names=c("bD.min", "bD.max"))

# on the group level
group.HDIs = group_HDIs(list(a.posterior, b0.posterior, bI.posterior, bD.posterior))

### Create matrix
## on the individual level
individual.results = cbind("n"=1:N,
                           a.HDIs, "a.BF10"=a.BF10[2:71],
                           b0.HDIs, "b0.BF10"=b0.BF10[2:71],
                           bI.HDIs, "bI.BF10"=bI.BF10[2:71],
                           bD.HDIs, "bD.BF10"=bD.BF10[2:71])

## on the group level
group.results = rbind(group.HDIs, c(a.BF10[1], b0.BF10[1], bI.BF10[1], bD.BF10[1]))
colnames(group.results) = c("a", "b0", "bI", "bD")

### Make plots
## of Bayes factors
plot_bayes_factors(a.BF10, label=expression(paste("Bayes factor: ", alpha)), limits=c(1,150), parameter="alpha")
plot_bayes_factors(b0.BF10, label=expression(paste("Bayes factor: ", beta^0)), limits=c(1,200000), parameter="b0")
plot_bayes_factors(bI.BF10, label=expression(paste("Bayes factor: ", beta^I)), limits=c(0.1,1), parameter="bI")
plot_bayes_factors(bD.BF10, label=expression(paste("Bayes factor: ", beta^D)), limits=c(1,25), parameter="bD")

## of posteriors
plot_posteriors(a_ggs, "alpha", N, 0, label=expression(alpha), limits=c(0,10), interval=1)
plot_posteriors(b0_ggs, "b0", N, 0, label=expression(beta^0), limits=c(-3,0), interval=0.5)
plot_posteriors(bI_ggs, "bI", N, 0, label=expression(beta^I), limits=c(-1,1.5), interval=0.5)
plot_posteriors(bD_ggs, "bD", N, 0, label=expression(beta^D), limits=c(-34,2), interval=4)