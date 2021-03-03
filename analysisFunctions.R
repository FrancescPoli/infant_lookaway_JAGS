### Set up workspace
## Set Libraries
library(R2jags)                   # JAGS
library(ggmcmc)                   # mcmc to ggplot2
library(ggthemes)                 # ggplot2 themes
library(logspline)                # logarithmic spline fitting
library(scales)                   # logarithmic scaling
source("DBDA2E-few-utilities.R")  # DBDA2E functions

### Function definitions
## Calculating Bayes factors
# On the individual level
individual_bayes_factors = function(posterior, prior_point, h0_value, ...) {
  BF10 = c()
  
  # Calculate individual Bayes factors
  for(i in 1:length(posterior[1,])) {
    fit = logspline(posterior[,i], error.action = 2, ...)
    post_point = dlogspline(h0_value, fit)
    BF10 = c(BF10, prior_point/post_point)
  }
  
  # Return individual Bayes factors
  return(BF10)
}

# On the group level
group_bayes_factor = function(posterior, prior_point, h0_value, rows = sample(1:length(posterior[,1]), min(length(posterior[,1]),10000)), ...) {
  # Subset the posterior
  group_posterior = posterior[rows,]
  
  # Calculate group Bayes Factor
  fit = logspline(group_posterior, error.action = 2, ...)
  post_point = dlogspline(h0_value, fit)
  BF10 = prior_point/post_point
  
  # Return group Bayes factors
  return(BF10)
}

# On all levels
bayes_factors = function(posterior, prior_point, h0_value, rows = sample(1:length(posterior[,1]), min(length(posterior[,1]),10000)), ...) {
  BF10 = c(individual_bayes_factors(posterior, prior_point, h0_value))
  BF10 = c(group_bayes_factor(posterior, prior_point, h0_value, rows = rows, ...), BF10)
  
  return(BF10)
}

## Calculating HDIs
# On the individual level
individual_HDIs = function(posterior, credMass=0.95, names=c("min", "max")) {
  HDI = t(apply(posterior, 2, HDIofMCMC, credMass=credMass))
  colnames(HDI) = names
  return(HDI)
}

# On the group level
group_HDIs = function(posteriors, credMass=0.95, names=NULL) {
  HDIs = sapply(posteriors, HDIofMCMC, credMass=credMass)
  rownames(HDIs) = names
  return(HDIs)
}

## Plotting Bayes factors
# On the individual level
plot_individual_bayes_factors = function(BF10, label, limits) {
  # Create dataframe
  N = length(BF10)
  BF10.df = data.frame("Participant"=1:N, "BF10"=BF10)
  
  # Create barplot with log transformed Bayes factors
  ggplot(BF10.df, aes(x=Participant, y=BF10)) +
    scale_x_continuous(name="Participant number", breaks=seq(0,N,5)) + 
    #scale_y_log10(name=parse(text=label)) +
    #scale_y_continuous(name=parse(text=label), limits=limits, breaks=seq(limits[1], limits[2], interval), trans="log") + 
    scale_y_continuous(name=label, limits=limits, trans="log", breaks=trans_breaks("log", function(x) exp(x)), labels=trans_format("log", math_format(e^.x))) + 
    geom_bar(stat="identity") + 
    theme_few()
}

# On the group level
plot_group_bayes_factor = function(BF10, label, limits) {
  # Create dataframe
  BF10.df = data.frame("BF10"=BF10)
  
  # Create barplot with log transformed Bayes factor
  ggplot(BF10.df, aes(x="", y=BF10)) +
    labs(x="Group") + 
    scale_y_continuous(name=label, limits=limits, trans="log", breaks=trans_breaks("log", function(x) exp(x)), labels=trans_format("log", math_format(e^.x))) + 
    #scale_y_continuous(name=parse(text=label), limits=limits, trans="log", breaks=pretty_breaks()) + 
    geom_bar(stat="identity") + 
    theme_few()
}

# On all levels
plot_bayes_factors = function(BF10, label, limits, parameter) {
  # Plot group Bayes factor
  print(plot_group_bayes_factor(BF10[1], label, limits))
  ggsave(path="Images", filename=paste0(parameter, "_group_BF_barplot.png"), scale=0.5, width=2.5, height=7, dpi=300)
  
  # Plot individual Bayes factors
  print(plot_individual_bayes_factors(BF10[2:length(BF10)], label, limits))
  ggsave(path="Images", filename=paste0(parameter, "_individuals_BF_barplot.png"), width=8.75, height=3.5, dpi=300)
}

## Plotting posteriors
# On the individual level
plot_individual_posteriors = function(df, parameter, N, line, label, limits, interval) {
  # Create dataframe
  X = data.frame(Parameter=sprintf(paste0(parameter,"[%d]"), 1:N), Participant=1:N)
  
  # Plot individual posteriors
  ggs_caterpillar(df, X=X, line=line, greek=TRUE, horizontal=FALSE, sort=FALSE) +
    scale_y_continuous(name="Participant number", breaks=seq(0,N,5)) +
    scale_x_continuous(name=label, limits=limits, breaks=seq(limits[1], limits[2], interval)) + 
    theme_few()
}

# On the group level
plot_group_posterior = function(df, parameter, N, line, label, limits, interval) {
  # Create dataframe
  df$Parameter = rep("l", length(df[,1]))
  
  # Plot group posterior
  ggs_caterpillar(df, line=line, greek=TRUE, horizontal=FALSE, sort=FALSE) +
    labs(y="Group") + 
    scale_x_continuous(name=label, limits=limits, breaks=seq(limits[1], limits[2], interval)) + 
    theme_few()
}

# On all levels
plot_posteriors = function(df, parameter, N, line, label, limits, interval) {
  print(plot_group_posterior(df, parameter, N, line, label, limits, interval))
  ggsave(path="Images", filename=paste0(parameter, "_group_posterior_catterpillarplot.png"), scale=0.5, width=2.5, height=7, dpi=300)
  print(plot_individual_posteriors(df, parameter, N, line, label, limits, interval))
  ggsave(path="Images", filename=paste0(parameter, "_individuals_posterior_catterpillarplot.png"), width=8.75, height=3.5, dpi=300)
}