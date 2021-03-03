### Set up workspace
## Clear graphics and memory
graphics.off()                      # Close all of R's graphics windows
rm(list=ls())                       # Clear all of R's memory
memory.limit(size = 10000000000000) # Enlargen memory

## Set libraries
library(R2jags)     # JAGS
library(optimbase)  # Matrices

### Preprocess datasets
## Import stimulus locations
X = read.csv("sequences_t.csv",head=F)
## Calculate max trial and sequence numbers
T = length(X[1,])
S = length(X[,1])
## Give it meaningfull names
names(X) = sprintf("t%d", 1:T)
row.names(X) = sprintf("s%d", 1:S)

## Import data
D = read.csv("LookAwayData.csv",head=F)
## Give it short (or long) meaningful names
names(D) = c("subject", "sequence", "first trial", "last trial", "look away")
names(D) = c("i", "s", "t_f", "t_l", "la")
## Exclude participants (<20 trials watched)
excluded = 0
for(i in 1:tail(D,1)$i) {
  # Take subset with sequences of that participant
  rows = D[D$i==i,]
  # Check if the participant has watched <20 trials
  if(sum(rows$t_l-(rows$t_f-1))<20) {
    # If so, delete the data of this participant
    D = D[!(D$i==i),]
    excluded = excluded + 1
  } else { # Else, Reassign participant numbers
    D$i[D$i==i] = i-excluded
  }
}; rm(i); rm(rows)
## Delete too short sequences (<4 trials watched)
D = D[D$t_l-D$t_f>2,]
## Calculate max participant number
P = tail(D,1)[[1]]

## Generate first trial number matrix
FT=matrix(-1, nrow = P, ncol = S)
# where FT[i,s] is the first trial watched
# where -1 means the sequence has been skipped
for(e in 1:length(D[,1])) {
  FT[D[e,1],D[e,2]] = D$t_f[e]
}; rm(e)

## Generate look away trial matrix
Y=matrix(16, nrow = P, ncol = S)
# where Y[i,s] is the trial of look away
for(e in 1:length(D[,1])) {
  Y[D[e,1],D[e,2]] = D$t_l[e] + abs(D$la[e]-1)
}; rm(e)

## Generate watched matrix
W=matrix(0, nrow = P*S, ncol = T)
for(e in 1:length(D[,1])) {
  i_idx=D[e,1]  # get participant number
  s_idx=D[e,2]  # get sequence number
  tf_idx=D[e,3] # get first_trial number
  # now create vector for this sequence where 0 is before first trial and 1 is from first trial on
  # last trial is not specified here but in X
  seq=c(rep(0,tf_idx-1),rep(1,T+1-tf_idx))
  # add vector to matrix W
  W[(i_idx-1)*S+s_idx,]=seq
}; rm(e, s_idx, i_idx, tf_idx, seq)

## Generate counts matrix
C=matrix(0, nrow = P*S, ncol = T*4)
for(i in 1:P) {
  for(s in 1:S) {
    for(t in 1:T) {
      for(l in 1:4) {
        if(t!=1) {
          C[(i-1)*S+s,(t-1)*4+l] = C[(i-1)*S+s,(t-2)*4+l]
        }
        if(W[(i-1)*16+s,t]==1 & l==X[s,t]) {
          C[(i-1)*S+s,(t-1)*4+l] = C[(i-1)*S+s,(t-1)*4+l] + 1
        }
      }; rm(l)
    }; rm(t)
  }; rm(s)
}; rm(i)

L=C
for(i in 1:P) {
  for(s in 1:S) {
    for(t in T:1) {
      for(l in 1:4) {
        if(t==1) {
          L[(i-1)*S+s,(t-1)*4+l] = 0
        } else {
          L[(i-1)*S+s,(t-1)*4+l] = L[(i-1)*S+s,(t-2)*4+l]
        }
      }; rm(l)
    }; rm(t)
  }; rm(s)
}; rm(i)


### Use JAGS
## Encapture data for JAGS
data <- list("FT", "C", "L", "P", "S", "T", "W", "X", "Y")

## Set starting values
myinits <- list(list(r_alpha = 2, l_alpha = 1,
                     mean_b0 = 0, prec_b0 = 0.0001,
                     mean_bI = 0, prec_bI = 0.0001,
                     mean_bD = 0, prec_bD = 0.0001))

## Set parameters to be monitored
parameters <- c("r_alpha", "l_alpha", "alpha",
                "mean_b0", "prec_b0", "b0",
                "mean_bI", "prec_bI", "bI",
                "mean_bD", "prec_bD", "bD",
                "pred_Y")

## Collect garbage and run JAGS
gc()
samples <- jags(data=data, parameters.to.save = parameters,
                model.file ="lookAway.txt",
                n.chains=8, n.iter=300000, n.burnin=150000, n.thin=100, DIC=T, progress.bar="text")
# Now the values for the monitored parameters are in the "samples" object,
# ready for inspection.

## Save the samples for later inspections
save(samples, file = "samples_c8_i300K_b150K_t100_01_06_2020")

### Model analysis
# See analysis.R