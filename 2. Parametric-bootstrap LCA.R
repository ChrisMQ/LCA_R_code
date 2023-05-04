#======================================================================
# Parametric Bootstrap LCA (PB-LCA) Using the Distance-based Approach 
#======================================================================

#=========== Load the required libraries ===========

library(poLCA)
library(dplyr)

#=========== Load user-defined functions ===========

source("R/User-defined functions.R")

#=========== Demonstration using a simulated data ===========

# Specify the population values
K <- 3   # number of classes in population
J <- 8   # number of indicators
N <- 500 # sample size
P <- c(0.2, 0.3, 0.5) # mixing proportions of the three classes
IRP <- matrix(c(rep(0.9,J), rep(c(0.9,0.1), c(J/2,J/2)), rep(0.1,J)), # item response probabilities for each of the classes
              nrow = K, ncol = J, byrow = TRUE)

# Simulate true class memberships and observed values for the individuals
set.seed(123)
trueclass <- sample(1:K, prob = P, size = N, replace = TRUE) # simulate true class membership for each individual
dat <- matrix(0, nrow = N, ncol = J)
for(i in 1:N) {
  dat[i,] <- rbinom(n = J, size = 1, p = IRP[trueclass[i],])
} 
dat <- as.data.frame(dat) + 1    # the simulated values were added by 1 because the poLCA function only allows for binary data of 1 and 2                             
colnames(dat) <- paste0('Y', 1:J) # assign variable names

# Run the PB-LCA
set.seed(123)
mod <- lca_pb(dat, nclass = K, verbose = T)

# Check parameter estimates
mod$P.est  # mixing proportion estimates
mod$ip.est # item response probability estimates
