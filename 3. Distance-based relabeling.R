#========================================================================
# Using Distance-based Approach to Address Label Switching in Simulation 
#========================================================================

#=========== Load the required library ===========

library(poLCA)

#=========== Load user-defined functions ===========

source("R/User-defined functions.R")

#=========== The demonstration starts here ===========

# Specify the population values
K <- 3   # number of classes in population
J <- 8   # number of indicators
N <- 200 # sample size
P <- c(0.34, 0.33, 0.33) # mixing proportions of the three classes 
IRP <- matrix(c(rep(0.8,J), rep(c(0.8,0.2), c(J/2,J/2)), rep(0.2,J)), # item response probabilities for each of the classes
              nrow = K, ncol = J, byrow = TRUE)

# Simulate true class memberships and observed values for the individuals
set.seed(123)
trueclass <- sample(1:K, prob = P, size = N, replace = TRUE) 
dat <- matrix(0, nrow = N, ncol = J)
for(i in 1:N) {
  dat[i,] <- rbinom(n = J, size = 1, p = IRP[trueclass[i],])
} 
dat <- as.data.frame(dat) + 1
colnames(dat) <- paste0('Y', 1:J) 

# Fit the model
f <- as.formula(paste("cbind(",paste(paste0("Y",1:J),collapse=","),")","~1")) # specify model formula
mod <- poLCA(f, dat, nclass = K, maxiter = 10000, nrep = 20, verbose = FALSE) # fit the specified model

# Extract the estimated IRPs
# Before relabeling, the order is 2, 3, 1
# We want to make the order to be 1, 2, 3, which is consistent with the true order
(IRP.raw <- convertProb(mod$probs, J, K))

# Identify the new order using the distance-based approach
cprobs.raw <- mod$posterior # extract the posterior classification probabilities
neworder <- compLab(trueclass, cprobs.raw) # find the new order used for relabeling
IRP.new <- IRP.raw[neworder, ] # reorder the class-specific IRP estimates

# Compare the order of the parameter estimates with that of the true parameters 
# The two orders are consistent, indicating that label switching is solved
IRP.new # relabeled IRP estimates
IRP     # true parameters
