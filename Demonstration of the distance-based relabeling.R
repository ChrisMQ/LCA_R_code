# Load the required library
library(poLCA)

# Function used to create distinct permutations
permn <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permn(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}

# Function used to convert labels to latent indicators
dataPattern <- function(labels) {
  n <- length(labels)
  nc <- max(labels)
  patM <- matrix(0, nrow=n, ncol=nc)
  for(i in 1:n) {
    opt <- labels[i]
    patM[i,][opt] <- 1
  }
  return(patM)
}

# Function used to find the new order using the distance-based approach
compLab <- function(trueclass, cprobs) {
  pat <- dataPattern(trueclass)
  permus <- permn(max(trueclass))
  complik <- c()
  for(i in 1:nrow(permus)) {
    complik[i] <- sum(pat*cprobs[,permus[i,]])
  }
  opt <- which.max(complik)
  res <- permus[opt,]
  return(res)
}

# Function used to convert probs list yieleded by poLCA to a matrix of K by J
convertProb <- function(prob, J, K) {
  res <- matrix(0, nrow=K, ncol=J)
  for(j in 1:J) {
    res[,j] <- prob[[j]][,2]
  }
  rownames(res) <- paste0("Class ", 1:K)
  colnames(res) <- paste0("Y", 1:J)
  res <- round(res, 3)
  return(res)
}

#==================================
# The demonstration starts here 
#==================================

# Specify the population values
K <- 3
J <- 8
N <- 200
P <- c(0.34, 0.33, 0.33)
IRP <- matrix(c(rep(0.8,J), rep(c(0.8,0.2),c(J/2,J/2)), rep(0.2,J)),
              nrow=K, ncol=J, byrow=T)

# Simulate the true class label and observed values for the individuals
set.seed(123)
trueclass <- sample(1:K, prob = props, size = N, replace = TRUE) 
dat <- matrix(0, nrow=N, ncol=J)
for(i in 1:N) {
  dat[i,] <- rbinom(n=J, size=1, p=ip[trueclass[i],])
} 
dat <- as.data.frame(dat) + 1
colnames(dat) <- paste0('Y',1:J) 

# Fit the model
f <- as.formula(paste("cbind(",paste(paste0("Y",1:J),collapse=","),")","~1"))
mod <- poLCA(f, dat, nclass=K, maxiter=10000, nrep=20, verbose=FALSE)

# Extract the estimated IRPs
# Before relabeling, the order is 2, 3, 1
# We want to make the order to be 1, 2, 3, which is consistent with the true order
(IRP.raw <- convertProb(mod$probs, J, K))

# Identify the new order using the distance-based approach
cprobs.raw <- mod$posterior # extract the posterior classification probabilities
neworder <- compLab(trueclass, cprobs.raw) # find the new order used for relabeling
IRP.new <- IRP.raw[neworder, ] # reorder the class-specific IRPs

# Compare with the true parameter values; the two orders are consistent.
IRP.new
IRP
