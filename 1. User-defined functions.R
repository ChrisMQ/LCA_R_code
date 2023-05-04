#==========================
# User-defined Function 
#==========================

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

# Function used to compute ML CI
ml.ci <- function(theta.est, se.est, conf.level=.95) {
  alpha <- 1 - conf.level
  cv <- qnorm(p=alpha/2, mean=0, sd=1, lower.tail=F)
  ci <- c(theta.est - cv*se.est, theta.est + cv*se.est)
  if(ci[1] < 0) ci[1] <- 0
  if(ci[2] > 1) ci[2] <- 1
  labs <- c(alpha/2, 1-alpha/2)*100
  labs <- paste0(labs, "%")
  names(ci) <- labs
  return(ci)
}

# Function used to compute Student's t CI
st.ci <- function(theta.est, se.bt, N, conf.level=.95) {
  alpha <- 1 - conf.level
  tn <- qt(p=alpha/2, df=N-1, lower.tail=F)
  tp <- tn*se.bt
  low <- theta.est - tp
  up <- theta.est + tp
  ci <- c(low, up)
  if(ci[1] < 0) ci[1] <- 0
  if(ci[2] > 1) ci[2] <- 1
  labs <- c(alpha/2, 1-alpha/2)*100
  labs <- paste0(labs, "%")
  names(ci) <- labs
  return(ci)
}

# Function used to compute percentile CI
pt.ci <- function(theta.bt, conf.level=.95) {
  alpha <- 1 - conf.level
  low <- alpha/2
  up <- 1 - low
  ci <- quantile(theta.bt, probs = c(low, up))
  if(ci[1] < 0) ci[1] <- 0
  if(ci[2] > 1) ci[2] <- 1
  return(ci)
}

# Function used to compute bias-corrected CI 
bc.ci <- function(theta.est, theta.bt, conf.level=.95){
  low <- (1 - conf.level)/2
  up <- 1 - low
  sims <- length(theta.bt)
  z.inv <- length(theta.bt[theta.bt < theta.est])/sims
  z <- qnorm(z.inv)
  lower.inv <-  pnorm(z + z + qnorm(low))
  lower <- quantile(theta.bt, lower.inv, names=FALSE)
  upper.inv <-  pnorm(z + z + qnorm(up))
  upper <- quantile(theta.bt, upper.inv, names=FALSE)
  ci <- c(lower, upper)
  labs <- c(low, up)*100
  labs <- paste0(labs, "%")
  names(ci) <- labs
  if(ci[1] < 0) ci[1] <- 0
  if(ci[2] > 1) ci[2] <- 1
  return(ci)
}

# Function used to perform parametric bootstrap (PB) LCA based on the poLCA function
# dat: observations
# nclass: number of classes pre-specified
# alpha: significance level
# maxiter: maximum number of iterations through which the estimation algorithm will cycle
# nrep: number of times to estimate the model
# B: number of bootstrap replications
# verbose: whether the function should print the process bar

lca_pb <- function(dat, nclass, alpha=0.05, maxiter=10000, nrep=20, B=100, verbose=FALSE) {
  N <- nrow(dat)
  nvar <- ncol(dat)
  # model fit using the original sample
  f <- as.formula(paste("cbind(",paste(paste0("Y",1:nvar),collapse="," ),")","~1"))
  mod <- poLCA(f, dat, nclass, maxiter=maxiter, nrep=nrep, verbose=F, graphs=F)
  P0 <- mod$P
  P0.se <- mod$P.se
  ip0 <- convertProb(mod$probs, nvar, nclass)
  ip0.se <- convertProb(mod$probs.se, nvar, nclass)
  predclass0 <- mod$predclass
  posterior0 <- mod$posterior
  
  # containers of bootstrap parameter estimates
  P.pb <- matrix(0, nrow=B, ncol=nclass)
  ip.pb <- array(0, c(nclass, nvar, B))
  
  # parametric bootstrap starts here
  if(verbose)
    bar <- txtProgressBar(min=0, max=B, style=3, width=40, char="=")
  for(b in 1:B){
    truec.b <- sample(1:nclass, prob=P0, size=N, replace=TRUE) 
    dat.b <- matrix(0, nrow=N, ncol=nvar)
    for(i in 1:N) {
      dat.b[i,] <- rbinom(n=nvar, size=1, p=ip0[truec.b[i],])
    } 
    dat.b <- as.data.frame(dat.b) + 1
    colnames(dat.b) <- paste0('Y',1:nvar)
    pb.mod <- poLCA(f, dat.b, nclass, maxiter=maxiter, nrep=nrep, verbose=F, graphs=F)
    neworder.b <- compLab(truec.b, pb.mod$posterior)
    P.b.new <- pb.mod$P[neworder.b]
    ip.b <- convertProb(pb.mod$probs, nvar, nclass)
    ip.b.new <- ip.b[neworder.b, ]
    P.pb[b,] <- P.b.new
    ip.pb[,,b] <- ip.b.new
    if(verbose) 
      setTxtProgressBar(bar, b)
  }
  if(verbose) close(bar) 
  
  # convert the parameter estimates of each bootstrap replicate to vector format
  pb.paras <- matrix(0, nrow=B, ncol=(nclass+nvar*nclass))
  for(b in 1:B) {
    P <- P.pb[b,]
    ip <- c(t(ip.pb[,,b]))
    pb.paras[b, ] <- c(P, ip)
  }
  
  # create names for each parameter
  P.names <- paste0("P",1:nclass)
  K.ind <- rep(1:nclass, each=nvar)
  J.ind <- rep(1:nvar, nclass)
  ip.names <- c()
  for(i in 1:(nvar*nclass)) {
    tt <- paste0("ip.", J.ind[i], K.ind[i])
    ip.names <- c(ip.names, tt)
  }
  para.names <- c(P.names, ip.names)
  colnames(pb.paras) <- para.names 
  
  # confidence level
  conf.level <- 1-alpha 
  
  # compute required quantities
  ML.est <- c(P0, c(t(ip0)))
  ML.se <- c(P0.se, c(t(ip0.se)))
  pb.est <- colMeans(pb.paras)
  pb.se <- apply(pb.paras, 2, sd)
  
  # compute all CIs
  ci.ML <- t(mapply(ml.ci, ML.est, ML.se, conf.level))
  ci.st <- t(mapply(st.ci, pb.est, pb.se, N, conf.level))
  ci.pt <- t(apply(pb.paras, 2, pt.ci, conf.level))
  ci.bc <- matrix(0, nrow=(nclass+nvar*nclass), ncol=2)
  labs <- c(alpha/2, 1-alpha/2)*100
  labs <- paste0(labs, "%")
  colnames(ci.bc) <- labs
  for(d in 1:(nclass+nvar*nclass)) {
    ci.bc[d, ] <- bc.ci(ML.est[d], pb.paras[,d], conf.level=conf.level)
  }
  rownames(ci.ML) <- rownames(ci.st) <- rownames(ci.pt) <- rownames(ci.bc) <- para.names
  
  # put all CIs together and round the values to 4 decimal places
  all.CI <- list(ML=ci.ML, Student_t=ci.st, Percentile=ci.pt, Bias_corrected=ci.bc) 
  all.CI <- lapply(all.CI, round, 4)
  
  # aggregate the results
  # P.est: parameter estimates of mixing proportions produced by ML
  # ip.est: parameter estimates of item probabilities produced by ML
  # all.CI: four types of confidence intervals 
  # pb.est: parametric bootstrap estimates for each bootstrap replicate
  # predclass: class membership for each individual
  # posterior: posterior classification probabilities for each individual
  res <- list(P.est=P0, ip.est=ip0, 
              all.CI=all.CI,
              pb.est=pb.paras,
              predclass=predclass0, 
              posterior=posterior0)
  return(res) 
}
