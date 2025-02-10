
# Description: Return only the p-value of the spatial test


# 1. One sample -----------------------------------------------------------
spatial.z.test.os.fast <- function(data.x, mu0, model.x, cov.pars.x, nugget.x, alternative){
  
  # Co-located points
  
  # Data and coordinates
  x <- data.x[,3]
  coords <- data.x[,-3]
    
  # Covariance matrix
  dist.x  <- as.matrix(dist(coords))
  Sigma.x <- geoR::cov.spatial(dist.x, cov.model=model.x, cov.pars=cov.pars.x)
  diag(Sigma.x) <- diag(Sigma.x) + nugget.x
  
  # Precision matrix
  if(model.x=="gaussian"){
    # Eigen descomposition
    # eigen_decom <- eigen(Sigma.x, symmetric=T)
    # inv_diag <- 1/eigen_decom$values
    # ISigma.x <- eigen_decom$vectors %*% diag(inv_diag) %*% t(eigen_decom$vectors)
    
    # SVD descomposition
    svd_decomp <- svd(Sigma.x)
    inv_diag <- 1/svd_decomp$d
    ISigma.x <- svd_decomp$v %*% diag(inv_diag) %*% t(svd_decomp$u)
  } else {
    ISigma.x <- solve(Sigma.x)
  }
  
  # Statistic
  mx <- as.numeric((colSums(ISigma.x) %*% x)/sum(ISigma.x))
  vx <- 1/sum(ISigma.x)
  z.stat <- (mx - mu0) / sqrt(vx)
  
  # P-value
  if(alternative == "two.sided"){
    p.value = 2*min(pnorm(z.stat), pnorm(z.stat, lower.tail=F))
  }
  else if(alternative == "greater"){
    p.value = pnorm(z.stat, lower.tail=F)
  }
  else if(alternative == "less"){
    p.value = pnorm(z.stat)
  }
  
  # OUTPUT
  return(p.value)

} # Close function


# 2. Two sample --------------------------------------------------------------
spatial.z.test.ts.fast <- function(data.x, data.y, mu0, model.x, model.y, cov.pars.x, cov.pars.y, nugget.x, nugget.y, alternative){
  
  # Co-located points
  
  # Data and coordinates
  x <- data.x[,3]
  y <- data.y[,3]
  coords.x <- data.x[,-3]
  coords.y <- data.y[,-3]
    
  # Covariance matrix
  dist.x  <- as.matrix(dist(coords.x))
  dist.y  <- as.matrix(dist(coords.y))
  Sigma.x <- geoR::cov.spatial(dist.x, cov.model=model.x, cov.pars=cov.pars.x)
  Sigma.y <- geoR::cov.spatial(dist.y, cov.model=model.y, cov.pars=cov.pars.y)
  diag(Sigma.x) <- diag(Sigma.x) + nugget.x
  diag(Sigma.y) <- diag(Sigma.y) + nugget.y
  
  # Precision matrix
  if(model.x=="gaussian"){
    # SVD descomposition
    svd_decomp <- svd(Sigma.x)
    inv_diag <- 1/svd_decomp$d
    ISigma.x <- svd_decomp$v %*% diag(inv_diag) %*% t(svd_decomp$u)
  } else {
    ISigma.x <- solve(Sigma.x)
  }
  if(model.y=="gaussian"){
    # SVD descomposition
    svd_decomp <- svd(Sigma.y)
    inv_diag <- 1/svd_decomp$d
    ISigma.y <- svd_decomp$v %*% diag(inv_diag) %*% t(svd_decomp$u)
  } else {
    ISigma.y <- solve(Sigma.y)
  }
  
  # Statistic
  mx <- as.numeric((colSums(ISigma.x) %*% x) / sum(ISigma.x))
  my <- as.numeric((colSums(ISigma.y) %*% y) / sum(ISigma.y))
  vx <- 1/sum(ISigma.x)
  vy <- 1/sum(ISigma.y)
  z.stat <- (mx-my - mu0) / sqrt(vx + vy)
  
  # P-value
  if(alternative == "two.sided"){
    p.value = 2*min(pnorm(z.stat), pnorm(z.stat, lower.tail=F))
  }
  else if(alternative == "greater"){
    p.value = pnorm(z.stat, lower.tail=F)
  }
  else if(alternative == "less"){
    p.value = pnorm(z.stat)
  }
  
  # OUTPUT
  return(p.value)
  
} # Close function


# Examples -----------------------------------------------------------------

# # libraries
# library(geoR)
# source("random_coords.R")
# # parameters
# n=200
# model="gaussian"
# cov.pars=c(1,120)
# # data
# data.x <- as.data.frame(grf(n=n, grid=random_coords(n=n), nsim=1,
#                            cov.model=model, cov.pars=cov.pars, mean=0, messages=T))
# data.x <- as.matrix(data.x)
# data.y <- as.data.frame(grf(n=n, grid=random_coords(n=n), nsim=1,
#                            cov.model=model, cov.pars=cov.pars, mean=0, messages=F))
# data.y <- as.matrix(data.y)
# # Examples
# spatial.z.test.os.fast(data.x=data.x, mu0=0, model.x=model, cov.pars.x=cov.pars, alternative="greater")
# spatial.z.test.ts.fast(data.x=data.x, data.y=data.y, mu0=0, model.x=model, model.y=model, cov.pars.x=cov.pars, cov.pars.y=cov.pars, alternative="greater")
