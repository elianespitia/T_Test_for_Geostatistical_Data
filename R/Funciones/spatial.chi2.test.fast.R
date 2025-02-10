
# Description: Return only the p-value of the spatial test

spatial.chi2.test.os.fast <- function(data.x, mu0, model.x, cov.pars.x, nugget.x, alternative){
  
  # Co-located points
  
  # Data and coordinates
  x <- data.x[,3]
  coords <- data.x[,-3]
  
  # Parameter estimation
  lik_fit <- likfit(as.geodata(data.x), lik.method="ML", ini.cov.pars=cov.pars.x, cov.model=model.x,
                    fix.nugget=F, nugget=nugget.x, messages=F)
  
  # Covariance matrix
  dist.x  <- as.matrix(dist(coords))
  Sigma.x <- geoR::cov.spatial(dist.x, cov.model=model.x, cov.pars=lik_fit$cov.pars)
  diag(Sigma.x) <- diag(Sigma.x) + lik_fit$nugget
  
  # Precision matrix
  if(model.x=="gaussian"){
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
  chi2.stat <- (mx - mu0)^2 / vx
  
  # Efective sample size
  n <- nrow(data.x)
  ess <- n*sum(diag(Sigma.x))/sum(Sigma.x)
  
  # P-value
  if(alternative == "two.sided"){
    p.value = 2*min(pchisq(chi2.stat, df=1), pchisq(chi2.stat, df=1,lower.tail=F))
  }
  else if(alternative == "greater"){
    p.value = pchisq(chi2.stat, df=1, lower.tail=F)
  }
  else if(alternative == "less"){
    p.value = pchisq(chi2.stat, df=1)
  }
  
  # OUTPUT
  return(c(p.value, ess))
}

# Examples -----------------------------------------------------------------

# # Libraries
# library(geoR)
# source("Funciones/random_coords.R")
# 
# # Parameters
# n=300
# model.x="exponential"
# cov.pars.x=c(1,50)
# nugget.x = 0
# 
# pvalues <- replicate(n=1000,{
#   
#   # Simulate Data
#   grilla <- random_coords(n=n)
#   sim_data <- grf(n=n, grid=grilla, nsim=1, cov.model=model.x, cov.pars=cov.pars, mean=0, messages=F)$data
#   data.x <- cbind(grilla,sim_data)
#   
#   # test
#   spatial.chi2.test.os.fast(data.x=data.x, mu0=0, model.x=model.x, cov.pars.x=cov.pars.x, nugget.x=nugget.x, alternative="greater")
# })
# 
# dim(pvalues)
# mean(pvalues[1,] < 0.05)
# mean(pvalues[2,])
# 
# summary(pvalues[1,])
# summary(pvalues[2,])
# 
# hist(pvalues[1,])
# hist(pvalues[2,])


# sims = 1000 y cov.pars.x=c(1,50)
# model.x="exponential" y n=50  => power=0.223 | ess=4.796742
# model.x="exponential" y n=100 => power=0.243 | ess=4.565333
# model.x="exponential" y n=300 => power=0.181 | ess=3.857043