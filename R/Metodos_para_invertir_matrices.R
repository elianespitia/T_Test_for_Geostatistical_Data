

# Metodos para de invertir matrices
library(geoR)

n=200
phi=50
model="gaussian"
sims=10000

# 1. Cholesky descomposition (solve function) --------------------------------

time.taken <- Sys.time()
for(i in 1:sims){

  # simulate data
  grilla <- random_coords(n)
  sim_data <- grf(n=n, grid=grilla, nsim=1,
                  cov.model=model, cov.pars=c(1, phi), nugget=0,
                  mean=0, messages=F, method="svd")$data
  dat.x <- cbind(grilla, sim_data)
  # Inverte Sigma Matrix
  dist.x  <- as.matrix(dist(grilla))
  Sigma.x <- geoR::cov.spatial(dist.x, cov.model=model, cov.pars=c(1,phi))
  ISigma.x <- solve(Sigma.x)
}
time.taken <- Sys.time() - time.taken


# 2. SVD descomposition ------------------------------------------------------

time.taken <- Sys.time()
for(i in 1:sims){
  
  # simulate data
  grilla <- random_coords(n)
  sim_data <- grf(n=n, grid=grilla, nsim=1,
                  cov.model=model, cov.pars=c(1, phi), nugget=0,
                  mean=0, messages=F, method="svd")$data
  dat.x <- cbind(grilla, sim_data)
  # Inverte Sigma Matrix
  dist.x  <- as.matrix(dist(grilla))
  Sigma.x <- geoR::cov.spatial(dist.x, cov.model=model, cov.pars=c(1,phi))
  svd_decomp <- svd(Sigma.x)
  inv_diag <- 1/svd_decomp$d
  ISigma.x2 <- svd_decomp$v %*% diag(inv_diag) %*% t(svd_decomp$u)
}
time.taken <- Sys.time() - time.taken


# 3. Eigen descomposition -------------------------------------------------

time.taken <- Sys.time()
for(i in 1:sims){
  
  # simulate data
  grilla <- random_coords(n)
  sim_data <- geoR::grf(n=n, grid=grilla, nsim=1,
                  cov.model=model, cov.pars=c(1, phi), nugget=0,
                  mean=0, messages=F, method="svd")$data
  dat.x <- cbind(grilla, sim_data)
  # Inverte Sigma Matrix
  dist.x  <- as.matrix(dist(grilla))
  Sigma.x <- geoR::cov.spatial(dist.x, cov.model=model, cov.pars=c(1,phi))
  eigen_decom <- eigen(Sigma.x, symmetric=T)
  inv_diag <- 1/eigen_decom$values
  ISigma.x2 <- eigen_decom$vectors %*% diag(inv_diag) %*% t(eigen_decom$vectors)
}
time.taken <- Sys.time() - time.taken


# Resumen -----------------------------------------------------------------

# for sims=10_000
# n=20, cov.pars=c(1,50), model="exponential"  | 1: 10.02867 secs , 2: 10.78464 secs 
# n=50, cov.pars=c(1,50), model="exponential"  | 1: 23.52541 secs , 2: 31.02812 secs 
# n=200, cov.pars=c(1,50), model="exponential" | 1: 7.808053 mins , 2: 12.39773 mins, 3: 11.12377 mins 

# n=100, cov.pars=c(1,100), model="exponential" | 1: 1.466246 mins , 2: 2.168137 mins
# n=100, cov.pars=c(1,150), model="exponential" | 1: 1.449176 mins , 2: 2.114958 mins

# n=100, cov.pars=c(1,30), model="gaussian" | 1: ERROR , 2: 2.209491 mins, 3: 1.976124 mins
# n=200, cov.pars=c(1,30), model="gaussian" | 1: ERROR , 2: 11.73166 mins, 3: 11.1823

# n=100, cov.pars=c(1,30), model="spherical" | 1: 1.41956 ,  2: 2.388227 mins, 3: 2.24925
# n=200, cov.pars=c(1,30), model="spherical" | 1: 9.248367 , 2: 14.77769 mins, 3: 11.75376


# Concluison --------------------------------------------------------------

# Para los modelso exponencial y esferico usar solve()
# Para el modelo gaussiano usar sdv() o eigen(), esta ultima es mas rapido

