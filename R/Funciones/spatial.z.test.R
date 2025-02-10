

#rm(list=ls())

# Librerias
library(geoR)

# # Para probar valores
# n=20
# model = "exponential"
# mu0 = 0
# ha = "two.sided"
# set.seed(123)
# coord.x = random_coords(n=n)
# coord.y = random_coords(n=n)
# dat.x <- grf(n=n, grid=coord.x, nsim=1, cov.model="exponential", cov.pars=c(1, 50), nugget=0,
#         mean=0, messages=F)$data
# dat.y <- grf(n=n, grid=coord.x, nsim=1, cov.model="exponential", cov.pars=c(1, 50), nugget=0,
#          mean=0, messages=F)$data
# dat.x = data.frame(coord.x, dat.x)
# dat.y = data.frame(coord.y, dat.y)
# ?t.test


# Funcion -----------------------------------------------------------------
spatial.z.test <- function(dat.x, dat.y=NULL, model.x, model.y, cov.pars.x=NULL, cov.pars.y=NULL, mu0=0, alternative="two.sided", alpha=0.05){
  
  # Verificar inputs
  if(is.null(dat.x)){stop("x no puede ser nulo")}
  
  
  # Crear objeto geodata para dat.x
  x <- dat.x[,3]
  nx <- nrow(dat.x)
  names(dat.x) <- c("x", "y", "data")
  dat.x <- as.geodata(dat.x, coords.col=c(1,2), data.col=3) 

  # Caso para 1 muestra ---------------------------------------------------

  # Estimar Sigma: Usar el modelo que ya se sabe
  # cov_model_fit <- geoR::likfit(geo.x, cov.model = model.x, ini.cov.pars=c(5,20), trend="cte", messages=F)
  # dist_matrix   <- as.matrix(dist(coord.x))
  # Sigma <- geoR::cov.spatial(dist_matrix, cov.pars = cov_model_fit$cov.pars, cov.model = model.x)
  
  # Matriz de distancias, de covarianza, y de Precision para X
  dist.x  <- as.matrix(dist(dat.x$coords))
  Sigma.x <- geoR::cov.spatial(dist.x, cov.model=model.x, cov.pars=cov.pars.x)
  ISigma.x <- solve(Sigma.x)
  
  # Estimador de la media y su varianza
  mx <- as.numeric((colSums(ISigma.x) %*% x) / sum(ISigma.x))
  vx <- 1/sum(ISigma.x)
  
  # Estadistica Z
  if(is.null(dat.y)){
    # Estadistica 1 muestra
    Z.stat <- (mx - mu0) / sqrt(vx)
  } else{
    # Crear objeto geodata para Y
    y <- dat.y[,3]
    ny <- nrow(dat.y)
    names(dat.y) <- c("x", "y", "data")
    dat.y <- as.geodata(dat.y, coords.col=c(1,2), data.col=3)
    # Matriz de distancias y de covarianza para Y
    dist.y   <- as.matrix(dist(dat.y$coords))
    Sigma.y  <- geoR::cov.spatial(dist.y, cov.model=model.y, cov.pars=cov.pars.y)
    ISigma.y <- solve(Sigma.y)
    # Estimador de la media y su varianza
    my <- as.numeric((colSums(ISigma.y) %*% x) / sum(ISigma.y))
    vy <- 1/sum(ISigma.y)
    # Estadistica 2 muestras
    Z.stat <- (mx-my - mu0) / (sqrt(vx + vy))
  }
  
  # Pvalor
  if(alternative == "two.sided"){
    p.value = 2*min(pnorm(Z.stat), pnorm(Z.stat, lower.tail = F))
    }
  else if(alternative == "greater"){
    p.value = pnorm(Z.stat, lower.tail = F)
    }
  else if(alternative == "less"){
    p.value = pnorm(Z.stat)
  }
  
  # Output
  if(is.null(dat.y)){
    output = list("mx"=mx,
                  "vx"=vx,
                  "Z.stat"=Z.stat, 
                  "p.value"=p.value, 
                  "mu0"=mu0, 
                  "alternative"=alternative)
  } else{
    output = list("mx"=mx,
                  "my"=my,
                  "vx"=vx,
                  "vy"=vy,
                  "Z.stat"=Z.stat, 
                  "p.value"=p.value, 
                  "mu0"=mu0, 
                  "alternative"=alternative)
  }
  
  # Output
  return(output)
  
} # Fin de la funcion

# Prueba ------------------------------------------------------------------

# Prueba 1 muestra
# n=20
# set.seed(123)
# coord.x = random_coords(n=n)
# x <- grf(n=n, grid=coord.x, nsim=1, cov.model="exponential", cov.pars=c(1, 20), nugget=0,
#          mean=0, messages=F)
# dat.x <- as.data.frame(x)
# spatial.z.test(dat.x, mu0=0, model.x="exponential", cov.pars.x=c(1,50), alternative="two.sided", alpha=0.05)








