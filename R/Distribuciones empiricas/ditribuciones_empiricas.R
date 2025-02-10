


# Librerias ---------------------------------------------------------------
library(geoR)
source("../Funciones/random_coords.R") # ../ sube un nivel



# Distribucion empirica -------------------------------------------------

# Parametros
nseq=c(10,30,50,100,200,350)
cov.model="gaussian"
short_name_model = c("exp","gauss","sphe"); names(short_name_model) <- c("exponential","gaussian","spherical")
label_model = short_name_model[cov.model]
cov.pars=c(1,30)
nugget=0
fix.nugget=T
sims=1000
# Simulation
for(n in nseq){
  # Almacenamiento
  SIGMA  <- vector(mode="numeric", length=sims)
  PHI    <- vector(mode="numeric", length=sims)
  STAT   <- matrix(NA, nrow=sims, ncol=3); colnames(STAT) <- paste0("stat_",c("gls","griffith","griffith_ess"))
  ESS    <- matrix(NA, nrow=sims, ncol=2); colnames(ESS) <- c("griffith","vallejos")
  PVALUE <- list(matrix(NA, nrow=sims, ncol=5, dimnames=list(c(1:sims), c("normal","griffith","griffith_1","vallejos","vallejos_1"))),
                 matrix(NA, nrow=sims, ncol=5, dimnames=list(c(1:sims), c("normal","griffith","griffith_1","vallejos","vallejos_1"))),
                 matrix(NA, nrow=sims, ncol=5, dimnames=list(c(1:sims), c("normal","griffith","griffith_1","vallejos","vallejos_1")))
  )
  names(PVALUE) <- colnames(STAT)
  # Simulations
  for(i in 1:sims){
    #::::::::::::::::::::::::::::::::::::::::
    # Simular campo aleatorio
    #::::::::::::::::::::::::::::::::::::::::
    # Get random coords
    grilla <- random_coords(n)
    # Simulate geostatistical data
    sim_data <- geoR::grf(n=n, grid=grilla, nsim=1,cov.model=cov.model,
                          cov.pars=cov.pars, nugget=nugget,
                          mean=0, messages=F, method="eigen")$data
    data.x <- cbind(grilla, sim_data)
    data.x <- geoR::as.geodata(data.x)
    #::::::::::::::::::::::::::::::::::::::::
    # Estimar Matriz de Covarianzas
    #::::::::::::::::::::::::::::::::::::::::
    # 1. Estimar el variograma empirico
    variog_empirical <- geoR::variog(data.x, uvec=10, messages=F) # posiblemente modificar
    # 2. Ajustar un modelo teorico
    variog_fit <- geoR::variofit(variog_empirical, cov.model=cov.model, fix.nugget=fix.nugget, nugget=nugget, messages=F)
    # 3. Contruir matriz de covarianzas
    dist.x  <- as.matrix(dist(data.x$coords))
    Sigma.x <- geoR::cov.spatial(dist.x, cov.model=cov.model, cov.pars=variog_fit$cov.pars)
    diag(Sigma.x) <- diag(Sigma.x) + variog_fit$nugget
    sg <- variog_fit$nugget + variog_fit$cov.pars[1] # varianza
    # 4. Invertir matriz de covarianzas
    if(cov.model=="gaussian"){
      # Eigen descomposition
      eigen_decom <- eigen(Sigma.x, symmetric=T)
      inv_diag <- 1/eigen_decom$values
      ISigma.x <- eigen_decom$vectors %*% diag(inv_diag) %*% t(eigen_decom$vectors)
    } else {
      ISigma.x <- solve(Sigma.x)
    }
    # Estimador GLS y su varianza
    mx <- as.numeric((colSums(ISigma.x) %*% sim_data) / sum(ISigma.x))
    vx <- 1/sum(ISigma.x)
    # Efective sample size ESS
    ess1 <- n*sum(diag(Sigma.x))/sum(Sigma.x) # Griffith 2005
    ess2 <- sum(sg*ISigma.x)                  # Vallejos 2014
    # Statistic
    stat1 <- (mx - 0) / sqrt(vx)                            # GLS estimator
    stat2 <- (mean(sim_data)- 0)/ sqrt(sum(Sigma.x)/n^2)    # Griffith 2005
    stat3 <- (mean(sim_data)- 0)/ sqrt(sum(Sigma.x)/ess1^2) # Griffith 2005 with n=n*
    stats <- c(stat1,stat2,stat3)
    # P-value
    for(j in 1:length(stats)){
      pv <- c(p.value1 = pnorm(stats[j], lower.tail=F),         # Distribucion normal
              p.value2 = pt(stats[j], df=ess1, lower.tail=F),   # Distribucion t con df = griffith
              p.value3 = pt(stats[j], df=ess1-1, lower.tail=F), # Distribucion t con df = griffith - 1
              p.value4 = pt(stats[j], df=ess2, lower.tail=F),   # Distribucion t con df = vallejos
              p.value5 = pt(stats[j], df=ess2-1, lower.tail=F)  # Distribucion t con df = vallejos - 1
      )
      PVALUE[[j]][i,] <- pv
    }
    # Almacenar resultados
    SIGMA[i]  <- variog_fit$cov.pars[1]
    PHI[i] <- variog_fit$cov.pars[2]
    STAT[i,] <- stats
    ESS[i,] <- c(ess1, ess2)
  }
  # save samples
  save_samples <- list("SIGMA"=SIGMA, "PHI"=PHI, "STAT"=STAT, "ESS"=ESS, "PVALUE"=PVALUE)
  saveRDS(save_samples, file=paste0("distri_",label_model,"_n",n,".rds")) # los guarda en la direccion donde este este script dentro del R project
}


# Load data
distri <- readRDS("distri_gauss_n350.rds")
# Resultados SIGMA
summary(distri$SIGMA)
sum(distri$SIGMA < 20 & distri$SIGMA > 10)
distri$SIGMA[which(distri$SIGMA < 20 & distri$SIGMA > 10)]

quantile(distri$SIGMA, seq(0.75,0.95,by=0.05))
mean(distri$SIGMA)
sd(distri$SIGMA)
make_hist(distri$SIGMA)
hist(distri$SIGMA, freq=T)
hist(distri$SIGMA[distri$SIGMA < 50], freq=T)
#Resultados PHI
summary(distri$PHI)
hist(distri$PHI, freq = F, ylim=c(0,0.06))
lines(density(distri$PHI))
quantile(distri$PHI, seq(0.75,0.95,by=0.05))
mean(distri$PHI)
sd(distri$PHI)
make_hist(distri$PHI, symbol="phi")
# Promedio de las stats
colMeans(STAT)
apply(distri$STAT, 2, summary)

hist(STAT[,1]); 
hist(STAT[,2])
hist(STAT[,3], freq=F)
make_hist(distri$STAT[,1])
make_hist(distri$STAT[,2])
make_hist(distri$STAT[,3])

# Encontrar los gl a partir de la distribucion empirica
q95 <- quantile(STAT[,1], 0.95)
dft <- function(x){return(abs(qt(0.95, df=x)-q95))}
optim(par=colMeans(ESS)[1], fn=dft)
qt(0.95, df=9)
# Promedio de los ess
colMeans(distri$ESS)
hist(distri$ESS[,1])
hist(distri$ESS[,2])
qt(0.95, df=colMeans(ESS)[1]-2)

# Promedio pvalores
apply(distri$PVALUE$stat_gls,2,function(x){mean(x<0.05)})
apply(distri$PVALUE$stat_griffith,2,function(x){mean(x<0.05)})
apply(distri$PVALUE$stat_griffith_ess,2,function(x){mean(x<0.05)})

plot(grilla)
plot(variog_empirical)
lines(variog_fit)
curve(1 - cov.spatial(x,cov.model = cov.model, cov.pars = cov.pars), add=T)
variog_fit$cov.pars

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Estimacion ML y REML ----------------------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Parametros
# nseq=c(20,30,50,100,200,350)
nseq=50
short_name_model = c("exp","gauss","sphe"); names(short_name_model) <- c("exponential","gaussian","spherical")
cov.model="exponential"
label_model = short_name_model[cov.model]
cov.pars=c(1,100)
nugget=0
fix.nugget=F
sims=1000

# Simulation
set.seed(123)
for(n in nseq){
  # Almacenamiento
  SIGMA  <- vector(mode="numeric", length=sims)
  PHI    <- vector(mode="numeric", length=sims)
  STAT   <- matrix(NA, nrow=sims, ncol=3); colnames(STAT) <- paste0("stat_",c("gls","griffith","griffith_ess"))
  ESS    <- matrix(NA, nrow=sims, ncol=2); colnames(ESS) <- c("griffith","vallejos")
  PVALUE <- list(matrix(NA, nrow=sims, ncol=5, dimnames=list(c(1:sims), c("normal","griffith","griffith_1","vallejos","vallejos_1"))),
                 matrix(NA, nrow=sims, ncol=5, dimnames=list(c(1:sims), c("normal","griffith","griffith_1","vallejos","vallejos_1"))),
                 matrix(NA, nrow=sims, ncol=5, dimnames=list(c(1:sims), c("normal","griffith","griffith_1","vallejos","vallejos_1")))
  )
  names(PVALUE) <- colnames(STAT)
  # Simulations
  for(i in 1:sims){
    #::::::::::::::::::::::::::::::::::::::::
    # Simular campo aleatorio
    #::::::::::::::::::::::::::::::::::::::::
    # Get random coords
    grilla <- random_coords(n=n)
    # Simulate geostatistical data
    sim_data <- geoR::grf(grid=grilla, nsim=1,cov.model=cov.model,
                          cov.pars=cov.pars, nugget=nugget,
                          mean=0, messages=F, method="eigen")$data
    data.x <- cbind(grilla, sim_data)
    data.geo <- geoR::as.geodata(data.x)
    #::::::::::::::::::::::::::::::::::::::::
    # Estimar Matriz de Covarianzas
    #::::::::::::::::::::::::::::::::::::::::
    # 1. Estimar el variograma empirico
    variog_empirical <- geoR::variog(data.geo, messages=F) # posiblemente modificar
    # 2. Ajustar modelo teorico
    variog_fit <- geoR::variofit(variog_empirical, weights="equal", cov.model=cov.model, fix.nugget=F, nugget=nugget, messages=F)
    #ini.cov <- c(variog_fit$cov.pars[1] + variog_fit$nugget, variog_fit$cov.pars[2]) # c(var(data.x$data),10)
    ini.cov <- c(var(sim_data), 50)
    lik_fit <- likfit(data.geo, lik.method="ML", ini.cov.pars=ini.cov, cov.model=cov.model, fix.nugget=F, messages=F)
    #lik_fit <- likfit(as.geodata(data), lik.method="REML", ini.cov.pars=variog_fit, cov.model=cov.model, fix.nugget=F)
    # 3. Contruir matriz de covarianzas
    dist.x  <- as.matrix(dist(grilla))
    Sigma.x <- geoR::cov.spatial(dist.x, cov.model=cov.model, cov.pars=lik_fit$cov.pars)
    diag(Sigma.x) <- diag(Sigma.x) + lik_fit$nugget
    sg <- lik_fit$nugget + lik_fit$cov.pars[1] # varianza
    # 4. Invertir matriz de covarianzas
    if(cov.model=="gaussian"){
      # Eigen descomposition
      eigen_decom <- eigen(Sigma.x, symmetric=T)
      inv_diag <- 1/eigen_decom$values
      ISigma.x <- eigen_decom$vectors %*% diag(inv_diag) %*% t(eigen_decom$vectors)
    } else {
      ISigma.x <- solve(Sigma.x)
    }
    # Estimador GLS y su varianza
    mx <- as.numeric((colSums(ISigma.x) %*% sim_data) / sum(ISigma.x))
    vx <- 1/sum(ISigma.x)
    # Efective sample size ESS
    ess1 <- n*sum(diag(Sigma.x))/sum(Sigma.x) # Griffith 2005
    ess2 <- sum(sg*ISigma.x)                  # Vallejos 2014
    # Statistic
    stat1 <- (mx - 0) / sqrt(vx)                            # GLS estimator
    stat2 <- (mean(sim_data)- 0)/ sqrt(sum(Sigma.x)/n^2)    # Griffith 2005
    stat3 <- (mean(sim_data)- 0)/ sqrt(sum(Sigma.x)/ess1^2) # Griffith 2005 with n=n*
    stats <- c(stat1,stat2,stat3)
    # P-value
    for(j in 1:length(stats)){
      pv <- c(p.value1 = pnorm(stats[j], lower.tail=F),         # Distribucion normal
              p.value2 = pt(stats[j], df=ess1, lower.tail=F),   # Distribucion t con df = griffith
              p.value3 = pt(stats[j], df=ess1-1, lower.tail=F), # Distribucion t con df = griffith - 1
              p.value4 = pt(stats[j], df=ess2, lower.tail=F),   # Distribucion t con df = vallejos
              p.value5 = pt(stats[j], df=ess2-1, lower.tail=F)  # Distribucion t con df = vallejos - 1
      )
      PVALUE[[j]][i,] <- pv
    }
    # Almacenar resultados
    SIGMA[i]  <- sg
    PHI[i] <- lik_fit$cov.pars[2]
    STAT[i,] <- stats
    ESS[i,] <- c(ess1, ess2)
  }
  # save samples
  save_samples <- list("SIGMA"=SIGMA, "PHI"=PHI, "STAT"=STAT, "ESS"=ESS, "PVALUE"=PVALUE)
  saveRDS(save_samples, file=paste0("distri_ml_",label_model,"_n",n,"_phi",cov.pars[2],".rds")) # los guarda en la direccion donde este este script dentro del R project
}



# Analisis de muestras ----------------------------------------------------


# Load ata
DATA <- readRDS(file="Distribuciones empiricas/distri_ml_exp_n100_phi60.rds")
distri <- readRDS("sigma2_1_phi_30_OLS/distri_exp_n100.rds")
#list2env(distri, globalenv()) # Extrae los objectos y los pone en el eviroment
cov.model="exponential"
cov.pars=c(1,60)
nugget=0
fix.nugget=F

# Muestras atipicas
summary(distri$SIGMA)
quantile(distri$SIGMA, probs=c(seq(0.75,0.95,by=0.05), 0.97, 0.98,0.99))
which(distri$SIGMA > quantile(distri$SIGMA, probs=0.95))

# Datos
ind = 40 
dat <- DATA[[ind]]
View(distri$SIGMA[ind])
plot(as.geodata(dat))
cor(dat)
# Variograma empirico
variog_empirical <- variog(as.geodata(dat), estimator.type="classical")
plot(variog_empirical)
variog_empirical$n
variog_empirical$sd
eyefit(variog_empirical)
# Estimar variograma teorico (OLS)
vario_fit_equal <- variofit(variog_empirical, weights="equal", cov.model=cov.model, fix.nugget=fix.nugget, nugget=nugget)
vario_fit_npairs <- variofit(variog_empirical, weights="npairs", cov.model=cov.model, fix.nugget=fix.nugget, nugget=nugget)
vario_fit_cressie <- variofit(variog_empirical, weights="cressie", cov.model=cov.model, fix.nugget=fix.nugget, nugget=nugget)

vario_fit_equal$cov.pars
vario_fit_npairs$cov.pars
vario_fit_cressie$cov.pars
# Graficos
plot(variog_empirical)
lines(vario_fit_equal, col="red")   # equal
lines(vario_fit_npairs, col="blue") # npairs
lines(vario_fit_cressie, col="purple") # cressie
curve(1-cov.spatial(x,cov.model=cov.model, cov.pars=cov.pars), lwd=1.5, col=1, add=T)

# Estimar variograma teorico (ML y REML)
#ini.cov <- c(var(dat[,3]), 50)
ini.cov <- c(vario_fit_cressie$cov.pars[1] + vario_fit_cressie$nugget, vario_fit_cressie$cov.pars[2])
likfit_ml <- likfit(as.geodata(dat), lik.method="ML", ini.cov.pars=ini.cov, cov.model=cov.model, fix.nugget=F)
likfit_reml <- likfit(as.geodata(dat), lik.method="REML", ini.cov.pars=ini.cov, cov.model=cov.model, fix.nugget=F)
# fijar ml fijar nugget=0

print(paste("sigmaq:", round(likfit_ml$cov.pars[1],3),"tau:", round(likfit_ml$nugget,3),"phi:",round(likfit_ml$cov.pars[2],3)))
print(paste("sigmaq:", round(likfit_reml$cov.pars[1],3) ,"phi:",round(likfit_reml$cov.pars[2],3),"tau:", round(likfit_reml$nugget,3)))

plot(variog_empirical)
lines(likfit_ml, col="darkorange", lwd=2)
lines(likfit_reml, col="darkred", lwd=2)
curve(1-cov.spatial(x,cov.model=cov.model, cov.pars=cov.pars), lwd=1.5, col=1, add=T)
plot(variog_empirical, add=T, ylim=c(0,1))

?eyefit
?variofit
?variomodel

# muestra 15:
# Hay una alta correlacion entre el eje Y y los datos
# Estimar el variograma empirico con la funcion de perdida de Creassie mejora el ajuste respecto al verdadero covariograma



# Distribuciones empiricas ML ---------------------------------------------

distri <- readRDS("Distribuciones empiricas/distri_ml_exp_n100_phi60.rds")

# Resultados SIGMA
summary(distri$SIGMA)
quantile(distri$SIGMA, seq(0.75,0.95,by=0.05))
mean(distri$SIGMA)
sd(distri$SIGMA)
make_hist(distri$SIGMA)
hist(distri$SIGMA, freq=F)
#Resultados PHI
summary(distri$PHI)
quantile(distri$PHI, c(seq(0.75,0.95,by=0.05),0.97,0.98,0.99))
hist(distri$PHI, freq = F)
lines(density(distri$PHI))
mean(distri$PHI)
sd(distri$PHI)
make_hist(distri$PHI, symbol="phi")
# Promedio de las stats
colMeans(STAT)
apply(distri$STAT, 2, summary)

hist(STAT[,1]) 
hist(STAT[,2])
hist(STAT[,3], freq=F)
make_hist(distri$STAT[,1])
make_hist(distri$STAT[,2])
make_hist(distri$STAT[,3])

# # Encontrar los gl a partir de la distribucion empirica
# q95 <- quantile(STAT[,1], 0.95)
# dft <- function(x){return(abs(qt(0.95, df=x)-q95))}
# optim(par=colMeans(ESS)[1], fn=dft)
# qt(0.95, df=9)

# Promedio de los ess
colMeans(distri$ESS)
hist(distri$ESS[,1])
hist(distri$ESS[,2])
qt(0.95, df=colMeans(ESS)[1]-2)

# Promedio pvalores
apply(distri$PVALUE$stat_gls,2,function(x){mean(x<0.05)})
apply(distri$PVALUE$stat_griffith,2,function(x){mean(x<0.05)})
apply(distri$PVALUE$stat_griffith_ess,2,function(x){mean(x<0.05)}) # descartar


# Usara gl con ESS estimados
colMeans(distri$ESS)

hist(distri$STAT[,"stat_gls"], freq=F, breaks=30, ylim=c(0,0.4))
curve(dt(x, df=1), from=-5, to=5, add=T)

hist(distri$STAT[,"stat_griffith"], freq=F, breaks=30, ylim=c(0,0.4))
curve(dt(x, df=4), from=-5, to=5, add=T)


