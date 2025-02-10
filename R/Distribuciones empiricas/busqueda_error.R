


# Libraries
library(geoR)
library(scatterplot3d)
source("../Funciones/random_coords.R") # ../ sube un nivel
source("../Funciones/make_hist.R")

# Parametros
n=20
#nseq=c(10,30,50,100,200,350)
cov.model="exponential"
# short_name_model = c("exp","gauss","sphe"); names(short_name_model) <- c("exponential","gaussian","spherical")
# label_model = short_name_model[cov.model]
cov.pars=c(1,30)
nugget=0
fix.nugget=F
sims=1000

# Almacenamiento
SIGMA  <- vector(mode="numeric", length=sims)
PHI    <- vector(mode="numeric", length=sims)
DATA   <- vector(mode="list", length=sims)
# STAT   <- matrix(NA, nrow=sims, ncol=3); colnames(STAT) <- paste0("stat_",c("gls","griffith","griffith_ess"))
# ESS    <- matrix(NA, nrow=sims, ncol=2); colnames(ESS) <- c("griffith","vallejos")
# PVALUE <- list(matrix(NA, nrow=sims, ncol=5, dimnames=list(c(1:sims), c("normal","griffith","griffith_1","vallejos","vallejos_1"))),
#                matrix(NA, nrow=sims, ncol=5, dimnames=list(c(1:sims), c("normal","griffith","griffith_1","vallejos","vallejos_1"))),
#                matrix(NA, nrow=sims, ncol=5, dimnames=list(c(1:sims), c("normal","griffith","griffith_1","vallejos","vallejos_1")))
# )
# names(PVALUE) <- colnames(STAT)


# Simulations
set.seed(123)
for(i in 1:sims){
  #::::::::::::::::::::::::::::::::::::::::
  # Simular campo aleatorio
  #::::::::::::::::::::::::::::::::::::::::
  # Get random coords
  grilla <- random_coords(n)
  # Simulate geostatistical data
  sim_data <- geoR::grf(grid=grilla, nsim=1,cov.model=cov.model,
                        cov.pars=cov.pars, nugget=nugget,
                        mean=0, messages=F, method="eigen")$data
  data.x <- cbind(grilla, sim_data)
  DATA[[i]] <- data.x
  data.x <- geoR::as.geodata(data.x)
  #::::::::::::::::::::::::::::::::::::::::
  # Estimar Matriz de Covarianzas
  #::::::::::::::::::::::::::::::::::::::::
  # 1. Estimar el variograma empirico
  variog_empirical <- geoR::variog(data.x, messages=F) # posiblemente modificar
  
  # 2. Ajustar un modelo teorico
  variog_fit <- geoR::variofit(variog_empirical, cov.model=cov.model, fix.nugget=fix.nugget, nugget=nugget, messages=F)
  
  # 3. Contruir matriz de covarianzas
  # dist.x  <- as.matrix(dist(data.x$coords))
  # Sigma.x <- geoR::cov.spatial(dist.x, cov.model=cov.model, cov.pars=variog_fit$cov.pars)
  # diag(Sigma.x) <- diag(Sigma.x) + variog_fit$nugget
  # sg <- variog_fit$nugget + variog_fit$cov.pars[1] # varianza
  # 4. Invertir matriz de covarianzas
  # if(cov.model=="gaussian"){
  #   # Eigen descomposition
  #   eigen_decom <- eigen(Sigma.x, symmetric=T)
  #   inv_diag <- 1/eigen_decom$values
  #   ISigma.x <- eigen_decom$vectors %*% diag(inv_diag) %*% t(eigen_decom$vectors)
  # } else {
  #   ISigma.x <- solve(Sigma.x)
  # }
  # # Estimador GLS y su varianza
  # mx <- as.numeric((colSums(ISigma.x) %*% sim_data) / sum(ISigma.x))
  # vx <- 1/sum(ISigma.x)
  # # Efective sample size ESS
  # ess1 <- n*sum(diag(Sigma.x))/sum(Sigma.x) # Griffith 2005
  # ess2 <- sum(sg*ISigma.x)                  # Vallejos 2014
  # # Statistic
  # stat1 <- (mx - 0) / sqrt(vx)                            # GLS estimator
  # stat2 <- (mean(sim_data)- 0)/ sqrt(sum(Sigma.x)/n^2)    # Griffith 2005
  # stat3 <- (mean(sim_data)- 0)/ sqrt(sum(Sigma.x)/ess1^2) # Griffith 2005 with n=n*
  # stats <- c(stat1,stat2,stat3)
  # # P-value
  # for(j in 1:length(stats)){
  #   pv <- c(p.value1 = pnorm(stats[j], lower.tail=F),         # Distribucion normal
  #           p.value2 = pt(stats[j], df=ess1, lower.tail=F),   # Distribucion t con df = griffith
  #           p.value3 = pt(stats[j], df=ess1-1, lower.tail=F), # Distribucion t con df = griffith - 1
  #           p.value4 = pt(stats[j], df=ess2, lower.tail=F),   # Distribucion t con df = vallejos
  #           p.value5 = pt(stats[j], df=ess2-1, lower.tail=F)  # Distribucion t con df = vallejos - 1
  #   )
  #   PVALUE[[j]][i,] <- pv
  # }
  
  # Almacenar resultados
  SIGMA[i]  <- variog_fit$cov.pars[1]
  PHI[i] <- variog_fit$cov.pars[2]
  # STAT[i,] <- stats
  # ESS[i,] <- c(ess1, ess2)
}

# save rds
saveRDS(DATA, file="grilla_exp_n20_phi30.rds") # los guarda en la direccion donde este este script dentro del R project
# load rds
DATA <- readRDS(file="grilla_exp_n20.rds")

# Identificar estimaciones extrañas
summary(SIGMA)
quantile(SIGMA, probs=c(seq(0.75,0.95,by=0.05), 0.97, 0.98,0.99))
sum(SIGMA>quantile(SIGMA,0.97))
hist(SIGMA)
hist(SIGMA[SIGMA < quantile(SIGMA, probs=0.95)])

which(SIGMA > quantile(SIGMA, probs=0.95))


# Datos
data <- DATA[[49]]
plot(as.geodata(data))
# Regresion con las covariables
my.lm <- lm(data[,3] ~ data[,1] + data[,2]); summary(my.lm)
s3d$plane3d(my.lm)
s3d <- scatterplot3d(data, type = "h")
s3d$plane3d(my.lm)
# Correlacion
cor(data)
max(dist(data[,c(1,2)]))
# Variograma empirico
variog_empirical <- variog(as.geodata(data), estimator.type="classical")
plot(variog_empirical)
variog_empirical$n
variog_empirical$sd
eyefit(variog_empirical)
# Estimar variograma teorico (OLS)
vario_fit_equal <- variofit(variog_empirical, weights="equal", cov.model=cov.model, ini.cov.pars = c(1,30) , fix.nugget=fix.nugget, nugget=nugget)
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
likfit_ml <- likfit(as.geodata(data), lik.method="ML", ini.cov.pars=c(var(data[,3]), 70) , cov.model=cov.model, fix.nugget=F)
likfit_ml2 <- likfit(as.geodata(data), lik.method="ML", ini.cov.pars=c(var(data[,3]), 70), cov.model=cov.model, fix.nugget=F)
likfit_reml <- likfit(as.geodata(data), lik.method="REML", ini.cov.pars=c(50, 100), cov.model=cov.model, fix.nugget=F)

likfit_ml$cov.pars
likfit_ml2$cov.pars
likfit_reml$cov.pars

#plot(variog_empirical)
lines(likfit_ml, col="darkorange", lwd=2)
lines(likfit_ml2, col="darkgreen", lwd=2 )
lines(likfit_reml, col="darkred", lwd=2)
curve(1-cov.spatial(x,cov.model=cov.model, cov.pars=cov.pars), lwd=1.5, col=1, add=T)


?eyefit
?variofit
?variomodel

# muestra 15:
# Hay una alta correlacion entre el eje Y y los datos
# Estimar el variograma empirico con la funcion de perdida de Creassie mejora el ajuste respecto al verdadero covariograma




