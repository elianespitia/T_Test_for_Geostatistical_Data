

#' Potencia de la prueba t en presencia de correlacion espacial
#'
#' @param mu1 
#' @param mu2 
#' @param n1 
#' @param n2 
#' @param model 
#' @param sigma 
#' @param phi 
#' @param pepita 
#' @param k 
#' @param sims 
#' @param alpha 
#' @param ha a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#' @param varequal a logical variable indicating whether to treat the two variances as being equal. If TRUE then the pooled variance is used to estimate the variance otherwise the Welch (or Satterthwaite) approximation to the degrees of freedom is used.
#'
#' @return
#' @export
#'
#' @examples

power.t.test.two.sample <- function(mu1, mu2, n1, n2, model, sigma, phi, pepita=0, k=0.5,
                                    sims=1000, alpha=0.05, ha="two.sided", varequal=T){
  
  # Mensaje para que el usuario sepa que para que modelos el parametro kappa es necesario.
  if(model %in% c("matern","powered.exponential","cauchy","gencauchy","gneiting.matern")){
    message(paste0('El modelo de covarianza que suministro (', model, ') necesita del parametro kappa (k).\nPor defecto kappa = 0.5'))
  }
  
  # Change the interpretation of altenative in t.test function
  change_ha = list(greater="less", less="greater")
  ha = change_ha[[ha]]
  # Now "greater" is that y has a larger mean than x
  # Now "less"    is that x has a larger mean than y
  
  # Vetor para almacenar los pvalues
  pvalues <- vector(mode="numeric", length=sims)
  
  # Simulacion
  for(i in 1:sims){
    
    # Obtener coordenadas aleatorias
    grilla1 <- random_coords(n1)
    grilla2 <- random_coords(n2)
    
    # Simular muestras geoestadisticas de 2 poblaciones
    sim_data1 <- grf(n=n1, grid=grilla1, nsim=1,
                     cov.model=model, cov.pars=c(sigma^2, phi), nugget=pepita, kappa=k,
                     mean=mu1, messages=F)$data
    
    sim_data2 <- grf(n=n2, grid=grilla2, nsim=1,
                     cov.model=model, cov.pars=c(sigma^2, phi), nugget=pepita, kappa=k,
                     mean=mu2, messages=F)$data
    
    # Prueba-t para 2 muestras
    pvalues[i] <- t.test(x=sim_data1, y=sim_data2, var.equal=varequal, conf.level=1-alpha, alternative=ha)$p.value
  }
  
  # Output: Potencia
  return(mean(pvalues < alpha))
  
} # Fin  de la funcion!!



# Prueba ------------------------------------------------------------------

# set.seed(123)
# power.t.test.two.sample(mu1=0, mu2=0.5, n1=20, n2=20, model="pure.nugget", sigma=1, phi=10, sims=1000)


# Parametros --------------------------------------------------------------

# mu1=0
# mu2=0
# n1=20
# n2 = 20
# model="exponential"
# sigma=1
# phi=10
# pepita=0
# k=0.5
# sims=1000
# alpha=0.05
# ha="greater"
# varequal=T


