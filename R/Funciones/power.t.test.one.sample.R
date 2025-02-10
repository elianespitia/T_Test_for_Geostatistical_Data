


#' Potencia de la pruba-t para una muestra de datos geoespaciales
#'
#' @param mu Valor para el cual calcular la funcion de potencia.
#' @param mu0 Valor de contraste de la prueba
#' @param n Tamaño de muestra.
#' @param model Modelo covarianza.
#' @param sigma Varianza
#' @param phi Rango
#' @param pepita Efecto pepita.
#' @param k Parametro Kappa para los modelos
#' @param sims Numero de simulaciones
#' @param alpha Nivel de significancia de la prueba t
#' @param ha a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less".
#'
#' @return Valor de la potencia para mu
#' @export
#'
#' @examples

power.t.test.one.sample <- function(mu, mu0, n, model, sigma, phi, pepita=0, k=0.5,
                                    sims=1000, alpha=0.05, ha="two.sided"){
  
  # Mensaje para que el usuario sepa que para que modelos el parametro kappa es necesario.
  if(model %in% c("matern","powered.exponential","cauchy","gencauchy","gneiting.matern")){
    message(paste0('El modelo de covarianza que suministro (', model, ') necesita del parametro kappa (k).\nPor defecto k=0.5'))
  }
  
  # Vetor para almacenar los pvalues
  pvalues <- vector(mode="numeric", length=sims)
  
  # Simulaciones
  for(i in 1:sims){
    
    # Obtener coordenadas aleatorias
    grilla <- random_coords(n)
    
    # Simular datos geoestadisticos
    sim_data <- grf(n=n, grid=grilla, nsim=1,
                    cov.model=model, cov.pars=c(sigma^2, phi), nugget=pepita, kappa=k,
                    mean=mu, messages=F)$data
    
    # Prueba-t usual
    pvalues[i] <- t.test(x=sim_data, mu0=mu0, alternative=ha)$p.value
    
  } # cierra for
  
  # Calcular potencia
  power <- mean(pvalues < alpha)
  
  # Output
  return(power)
  
} # Fin de la funcion!!

# Librerias necesarias ---------------------------------------------------------------

#library(geoR)


# Ejemplo ------------------------------------------------------------------

# set.seed(123)
# power.t.test.one.sample(mu=0, mu0=0, n=20, model="pure.nugget", sigma=1, phi=10, sims=1000)

# Parametros --------------------------------------------------------------

# mu=0
# mu0=0
# n=20
# model="exponential"
# sigma=1
# phi=10
# pepita=0
# k=0.5
# sims=1000
# alpha=0.05
# ha="two.sided"

