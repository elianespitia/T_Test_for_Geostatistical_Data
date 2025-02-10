
#:::::::::::::::::::::::::::::::
# Power for one-sample test 
#:::::::::::::::::::::::::::::::

power_spatial_os <- function(mu, mu0, n, model, cov.pars, nugget=0, k=0.5,
                             sims=1000, alpha=0.05, ha="greater", test){
  
  # Select the test function
  test_function <- switch(test,
                          "spatial.z.test" = spatial.z.test.os.fast,
                          "spatial.t.test" = spatial.t.test.os.fast,
                          stop("Invalid test selected")
  )
  
  # Simulations
  set.seed(1) # con este semilla sal NaN en la primera simulacion, dejar o no??
  pvalues <- replicate(n=sims, {
    
    # Generate random coordinates
    grilla <- random_coords(n)
    
    # Simulate data
    sim_data <- geoR::grf(n = n, grid = grilla, nsim = 1, cov.model = model,
                          cov.pars = cov.pars, nugget = nugget, kappa = k,
                          mean = mu, messages=F, method="svd")$data
    data.x <- cbind(grilla, sim_data)
    
    # Apply the selected test function
    test_function(data.x=data.x, mu0=mu0, model.x=model, cov.pars.x=cov.pars, nugget.x=nugget, alternative=ha)
  })
  
  # OUTPUT
  power <- mean(pvalues < alpha)
  return(power)
  
} # Close function

#::::::::::
# Examples
#::::::::::

# source("random_coords.R")
# source("spatial.z.test.fast.R")
# source("spatial.t.test.fast.R")
# mu=0
# mu0=0
# n=100
# model="gaussian"
# cov.pars = c(1, 50)
# nugget = 0
# k=0.5
# sims=1000
# alpha=0.05
# ha="greater"
# test="spatial.z.test"

# time.taken <- Sys.time()
# power_spatial_os(mu=0, mu0=0, n=20, model="gaussian", cov.pars=c(1, 50),
#                  nugget=0, sims=10000, test="spatial.z.test")
# time.taken <- Sys.time() - time.taken

#:::::::::::::::::::::
# Resumen de corridas
#:::::::::::::::::::::

# sims=10_000 y test="spatial.z.test"
# for model="gaussian", phi=50,  n=20  => 12.31 sec
# for model="gaussian", phi=50,  n=50  => 32.86 sec
# for model="gaussian", phi=50,  n=100 => 2.26  min


