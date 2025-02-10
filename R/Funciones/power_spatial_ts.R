
#:::::::::::::::::::::::::::::::
# Power for two-sample test 
#:::::::::::::::::::::::::::::::

power_spatial_ts <- function(mux, muy, nx, ny, model.x, model.y, cov.pars.x, cov.pars.y,
                             nugget.x, nugget.y, kx=0.5, ky=0.5,
                             sims=1000, alpha=0.05, ha="greater", test){
  
  # Select the test function
  test_function <- switch(test,
                          "spatial.z.test" = spatial.z.test.ts.fast,
                          "spatial.t.test" = spatial.t.test.ts.fast,
                          stop("Invalid test selected")
  )
  
  # Simulations
  set.seed(1) # con este semilla sal NaN en la primera simulacion. dejar o no??
  pvalues <- replicate(n=sims, {
    
    # Generate random coordinates
    grilla.x <- random_coords(nx)
    grilla.y <- random_coords(ny)
    
    # Simulate data
    sim_data_x <- geoR::grf(n = nx, grid = grilla.x, nsim = 1, cov.model = model.x,
                            cov.pars = cov.pars.x, nugget = nugget.x, kappa = kx, 
                            mean = mux, messages=F, method="svd")$data
    sim_data_y <- geoR::grf(n = ny, grid = grilla.y, nsim = 1, cov.model = model.y,
                            cov.pars = cov.pars.y, nugget = nugget.y, kappa = ky,
                            mean = muy, messages=F, method="svd")$data
    data.x <- cbind(grilla.x, sim_data_x)
    data.y <- cbind(grilla.y, sim_data_y)
    
    # Apply the selected test function
    test_function(data.x=data.x, data.y=data.y, mu0=mu0, model.x=model.x, model.y=model.y,
                  cov.pars.x=cov.pars.x, cov.pars.x=cov.pars.y, nugget.x=nugget.x, nugget.y=nugget.y, alternative=ha)
  })
  
  # OUTPUT
  power <- mean(pvalues < alpha)
  return(power)
  
} # Close function