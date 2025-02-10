
n=20
sims=50
DATAX <- vector(mode="list", l=sims)
PVALUE <- vector(mode="numeric", l=sims)
set.seed(82) # con este semilla sal NaN en la primera simulacion
for(i in 1:sims){
  
  # Select the test function
  test_function <- switch(test,
                          "spatial.z.test" = spatial.z.test.os.fast,
                          "spatial.t.test" = spatial.t.test.os.fast,
                          stop("Invalid test selected")
  )
  
  # Generate random coordinates
  grilla <- random_coords(n)
  
  # Simulate data
  sim_data <- geoR::grf(n = n, grid = grilla, nsim = 1, cov.model = model,
                        cov.pars = cov.pars, nugget = nugget, kappa = k,
                        mean = mu, messages=F, method="svd")$data
  data.x <- cbind(grilla, sim_data)
  DATAX[[i]] <- data.x 
  
  # Apply the selected test function
  PVALUE[i] <- test_function(data.x=data.x, mu0=mu0, model.x=model, cov.pars.x=cov.pars, alternative=ha)
}

# Ver que sucede con la muestra 37
data.x <- DATAX[[2]]
plot(as.geodata(data.x))
dist.x <- dist(data.x[,c(1,2)]); View(as.matrix(dist.x))
min(as.vector(dist.x))

dist.x[lower.tri(x, diag = FALSE)]

upper.tri(x, diag = FALSE)
