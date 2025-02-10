
# Clean enviroment
rm(list=ls())

# Current directory
getwd()

# 0. Load Functions ---------------------------------------------------------------

# Libraries
library(geoR)       
library(parallel)
library(readxl)
library(openxlsx)

# Funtions
source("Funciones/random_coords.R")
source("Funciones/power_spatial_os.R")
source("Funciones/power_spatial_ts.R")
source("Funciones/spatial.z.test.fast.R")
source("Funciones/spatial.t.test.fast.R")



# 1. Set up Cluster ----------------------------------------------------------

# In this section we set up a cluster for running in parallel

detectCores()
num_cl = 4
cl = makeCluster(getOption("cl.cores", num_cl))
clusterEvalQ(cl, {
  
  # Libraries
  library(geoR)
  
  # Funtions
  source("Funciones/random_coords.R")
  source("Funciones/power_spatial_os.R")
  #source("Funciones/power_spatial_ts")
  source("Funciones/spatial.z.test.fast.R")
  #source("Funciones/spatial.t.test.fast.R")
})


# 2. Simulation Parameters ---------------------------------------------------

sims = 1000
n_seq   = c(20, 50, 100)
phi_seq = c(10, 50, 80, 100)
tau_sqe = c(0, 0.1, 1, 10)
models  = c("exponential", "gaussian", "spherical")

ha = "greater"
mu_seq <- if (ha == "greater"){
  seq(0, 2, by = 0.1)
  } else if (ha == "two.sided"){
    seq(-2, 2, by = 0.1)
    } else {
      stop("Invalid ha selected")
    }


# 3. T test ---------------------------------------------------------------



# 4. Spatial Z test -------------------------------------------------------

## 4.1 One sample ----------------------------------------------------------

test = "spatial.z.test" # type of test
output <- list() # List to save results

time.start <- Sys.time()
for(model in models){
  for(n in n_seq){
    for(tau in tau_sqe){
      
      # Data output
      out <- data.frame(x=mu_seq)
      for(phi in phi_seq){
        out[paste0("phi", phi)] <- parSapply(cl=cl, mu_seq, power_spatial_os, mu0=0, n=n,
                                             model=model, cov.pars=c(1,phi), nugget=tau, 
                                             sims=sims, ha=ha, test=test)
      }
      
      # Save data
      lab <- paste0(substr(model, start=1, stop=3),"_n",n,"_tau",tau)
      output[[lab]] <- out
      print(paste("Finish", lab))
    }
  }
}
time.final <- Sys.time() - time.start
beepr::beep()

# sims=1000  y ha="greater" => 1.70 hours
# sims=5000  y ha="greater" =>
# sims=10000 y ha="greater" =>

# Save results
saveRDS(output, file=paste0("Simulaciones/ztest_os_",ha,".rds"))


## 4.2 Two samples ----------------------------------------------------------

test = "spatial.z.test" # type of test
output <- list() # List to save results



# 5. Spatial T test -------------------------------------------------------

## 5.1 One sample ----------------------------------------------------------

test = "spatial.t.test" # type of test
output <- list() # List to save results

time.start <- Sys.time()
for(model in models){
  for(n in n_seq){
    for(tau in tau_sqe){
      
      # Data output
      out <- data.frame(x=mu_seq)
      for(phi in phi_seq){
        out[paste0("phi", phi)] <- parSapply(cl=cl, mu_seq, power_spatial_os, mu0=0, n=n,
                                             model=model, cov.pars=c(1,phi), nugget=tau, 
                                             sims=sims, ha=ha, test=test)
      }
      
      # Save data
      lab <- paste0(substr(model, start=1, stop=3),"_n",n,"_tau",tau)
      output[[lab]] <- out
      print(paste("Finish", lab))
    }
  }
}
time.final <- Sys.time() - time.start
beepr::beep()

# sims=1000  y ha="greater" => 

# Save results
saveRDS(output, file=paste0("Simulaciones/ttest_os_",ha,".rds"))


## 5.2 Two samples ----------------------------------------------------------

test = "spatial.t.test" # type of test
output <- list() # List to save results



