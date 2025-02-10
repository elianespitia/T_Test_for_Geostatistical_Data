
# Librerias
library(geoR)
library(ggplot2)
library(parallel)
library(pbapply)
library(openxlsx)

# Directorios
wd  = "C:/Users/driao/Documents/ARTICULOS/One and two sample test for geostatistical data"
wdf = paste0(wd, "/Funciones")
wdr = paste0(wd, "/Resultados")

# Crear cluster
cl <- makeCluster(getOption("cl.cores", 8))
clusterEvalQ(cl, {
  
  # Librerias
  library(geoR) 
  library(ggplot2)
  library(parallel)
  library(pbapply)
  
  # Funciones
  source(file=paste0(wdf,"/random_coords.R"))
  source(file=paste0(wdf,"/power.t.test.one.sample.R"))
  source(file=paste0(wdf,"/power.t.test.two.sample.R"))
  source(file=paste0(wdf,"/make_plot.R"))
  # list.files(path=wdf, full.names=T)
  
  })

stopCluster(cl)
rm(cl)

mu_range=seq(from=-2, to=2, by=0.05)
type = "one.sample"
phi = c(10,30,60,90)
size = c(20,50,100,200)
export_path = wdr

# Funcion
power.curve.t.test <- function(mu_range=seq(from=-2, to=2, by=0.05), type=NULL, phi=NULL, size=NULL, export_path){
  
  # Verificar las entradas de los parametros
  if(type ==)
  type <- match.arg(tipo, choices = c("one.sample", "two.sample"))
  if(is.null(phi) & is.null(size)){stop("Los parametros phi y size no pueden ser ambos nulos")}
  
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # phi y size no nulos --------------------------------------------------------
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  OUT <- list() # Create a output list
  
  if( !is.null(phi) && !is.null(size) && length(phi) > 1 ){
    
    if(type == "one.sample"){
      
      for(i in size){
        n <- size[i] # Tamño de muestra
        out <- data.frame(x=mu_range) # Create a dataframe
        
        for(j in phi){
          out[paste0("phi",phi[j])] <- pbsapply(cl=cl, mu_range, power.t.test.one.sample,
                                                mu0=mu0, n, model, sigma, phi=phi[i], pepita=0, k=0.5,
                                                sims=1000, alpha=0.05, ha="two.sided")
          }
        
        # Save and export data
        OUT[[paste0("n", n)]] <- out
        write.xlsx(out, file=paste0(export_path, "/", paste("os", model, paste0("n",n), ha, sep="_"), ".xlsx"))
      } # cierra for size
    
    } else {
  
      for(i in size){
        n <- size[i] # Tamño de muestra
        out <- data.frame(x=mu_range) # Create a dataframe
        
        for(j in phi){
          out[paste0("phi",phi[j])] <- pbsapply(cl=cl, mu_range, power.t.test.two.sample,
                                                mu1, mu2, n1, n2, model, sigma, phi, pepita=0, k=0.5,
                                                sims=1000, alpha=0.05, ha="two.sided", varequal=T)
        }
        
        # Save and export data
        OUT[[paste0("n", n)]] <- out
        write.xlsx(out, file=paste0(export_path, "/", paste("ts", model, paste0("n",n), ha, sep="_"), ".xlsx"))
      } # cierra for size
        
    } # close else type
  
  }
  
  if( !is.null(phi) && !is.null(size) && length(phi) == 1 ){

    if(type == "one.sample"){
      
      for(i in size){
        out <- data.frame(x=mu_range) # Create a dataframe
        
        out[paste0("n", size[i])] <- pbsapply(cl=cl, mu_range, power.t.test.one.sample,
                                        mu0=mu0, n=size[i], model, sigma, phi=phi, pepita=0, k=0.5,
                                        sims=1000, alpha=0.05, ha="two.sided")
      }
      
      # Save and export data
      OUT[[paste0("phi", phi)]] <- out
      write.xlsx(out, file=paste0(export_path, "/", paste("os", model, paste0("phi",phi), ha, sep="_"), ".xlsx"))
    
    } else {
      
      out <- data.frame(x=mu_range) # Create a dataframe
        
        for(i in size){
          out[paste0("n", n)] <- pbsapply(cl=cl, mu_range, power.t.test.two.sample,
                                                mu1, mu2, n1, n2, model, sigma, phi, pepita=0, k=0.5,
                                                sims=1000, alpha=0.05, ha="two.sided", varequal=T)
        }
        
        # Save and export data
        OUT[[paste0("n", n)]] <- out
        write.xlsx(out, file=paste0(export_path, "/", paste("ts", model, paste0("n",n), ha, sep="_"), ".xlsx"))
      } # cierra for size
      
    } # close else type
  
  
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # phi nulo -------------------------------------------------------------------
  #:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  if(is.null(phi)){
    phi = 10
    names_out = paste0("n", size)
    message("Parametro phi no ingresado. Por defecto se toma phi=10")

    

  }

  # OUTPUT
  return(OUT)

}



# Ejemplo
power.curve.t.test(mu_range=seq(from=-2, to=2, by=0.05), type, phi=NULL, size=c(20,50,100,200))





