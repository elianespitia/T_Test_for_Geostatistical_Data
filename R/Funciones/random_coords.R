
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Funcion para coordenadas aleatorias -------------------------------------
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Title Simulate random coords on a rectangular area.
#'
#' @param n Number of simulated points.
#' @param xrange Range of x-axis. 
#' @param yrange Range of y-axis.
#'
#' @return nx2 matriz with simulated coordinates
#' @export
#'
#' @examples random_coords(n=10)

random_coords <- function(n, xrange=c(0,100), yrange=c(0,100)){
  
  # Muetreas coordenadas de una uniforme sobre el area
  x <- runif(n, min=xrange[1], max=xrange[2])
  y <- runif(n, min=yrange[1], max=yrange[2])
  
  # Output
  return(cbind(x,y))
}
