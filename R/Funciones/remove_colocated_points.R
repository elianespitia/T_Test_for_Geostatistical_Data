
# Function to remove co-localed points

remove_colocated_points <- function(coords, d) {
  
  # Ensure input is a matrix with two columns
  if (!is.matrix(coords) || ncol(coords) != 2) {
    stop("Input must be a matrix with two columns representing coordinates.")
  }
  
  # Compute pairwise distances
  dist_matrix <- as.matrix(dist(coords))
  
  # Identify points that are within distance d
  to_remove <- rep(FALSE, nrow(coords))
  
  for (i in 1:(nrow(coords) - 1)) {
    if (!to_remove[i]) {  # If the point has not already been marked for removal
      close_points <- which(dist_matrix[i, ] < d & !to_remove)
      close_points <- close_points[close_points != i] # Exclude itself
      
      if (length(close_points) > 0) {
        to_remove[close_points] <- TRUE # Mark duplicates for removal
      }
    }
  }
  
  # Return filtered coordinates
  return(coords[!to_remove, , drop = FALSE])
}

#:::::::::
# Example
#:::::::::

# source("random_coords.R")
# grilla <- random_coords(n=20)
# plot(grilla)
# dist.x <- as.matrix(dist(grilla))
# 
# new_grilla <- remove_colocated_points(coords=grilla, d=5)
# new_grilla
# 
# par(mfrow=c(1,2))
# plot(grilla); plot(new_grilla)
