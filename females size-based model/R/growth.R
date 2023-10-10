# Growth matrix for a single linear growth increment:
growth.matrix <- function(x, theta, ymax, ...){
  ux <- sort(unique(x))
  xmax <- max(ux)
  dx <- min(diff(ux))
  x0 <- seq(min(ux), max(ux), by = dx)
  
  # Define growth means and standard errors:
  m <- theta["intercept"] + theta["slope"] * x0 
  s <- exp(theta["log.sigma"]) * m
  
  
  # Calculate corresponding gamma parameters:
  phi <- s^2 / m # Scale parameter.
  k <- m^2 / s^2 # Shape parameter.
  
  # Define growth output vector:
  if (missing(ymax)) ymax <- xmax + max(m + 3 * s) 
  
  # Map growth increments onto growth matrix:
  y0 <- seq(min(x0), ymax, by = dx)
  ymax <- y0[length(y0)]
  G <- matrix(0, nrow = length(x0), ncol = length(y0))
  dimnames(G) <- list(x = x0, y = y0)
  for (i in 1:length(x0)){
    z <- seq(x0[i], as.numeric(y0[length(y0)-1]), by = dx)
    G[i,as.character(z)] <- pgamma(z-x0[i]+dx/2, k[i], 1/phi[i]) - pgamma(z - x0[i] - dx/2, k[i], 1/phi[i]) 
    G[i,as.character(y0[length(y0)])] <- 1 - pgamma(y0[length(y0)] - x0[i] - dx/2, k[i], 1/phi[i]) 
  }
  
  return(G)
}


