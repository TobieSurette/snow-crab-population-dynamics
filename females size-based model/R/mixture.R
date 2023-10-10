mixture <- function(x, f){
   # Parse data observations:
   if (missing(f)){
      f <- NULL
      if (!is.null(names(x))){
         if (all(gsub("[0-9.]", "", names(x)) == "") & !any(names(x) != "")){
            f <- as.numeric(x)
            x <- as.numeric(names(x))
         }
      }
      if (is.null(f)){
         x <- table(x)
         f <- as.numeric(x)
         x <- as.numeric(names(x))
      }
   }else{
      if (length(x) != length(f)) stop("'x' and 'f' must be the same size.") 
      r <- aggregate(f, by = list(x), sum)
      f <- r[,2]
      x <- r[,1]
   } 
  
   v <- list()
   v$data$x <- x
   v$data$f <- f
   
   class(v) <- c("mixture", "list")
   
   return(v)
}

parameters.mixture <- function(x, ...){
   # Instar proportions:
   p <- theta[grep("logit", names(theta))]
   p <- exp(p) / (1 + sum(exp(p)))  
   k <- length(p) + 1  # Number of instars.
   p[k] <- 1 - sum(p)
  
   # Instar mean sizes:
   mu <- theta["mu0"]          
   for (i in 2:k) mu[i] <- theta["intercept"] + theta["slope"] * mu[i-1]
  
   # Instar standard errors:
   sigma <- exp(theta["log.sigma"]) * mu
   
   v <- data.frame(p = p, mu = mu, sigma = sigma)
   
   return(v)
}

density.mixture <- function(x, x0, ...){
   
   # Calculate mixture density:
   ll <- rep(0, length(x))
   for (i in 1:k){
     ll <- ll + p[i] * dnorm(x, mu[i], sigma[i])
   }
}

# Mixture stats:
summary.mixture <- function(x, ...){
  
  
}
  
# Fit Hiatt mixture density to immatures:
loglike.mixture <- function(theta, x, fixed){
  if (!missing(fixed)) theta <- c(theta, fixed)
  
  # Instar proportions:
  p <- theta[grep("logit", names(theta))]
  p <- exp(p) / (1 + sum(exp(p)))  
  k <- length(p) + 1  # Number of instars.
  p[k] <- 1 - sum(p)
  
  # Instar mean sizes:
  mu <- theta["mu0"]          
  for (i in 2:k) mu[i] <- theta["intercept"] + theta["slope"] * mu[i-1]
  
  # Instar standard errors:
  sigma <- exp(theta["log.sigma"]) * mu
  
  # Prepare data:
  f <- as.numeric(x)
  x <- as.numeric(names(x))
  
  # Calculate mixture density:
  ll <- rep(0, length(x))
  for (i in 1:k){
    ll <- ll + p[i] * dnorm(x, mu[i], sigma[i])
  }
  
  return(sum(- f * log(ll)))
}


# Instar proportions:
p <- theta[grep("logit", names(theta))]
p <- exp(p) / (1 + sum(exp(p)))  
k <- length(p) + 1  # Number of instars.
p[k] <- 1 - sum(p)

# Instar mean sizes:
mu <- theta["mu0"]          
for (i in 2:k) mu[i] <- theta["intercept"] + theta["slope"] * mu[i-1]
names(mu) <- as.roman(4:9)

# Instar standard errors:
sigma <- exp(theta["log.sigma"]) * mu