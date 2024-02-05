# Tools for fitting a growth-structured mixture model to crab size data.

# Find initial parameter values:
init.smix <- function(x, k = 6, ...){
   theta <- c(slope = 1.25, 
              intercept = 2.7, 
              mu0 = 9, 
              log.sigma0 = -2, 
              log.sigma = -0.5, 
              logit.p = rep(1, k))
 
   return(theta)
}

# Parse structured mixture parameter vector:
parse.theta <- function(theta){
   # Instar proportions:
   p <- theta[grep("logit", names(theta))]
   p <- exp(p) / (1 + sum(exp(p)))  
   k <- length(p) + 1  # Number of instars.
   p[k] <- 1 - sum(p)
   
   # Instar mean sizes:
   mu <- theta["mu0"]          
   for (i in 2:k) mu[i] <- theta["intercept"] + theta["slope"] * mu[i-1]
   
   # Instar standard errors:
   #sigma <- exp(theta["log.sigma"]) * mu
   sigma <- exp(theta["log.sigma0"])
   for (i in 2:k) sigma[i] <- sqrt((theta["slope"] * sigma[i-1])^2 + exp(theta["log.sigma"])^2)
   
   # Create list of mixture proportions and moments:
   v <- data.frame(p = p, mu = mu, sigma = sigma)
   rownames(v) <- NULL
   
   # Label rows with instar numbers:
   rownames(v) <- as.roman(4:(3+nrow(v)))
   
   return(v) 
}

# Structured mixture model object:
smix <- function(x, f, theta = init.smix()){
   if ("smix" %in% class(x)){
      if (!missing(theta)){
         if (is.null(theta)) x$theta <- theta else x$theta[names(theta)] <- theta
      }
      return(x)
   }
   
   # Parse data observations:
   if (missing(f)){
      f <- NULL
      if (!is.null(names(x))){
         if (all(gsub("[0-9.]", "", names(x)) == "") & !any(names(x) == "")){
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
   v$theta <- theta
   
   class(v) <- c("smix", "list")
   
   return(v)
}

# Structured mixture density function:
dsmix <- function(x, x0, collapse = TRUE, ...){
   if (missing(x0)) x0 <- x$data$x 
  
   # Parse parameter values:
   v <- parse.theta(x$theta)

   # Calculate mixture density:
   ll <- matrix(0, nrow = length(x0), ncol = nrow(v))
   
   for (j in 1:nrow(v)) ll[,j] <- v$p[j] * dnorm(x0, v$mu[j], v$sigma[j])
  
   if (collapse) ll <- apply(ll, 1, sum)
   
   return(ll)
}

# Structured mixture cumulative density function:
psmix <- function(x, x0, collapse = TRUE, ...){
   if (missing(x0)) x0 <- x$data$x 
   
   # Parse parameter values:
   v <- parse.theta(x$theta)
   
   # Calculate mixture density:
   ll <- matrix(0, nrow = length(x0), ncol = nrow(v))
   
   for (j in 1:nrow(v)) ll[,j] <- v$p[j] * pnorm(x0, v$mu[j], v$sigma[j])
   
   if (collapse) ll <- apply(ll, 1, sum)
   
   return(ll)
}

# Fit Hiatt mixture density to immatures:
lsmix <- function(theta, x, f, discrete = TRUE){
   if ("smix" %in% class(theta)) {
      x <- theta
      theta <- x$theta
   }

   # Convert to 'smix' object:
   x <- smix(x, theta = theta)

   # Calculate log-likelihood:
   if (discrete){
      xm <- xp <- x
      xp$data$x <- xp$data$x + 0.5
      xm$data$x <- xm$data$x - 0.5
      ll <- x$data$f * log(psmix(xp) - psmix(xm))
      ll[x$data$f == 0] <- 0
   }else{
      ll <- x$data$f * log(dsmix(x))
      ll[x$data$f == 0] <- 0
   }

   return(-sum(ll))
}

# Plot model:
plot.smix <- function(x, ...){
   gbarplot(x$data$f, x$data$x, grid = TRUE, col = "grey85", border = "grey65", ylim = c(0, 1.1 * max(x$data$f)))
   x0 <- seq(0, 90, len = 1000)
   
   ll <- dsmix(x, x0, collapse = FALSE)
   for (j in 1:ncol(ll)){
      lines(x0, sum(x$data$f) * ll[,j], lwd = 1, lty = "dashed")
   }
   lines(x0, sum(x$data$f) * apply(ll, 1, sum), lwd = 2, col = "grey30")

   # Display mean sizes:
   mu <- parse.theta(x$theta)$mu
   dy <- approx(x0, sum(x$data$f) * apply(ll, 1, sum), mu)$y
   vline(mu, lower = 0, upper = dy, lwd = 2, col = "grey30")
   text(mu, dy, rownames(parse.theta(x$theta)), pos = 3, font = 2, cex = 1.0)
   
   mtext("Frequency", 2, 2.5, font = 2, cex = 1.25)
   mtext("Carapace width (mm)", 1, 2.5, font = 2, cex = 1.25)

   box(col = "grey50")
}

# Generate random variates from structured mixture model:
rsmix <- function(x, n, ...){
   v <- parse.theta(x$theta)
   N <- rmultinom(1, n, v$p)
   
   r <- NULL
   for (i in 1:length(N)) r <- c(r, rnorm(N[i], v$mu[i], v$sigma[i]))
   
   return(r)
}

theta <- init.smix(k = 6)

theta["log.sigma"] <- -1
theta["log.sigma0"] <- -1
x <- smix(mi, theta = theta)

# Fit proportions:
vars <- names(x$theta)[grep("logit", names(x$theta))]
x$theta[vars] <- optim(x$theta[vars], lsmix, x = x, control = list(trace = 3))$par

# Fit growth parameters:
vars <- c("slope", "intercept", "mu0", "log.sigma0", "log.sigma")
x$theta[vars] <- optim(x$theta[vars], lsmix, x = x, control = list(trace = 3))$par

# Fit complete model:
x$theta <- optim(x$theta, lsmix, x = x, control = list(trace = 3, maxit = 5000))$par

plot.smix(x)
