library(gulf.data)
library(gulf.graphics)

# Read average female density data from the SCS:
x <- read.csv("summer 2024/female SF averages by maturity 1997-2023.csv")

plot(log(x$carapace.width), log(x$immature), type = "l", xlim = c(2, 5))
vline(log(c(10, 14.5, 20.5, 28.5, 37.5, 49)))

# Parse parameter vector:
parse.theta <- function(theta){
   # Mixture proportions:
   logit.p <- theta[grep("logit.p", names(theta))]
   p <- exp(logit.p) / (1 + sum(exp(logit.p)))
   p[length(p)+1] <- 1 - sum(p)
  
   # Number of mixture components:
   k <- length(p)  
   
   # Calculate mean:
   mu <- theta[["mu0"]]
   for (i in 2:k) mu[i] <- theta[["a"]] * mu[i-1] + theta[["b"]]
  
   # Calculate sigma:
   sigma <- exp(theta[["log.sigma0"]])
   for (i in 2:k) sigma[i] <- theta[["a"]] * sigma[i-1] + exp(theta[["log.epsilon"]])
  
   res <- data.frame(p = p, mu = mu, sigma = sigma)
   rownames(res) <- as.roman(4:(4+k-1))
   
   return(res)
}

# Immature mixture likelihood:
loglike <- function(theta, x, f, fixed){
   # Append fixed parameters:
   if (!missing(fixed)) theta <- c(theta, fixed) 
   
   # Parse parameter vector:
   z <- parse.theta(theta)
   
   # Calculate mixture density:
   d <- rep(0, length(x))
   for (i in 1:length(z$p)){
      d <- d + z$p[i] * (pnorm(x+0.5, z$mu[i], z$sigma[i]) - pnorm(x-0.5, z$mu[i], z$sigma[i]))
   }
   
   return(-sum(f[f > 0] * log(d[f > 0])))
}

theta <- c(mu0 = 10, log.sigma0 = 0, 
           a = 1.27, b = 2, log.epsilon = -2,
           logit.p = c(0,0,0,0,0))

loglike(theta, x = x$carapace.width, f = x$immature)
fixed <- theta[-grep("logit.p", names(theta))]
theta <- theta[setdiff(names(theta), names(fixed))]

loglike(theta, x = x$carapace.width, f = x$immature, fixed = fixed)

theta <- optim(theta, fn = loglike, x = x$carapace.width, f = x$immature, fixed = fixed, control = list(trace = 3))$par
theta <- c(theta, fixed)

theta <- optim(theta, fn = loglike, x = x$carapace.width, f = x$immature, fixed = fixed, control = list(maxit = 2000, trace = 3))$par

z <- parse.theta(theta)

# Calculate mixture density:
xx <- seq(0, 95, len = 1000)
dd <- NULL
for (i in 1:k){
  dd <- cbind(dd, z$p[i] * dnorm(xx, z$mu[i], z$sigma[i]))
}
d <- apply(dd, 1, sum)

gbarplot(x$immature, x$carapace.width, xlim = c(0, 70))
gbarplot(x$adolescent, x$carapace.width, col = fade("green"), add = TRUE)
for (i in 1:ncol(dd)){
  lines(xx, sum(x$immature) * dd[,i], col = "red2", lwd = 1, lty = "dashed")
}
lines(xx, sum(x$immature) * d, col = "red2", lwd = 2)
vline(mu, lty = "dashed", col = "red2")

theta.immature <- theta

z <- parse.theta(theta)

mu    <- z["VII","mu"]    # Mean size of instar to be split.
sigma <- z["VII","sigma"] # Standard deviation of instar to be split.

logit.p.split.VII <- 0
logit.t.split.VII <- 0
p.split.VII <- 1 / (1 + exp(-logit.p.split.VII))
t.split.VII <- 1 / (1 + exp(-logit.t.split.VII))

# Split the means of instar VII:
mu.a <- z["VII","mu"] - t.split.VII * sqrt((1-p.split.VII) / p.split.VII) * z["VII","sigma"]  # Mean size of first instar component (immature).
mu.b <- (z["VII","mu"] - p.split.VII * mu.a)/(1-p.split.VII)                                  # Mean size of second instar component (adolescent)..

# Split the variance of instar VII, assuming that the split variances are equal:
V <- z["VII","sigma"]^2 - p.split.VII * (1-p.split.VII) * (mu.b-mu.a)^2
sigma.a <- sigma.b <- sqrt(V)

# Calculate instar VIII values:
mu.a <- theta[["a"]] * mu.a + theta[["b"]]
mu.b <- theta[["a"]] * mu.b + theta[["b"]]
sigma.a <- theta[["a"]] * sigma.a + exp(theta[["log.epsilon"]])
sigma.b <- theta[["a"]] * sigma.b + exp(theta[["log.epsilon"]])

logit.p.split.VIII <- 2
logit.t.split.VIII <- 0
p.split.VIII <- 1 / (1 + exp(-logit.p.split.VIII))
t.split.VIII <- 1 / (1 + exp(-logit.t.split.VIII))

# Split the means of instar VIII:
mu.a <- z["VII","mu"] - t.split.VIII * sqrt((1-p.split.VIII) / p.split.VIII) * z["VII","sigma"]  # Mean size of first instar component (immature).
mu.b <- (z["VII","mu"] - p.split.VIII * mu.a)/(1-p.split.VIII)                                  # Mean size of second instar component (adolescent)..

# Split the variance of instar VII, assuming that the split variances are equal:
V <- z["VII","sigma"]^2 - p.split.VII * (1-p.split.VII) * (mu.b-mu.a)^2
sigma.a <- sigma.b <- sqrt(V)

# Calculate instar VIII values:

mu.a <- theta[["a"]] * mu.a + theta[["b"]]
mu.b <- theta[["a"]] * mu.b + theta[["b"]]
sigma.a <- theta[["a"]] * sigma.a + exp(theta[["log.epsilon"]])
sigma.b <- theta[["a"]] * sigma.b + exp(theta[["log.epsilon"]])

# Immature and adolescent:
loglike <- function(theta, x, fi, fa, fixed){
  # Append fixed parameters:
  if (!missing(fixed)) theta <- c(theta, fixed) 
  
  # Mixture proportions:
  logit.p <- theta[grep("logit.p", names(theta))]
  p <- exp(logit.p) / (1 + sum(exp(logit.p)))
  p[length(p)+1] <- 1 - sum(p)
  
  # Number of mixture components:
  k <- length(p)  
  
  # Calculate mean:
  mu <- theta[["mu0"]]
  for (i in 2:k) mu[i] <- theta[["a"]] * mu[i-1] + theta[["b"]]
  
  # Calculate sigma:
  sigma <- exp(theta[["log.sigma0"]])
  for (i in 2:k) sigma[i] <- theta[["a"]] * sigma[i-1] + exp(theta[["log.epsilon"]])
  
  # Calculate mixture density:
  d <- rep(0, length(x))
  for (i in 1:k){
    d <- d + p[i] * (pnorm(x+0.5, mu[i], sigma[i]) - pnorm(x-0.5, mu[i], sigma[i]))
  }
  
  return(-sum(f[f > 0] * log(d[f > 0])))
}


