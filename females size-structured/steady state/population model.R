# Steady-state population model:

loglike <- function(phi, x, fixed){
   # Base mixture parameters:
   v <- parse.theta(x$theta)
   
   if (!missing(fixed)) phi <- c(phi, fixed)
   
   # Selectivity:
   sel <- function(x, mu = 30, sigma = 10){
      shape <- mu^2 / sigma^2
      scale <- sigma^2 / mu
      v <- pgamma(x, shape = shape, scale = scale)
      return(v)
   }
   S <- sel(v$mu, mu = exp(phi["log.sel.mu"]), sigma = exp(phi["log.sel.sigma"])) 
   names(S) <- rownames(v)
  
   # Mortality:
   M <- rep(exp(phi["log.M"]), nrow(v))  # Probability of dying immature.
   names(M) <- rownames(v)
  
   # Maturation:
   pi <- rep(1, nrow(v))  # Probability of remaining immature.
   names(pi) <- rownames(v)
   logistic <- function(x) return(1 / (1 + exp(-x)))
   pi["VII"]  <- logistic(phi["logit.mat.VII"])
   pi["VIII"] <- logistic(phi["logit.mat.VIII"])
   pi["IX"] <- 0
   pi["X"] <- 0
   
   #pp <- (1 / S) * (1-M) * (1-pi) * v$p # Pubescent transition.
   #names(pp) <- as.roman(rownames(v)) + 1
   
   # Recruitment:
   R <- c(exp(phi["log.R"]))
   names(R) <- rownames(v)[1]
   
   # Calculate new population distribution for immatures:
   tmp <- (1 / S) * (1-M) * pi * v$p    # Apply reverse selectivity, mortality and molt to immaturity.
   tmp <- tmp[1:(length(tmp)-1)]        # Assume last instar matures to next phase.
   names(tmp) <- rownames(v)[2:nrow(v)] # Increment instar numbers.
   pnew <- S * c(R, tmp)                # Add recruitment and apply selectivity.
   
   # Calculate mixture density:
   ll <- matrix(0, nrow = length(x$data$x), ncol = nrow(v))
   for (j in 1:nrow(v)) ll[,j] <- pnew[j] * (pnorm(x$data$x+0.5, v$mu[j], v$sigma[j]) - pnorm(x$data$x-0.5, v$mu[j], v$sigma[j]))
   ll <- apply(ll, 1, sum)
  
   # Sum-of-squares:
   ss <- sum((sum(x$data$f) * ll - x$data$f)^2)
     
   return(ss)
}

# Initialize population model parameters:
phi <- c(log.sel.mu = log(25),    # Selectivity center size (log-scale).
         log.sel.sigma = log(3), # Selectivity error size (log-scale). 
         log.M = log(0.1),        # Mortality proportion (log-scale).
         logit.mat.VII = 0,       # Probability of instar VII remaining immature (logit-scale).
         logit.mat.VIII = -3,     # Probability of instar VIII remaining immature (logit-scale).
         log.R = log(1))        # Recruitment of first instar (log-scale).

v <- parse.theta(x$theta)

# Selectivity:
sel <- function(x, mu = 30, sigma = 10){
  shape <- mu^2 / sigma^2
  scale <- sigma^2 / mu
  v <- pgamma(x, shape = shape, scale = scale)
  return(v)
}
#sel <- function(x, mu = 30, sigma = 10){
#  v <- 1 / (1 + exp(-(x-mu)/sigma))
#  return(v)
#}
S <- sel(v$mu, mu = exp(phi["log.sel.mu"]), sigma = exp(phi["log.sel.sigma"])) 
names(S) <- rownames(v)

# Mortality:
M <- rep(exp(phi["log.M"]), nrow(v))  # Probability of dying immature.
names(M) <- rownames(v)

# Maturation:
pi <- rep(1, nrow(v))  # Probability of remaining immature.
names(pi) <- rownames(v)
logistic <- function(x) return(1 / (1 + exp(-x)))
pi["VII"]  <- logistic(phi["logit.mat.VII"])
pi["VIII"] <- logistic(phi["logit.mat.VIII"])
pi["IX"] <- 0
pi["X"] <- 0

# Recruitment:
R <- c(exp(phi["log.R"]))
names(R) <- rownames(v)[1]

# Calculate new population distribution:
tmp <- (1 / S) * (1-M) * pi * v$p   # Apply reverse selectivity, mortality and molt to immaturity.
tmp <- tmp[1:(length(tmp)-1)] # Assume last instar matures to next phase.
names(tmp) <- rownames(v)[2:nrow(v)]  # Increment instar numbers.
pnew <- S * c(R, tmp)          # Add recruitment and apply selectivity.

# Calculate mixture density:
ll <- matrix(0, nrow = length(x$data$x), ncol = nrow(v))
for (j in 1:nrow(v)) ll[,j] <- pnew[j] * (pnorm(x$data$x+0.5, v$mu[j], v$sigma[j]) - pnorm(x$data$x-0.5, v$mu[j], v$sigma[j]))
ll <- apply(ll, 1, sum)

gbarplot(x$data$f, x$data$x)
lines(x$data$x, sum(x$data$f) * ll)

plot(x0, sel(x0, mu = exp(phi["log.sel.mu"]), sigma = exp(phi["log.sel.sigma"]))) 


# Initialize population model parameters:
phi <- c(log.sel.mu = log(25),    # Selectivity center size (log-scale).
         log.sel.sigma = log(10), # Selectivity error size (log-scale). 
         log.M = log(0.1),        # Mortality proportion (log-scale).
         logit.mat.VII = 0,       # Probability of instar VII remaining immature (logit-scale).
         logit.mat.VIII = -3,     # Probability of instar VIII remaining immature (logit-scale).
         log.R = log(1))          # Recruitment of first instar (log-scale).


fixed <- c("log.sel.mu", "log.sel.sigma", "log.M", "logit.mat.VII", "logit.mat.VIII")
fixed <- c("log.M", "logit.mat.VII", "logit.mat.VIII")
fixed <- phi[fixed]
phi <- phi[setdiff(names(phi), names(fixed))]

loglike(phi, x, fixed = fixed)
  
phi <- optim(phi, loglike, x = x, fixed = fixed, control = list(trace = 3))$par

phi <- c(phi, fixed)
fixed <- c("log.M", "logit.mat.VII", "logit.mat.VIII")
fixed <- phi[fixed]
phi <- phi[setdiff(names(phi), names(fixed))]


phi <- optim(phi, loglike, x = x, fixed = fixed, control = list(trace = 3))$par
phi <- c(phi, fixed)

phi <- optim(phi, loglike, x = x, control = list(trace = 3, maxit = 5000))$par

