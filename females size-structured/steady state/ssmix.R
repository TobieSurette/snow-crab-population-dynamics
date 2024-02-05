# Steady-state mixture model:

# Initialize parameters:
init.ssmix <- function(x, ...){
   theta <- c(slope.immature = 1.23,             # Immature growth slope parameter.
              intercept.immature = 2.76,         # Immature growth intercept parameter.
              slope.pubescent = 1.25,            # Pubescent growth slope parameter.
              intercept.pubescent = 2.7,         # Pubescent growth intercept parameter.
              slope.mature = 1.25,               # Mature growth slope parameter.
              intercept.mature = 2.7,            # Mature growth intercept parameter.   
              mu0 = 9.6,                         # Size of first immature instar.
              log.sigma0 = -0.62,                # Size-variation of first immature instar.
              log.sigma = 0.19,                  # Extra instar variation growth.
              log.scale.immature = c(1.3, 4.5, 6.3, 7.2, 7.3, 4.1)  # Immature instar abundance scales.
              )
   
   # Initialize population model parameters:
   phi <- c(log.sel.mu = log(3.41),              # Selectivity center size (log-scale).
            log.sel.sigma = log(2.25),           # Selectivity error size (log-scale). 
            log.M = log(0.2),                    # Mortality proportion (log-scale).
            logit.mat.VII = 0,                   # Probability of instar VII remaining immature (logit-scale).
            logit.mat.VIII = -3,                 # Probability of instar VIII remaining immature (logit-scale).
            log.R = -3.0)                        # Recruitment of first instar (log-scale).
   
   return(c(theta, phi))
}

# Parse structured mixture parameter vector:
parse.parameters <- function(theta){
   # Instar abundances:
   p.immature  <- exp(theta[grep("log.scale.immature", names(theta))])

   # Number of immature instars:
   k <- length(p.immature)
  
   # Name instar stages:
   names(p.immature)  <- as.roman(4:(k+3))

   # Instar mean sizes:
   mu.immature <- c(theta["mu0"], rep(NA, k-1))
   names(mu.immature) <- names(p.immature)
   for (i in 2:k) mu.immature[i] <- theta["intercept.immature"] + theta["slope.immature"] * mu.immature[i-1]
   mu.pubescent <- theta["intercept.pubescent"] + theta["slope.pubescent"] * mu.immature[(k-2):k]
   mu.mature    <- theta["intercept.mature"] + theta["slope.mature"] * mu.pubescent 

   names(mu.immature)  <- names(p.immature)
   names(mu.pubescent) <- names(p.pubescent)
   names(mu.mature)    <- names(p.mature)
  
   # Instar standard errors:
   sigma.immature <- exp(theta["log.sigma0"])
   for (i in 2:k) sigma.immature[i] <- sqrt((theta["slope.immature"] * sigma.immature[i-1])^2 + exp(theta["log.sigma"])^2)
   names(sigma.immature)  <- names(p.immature)
   sigma.pubescent <- sqrt((theta["slope.pubescent"] * sigma.immature[(k-2):k])^2 + exp(theta["log.sigma"])^2)
   sigma.mature    <- sqrt((theta["slope.mature"] * sigma.pubescent)^2 + exp(theta["log.sigma"])^2)
   names(sigma.immature)  <- names(p.immature)
   names(sigma.pubescent) <- names(p.pubescent)
   names(sigma.mature)    <- names(p.mature)
   
   #============================================================================
   # Selectivity:
   sel <- function(x, mu = 30, sigma = 10){
     shape <- mu^2 / sigma^2
     scale <- sigma^2 / mu
     v <- pgamma(x, shape = shape, scale = scale)
     return(v)
   }
   S.immature  <- sel(mu.immature,  mu = exp(theta["log.sel.mu"]), sigma = exp(theta["log.sel.sigma"])) 
   S.pubescent <- sel(mu.pubescent, mu = exp(theta["log.sel.mu"]), sigma = exp(theta["log.sel.sigma"])) 
   S.mature    <- sel(mu.mature,    mu = exp(theta["log.sel.mu"]), sigma = exp(theta["log.sel.sigma"])) 
   
   # Mortality:
   M.immature  <- rep(exp(theta["log.M"]), length(mu.immature))   
   M.pubescent <- rep(exp(theta["log.M"]), length(mu.pubescent))  
   M.mature    <- rep(exp(theta["log.M"]), length(mu.mature))     
   names(M.immature)  <- names(mu.immature)
   names(M.pubescent) <- names(mu.pubescent)
   names(M.mature)    <- names(mu.mature)
   
   # Probability of remaining immature:
   p.moult.imm <- rep(1, length(mu.immature))  
   names(p.moult.imm) <- names(mu.immature)
   p.moult.imm["VII"]  <- 1 / (1 + exp(-theta["logit.mat.VII"])) 
   p.moult.imm["VIII"] <- 1 / (1 + exp(-theta["logit.mat.VIII"])) 
   p.moult.imm["IX"] <- 0
   
   # Calculate new population distribution for immatures:
   p.immature.new <- (1 / S.immature) * (1-M.immature) * p.moult.imm * p.immature
   p.immature.new <- p.immature.new[1:(k-1)]
   names(p.immature.new) <- as.roman(names(p.immature.new)) + 1                     
   p.immature.new <- c(IV = exp(theta[["log.R"]]), p.immature.new)
   p.immature.new <- S.immature * p.immature.new
   
   # Calculate new population distribution for pubescents:
   p.pubescent <- (1 / S.immature) * (1-M.immature) * (1-p.moult.imm) * p.immature
   p.pubescent <- p.pubescent[(k-2):k]
   names(p.pubescent) <- as.roman(names(p.pubescent)) + 1                     
   p.pubescent <- S.pubescent * p.pubescent
   
   # Calculate new population distribution for pubescents:
   p.mature <- (1 / S.pubescent) * (1-M.pubescent) * p.pubescent
   names(p.mature) <- as.roman(names(p.mature)) + 1                     
   p.mature <- S.mature * p.mature
   
   # Create list of mixture proportions and moments:
   v.immature  <- data.frame(stage = "immature",  instar = names(mu.immature), scale = p.immature.new, mu = mu.immature, sigma = sigma.immature)
   v.pubescent <- data.frame(stage = "pubescent", instar = names(mu.pubescent), scale = p.pubescent, mu = mu.pubescent, sigma = sigma.pubescent)
   v.mature    <- data.frame(stage = "mature", instar = names(mu.mature), scale = p.mature, mu = mu.mature, sigma = sigma.mature)
   v <- rbind(v.immature, v.pubescent, v.mature)
   rownames(v) <- NULL

   return(v)  
}

# Population likelihood function:
loglike.ssmix <- function(theta, x){
   # Add fixed parameters:
   theta <- c(theta, x$theta[setdiff(names(x$theta), names(theta))])
   
   v <- parse.parameters(theta)
   
   # Immature mean:
   ix <- which(v$stage == "immature")
   ll.immature <- matrix(0, nrow = nrow(x$data), ncol = length(ix))
   for (i in 1:length(ix)){
      ll.immature[,i] <- v$scale[ix[i]] * (pnorm(x$data$cw + 0.5, v$mu[ix[i]], v$sigma[ix[i]]) - 
                                           pnorm(x$data$cw - 0.5, v$mu[ix[i]], v$sigma[ix[i]]))
   } 
   
   # Pubescent mean:
   ix <- which(v$stage == "pubescent")
   ll.pubescent <- matrix(0, nrow = nrow(x$data), ncol = length(ix))
   for (i in 1:length(ix)){
     ll.pubescent[,i] <- v$scale[ix[i]] * (pnorm(x$data$cw + 0.5, v$mu[ix[i]], v$sigma[ix[i]]) - 
                                           pnorm(x$data$cw - 0.5, v$mu[ix[i]], v$sigma[ix[i]]))
   } 
   
   # Pubescent mean:
   ix <- which(v$stage == "mature")
   ll.mature <- matrix(0, nrow = nrow(x$data), ncol = length(ix))
   for (i in 1:length(ix)){
     ll.mature[,i] <- v$scale[ix[i]] * (pnorm(x$data$cw + 0.5, v$mu[ix[i]], v$sigma[ix[i]]) - 
                                        pnorm(x$data$cw - 0.5, v$mu[ix[i]], v$sigma[ix[i]]))
   } 
   
   # Sum-of-squares:
   ss <- sum((x$data$immature - apply(ll.immature, 1, sum))^2) + 
         sum((x$data$pubescent - apply(ll.pubescent, 1, sum))^2) + 
         sum((x$data$mature - apply(ll.mature, 1, sum))^2)
   
   return(ss)
}

x <- list()
x$data <- data.frame(cw = as.numeric(fvars), immature = mi, pubescent = mp, mature = mn)
x$theta <- init.ssmix()
rownames(x$data) <- NULL

theta <- x$theta[c("log.R", "log.sel.mu", "log.sel.sigma")]
loglike.ssmix(theta, x)
theta <- optim(theta, loglike.ssmix, x = x, control = list(trace = 3))$par
x$theta <- c(theta, x$theta[setdiff(names(x$theta), names(theta))])


theta <- x$theta[c("logit.mat.VII", "logit.mat.VIII")]
loglike.ssmix(theta, x)
theta <- optim(theta, loglike.ssmix, x = x, control = list(trace = 3))$par
x$theta <- c(theta, x$theta[setdiff(names(x$theta), names(theta))])

vars <- c("log.R", "log.sel.mu", "log.sel.sigma", "logit.mat.VII", "logit.mat.VIII")
theta <- x$theta[vars]
loglike.ssmix(theta, x)
theta <- optim(theta, loglike.ssmix, x = x, control = list(trace = 3))$par
x$theta <- c(theta, x$theta[setdiff(names(x$theta), names(theta))])

vars <- c("slope.pubescent", "intercept.pubescent", "slope.mature", "intercept.mature")
theta <- x$theta[vars]
loglike.ssmix(theta, x)
theta <- optim(theta, loglike.ssmix, x = x, control = list(trace = 3))$par
x$theta <- c(theta, x$theta[setdiff(names(x$theta), names(theta))])

vars <- c("log.R", "log.sel.mu", "log.sel.sigma", "logit.mat.VII", "logit.mat.VIII", "slope.pubescent", "intercept.pubescent", "slope.mature", "intercept.mature")
theta <- x$theta[vars]
loglike.ssmix(theta, x)
theta <- optim(theta, loglike.ssmix, x = x, control = list(trace = 3))$par
x$theta <- c(theta, x$theta[setdiff(names(x$theta), names(theta))])

vars <- c("log.R", "log.sel.mu", "log.sel.sigma", "logit.mat.VII", "logit.mat.VIII", 
          "slope.pubescent", "intercept.pubescent", "slope.mature", "intercept.mature",
          "slope.immature", "intercept.immature", "log.sigma")
theta <- x$theta[vars]
loglike.ssmix(theta, x)
theta <- optim(theta, loglike.ssmix, x = x, control = list(maxit = 1500, trace = 3))$par
x$theta <- c(theta, x$theta[setdiff(names(x$theta), names(theta))])

# Add mortality:
vars <- c("log.M", "log.R", "log.sel.mu", "log.sel.sigma", "logit.mat.VII", "logit.mat.VIII", 
          "slope.pubescent", "intercept.pubescent", "slope.mature", "intercept.mature",
          "slope.immature", "intercept.immature", "log.sigma")
theta <- x$theta[vars]
loglike.ssmix(theta, x)
for (i in 1:5) theta <- optim(theta, loglike.ssmix, x = x, control = list(maxit = 2500, trace = 3))$par
x$theta <- c(theta, x$theta[setdiff(names(x$theta), names(theta))])


loglike.ssmix(x$theta, x)

v <- parse.parameters(x$theta)

# Immature mean:
ix <- which(v$stage == "immature")
ll.immature <- matrix(0, nrow = nrow(x$data), ncol = length(ix))
for (i in 1:length(ix)) ll.immature[,i] <- v$scale[ix[i]] * (pnorm(x$data$cw + 0.5, v$mu[ix[i]], v$sigma[ix[i]]) - pnorm(x$data$cw - 0.5, v$mu[ix[i]], v$sigma[ix[i]]))

# Pubescent mean:
ix <- which(v$stage == "pubescent")
ll.pubescent <- matrix(0, nrow = nrow(x$data), ncol = length(ix))
for (i in 1:length(ix)) ll.pubescent[,i] <- v$scale[ix[i]] * (pnorm(x$data$cw + 0.5, v$mu[ix[i]], v$sigma[ix[i]]) - pnorm(x$data$cw - 0.5, v$mu[ix[i]], v$sigma[ix[i]]))

# Pubescent mean:
ix <- which(v$stage == "mature")
ll.mature <- matrix(0, nrow = nrow(x$data), ncol = length(ix))
for (i in 1:length(ix)) ll.mature[,i] <- v$scale[ix[i]] * (pnorm(x$data$cw + 0.5, v$mu[ix[i]], v$sigma[ix[i]]) - pnorm(x$data$cw - 0.5, v$mu[ix[i]], v$sigma[ix[i]]))

cols <- c("deepskyblue", "green3", "royalblue", "brown3")
plot(c(0, 85), c(0, 250), xaxs = "i", yaxs = "i", type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
grid()
lines(x$data$cw, apply(ll.immature, 1, sum), lwd = 2, col = cols[1])
points(x$data$cw, x$data$immature, pch = 21, bg = cols[1])

lines(x$data$cw, apply(ll.pubescent, 1, sum), lwd = 2, col = cols[2])
points(x$data$cw, x$data$pubescent, pch = 21, bg = cols[2])

lines(x$data$cw, apply(ll.mature, 1, sum), lwd = 2, col = cols[3])
points(x$data$cw, x$data$mature, pch = 21, bg = cols[3])

#points(x$data$cw, mo / 3, pch = 21, bg = cols[4])

axis(1)

plot(x0, sel(x0, exp(x$theta["log.sel.mu"]), exp(x$theta["log.sel.sigma"])))

