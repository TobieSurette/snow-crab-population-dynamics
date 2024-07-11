
# Parse parameter vector:
parse.theta <- function(theta){
  # Calculate mean:
  mu.immature <- theta[["mu0"]]
  for (i in 2:4) mu.immature[i] <- theta[["a"]] * mu.immature[i-1] + theta[["b"]]
  names(mu.immature) <- as.roman(4:7)
  
  # Calculate sigma:
  sigma.immature <- exp(theta[["log.sigma0"]])
  for (i in 2:4) sigma.immature[i] <- theta[["a"]] * sigma.immature[i-1] + exp(theta[["log.epsilon"]])
  names(sigma.immature) <- as.roman(4:7)
  
  # Split instar VII:
  p.split <- 1 / (1 + exp(-theta[grep("logit.p.split.VII", names(theta))]))
  t.split <- 1 / (1 + exp(-theta[grep("logit.t.split.VII", names(theta))]))
  names(p.split) <- as.roman(7:8)
  names(t.split) <- as.roman(7:8)
  
  # Split the means of instar VII:
  mu.a <- mu.immature["VII"] - t.split["VII"] * sqrt((1-p.split["VII"]) / p.split["VII"]) * sigma.immature["VII"] 
  mu.b <- (mu.immature["VII"] - p.split["VII"] * mu.a)/(1-p.split["VII"]) 
  mu.immature["VIII"] <- theta[["a"]] * mu.a + theta[["b"]]
  mu.adolescent <- NA * mu.immature
  mu.adolescent["VIII"] <- theta[["a"]] * mu.b + theta[["b"]]
  
  # Split the variance of instar VII:
  V <- sigma.immature["VII"]^2 - p.split["VII"] * (1-p.split["VII"]) * (mu.b-mu.a)^2
  sigma.a <- sigma.b <- sqrt(V)
  sigma.immature["VIII"] <- theta[["a"]] * sigma.a + exp(theta[["log.epsilon"]])
  sigma.adolescent <- NA * sigma.immature
  sigma.adolescent["VIII"] <- theta[["a"]] * sigma.b + exp(theta[["log.epsilon"]])
    
  # Split the means of instar VIII:
  mu.a <- mu.immature["VIII"] - t.split["VIII"] * sqrt((1-p.split["VIII"]) / p.split["VIII"]) * sigma.immature["VIII"] 
  mu.b <- (mu.immature["VIII"] - p.split["VIII"] * mu.a)/(1-p.split["VIII"]) 
  mu.immature["IX"]   <- theta[["a"]] * mu.a + theta[["b"]]
  mu.adolescent["IX"] <- theta[["a"]] * mu.b + theta[["b"]]
  mu.immature[c("X", "XI")]    <- NA
  mu.adolescent["X"]  <- theta[["a"]] * mu.immature["IX"] + theta[["b"]]
  mu.adolescent["XI"] <- NA
  
  # Split the variance of instar VIII:
  V <- sigma.immature["VIII"]^2 - p.split["VIII"] * (1-p.split["VIII"]) * (mu.b-mu.a)^2
  sigma.a <- sigma.b <- sqrt(V)
  sigma.immature["IX"] <- theta[["a"]] * sigma.a + exp(theta[["log.epsilon"]])
  sigma.adolescent["IX"] <- theta[["a"]] * sigma.b + exp(theta[["log.epsilon"]])
  sigma.immature[c("X", "XI")] <- NA
  sigma.adolescent["X"] <- theta[["a"]] * sigma.immature["IX"] + exp(theta[["log.epsilon"]])
  sigma.adolescent["XI"] <- NA
  
  # Primiparous females:
  mu.primiparous <- NA * mu.adolescent
  mu.primiparous[c("IX", "X", "XI")] <- theta[["a.mature"]] * mu.adolescent[c("VIII", "IX", "X")] + theta[["b.mature"]]
  sigma.primiparous <- NA * sigma.adolescent
  sigma.primiparous[c("IX", "X", "XI")] <- theta[["a"]] * sigma.adolescent[c("VIII", "IX", "X")] + exp(theta[["log.epsilon"]])
  
  # Mortality parameter:
  M <- 1 / (1 + exp(-theta["logit.M"]))
  
  # Population abundances:
  N.immature <- theta[["N.immature0"]]
  for (i in 2:4) N.immature[i] <- (1-M) * N.immature[i-1]
  names(N.immature) <- as.roman(4:7)
  
  # Apply instar maturation split:
  N.immature["VIII"]   <- (1-M) *    p.split["VII"]   * N.immature["VII"]
  N.adolescent <- 0 * N.immature
  N.adolescent["VIII"] <- (1-M) * (1-p.split["VII"])  * N.immature["VII"]
  N.immature["IX"]     <- (1-M) *    p.split["VIII"]  * N.immature["VIII"]
  N.adolescent["IX"]   <- (1-M) * (1-p.split["VIII"]) * N.immature["VIII"]
  N.adolescent["X"]    <- (1-M) * N.immature["IX"]
  N.immature["X"]      <- 0
  N.immature["XI"]     <- 0
  N.adolescent["XI"]   <- 0
  
  # Primiparous:
  N.primiparous <- 0 * N.immature
  N.primiparous["IX"] <- (1-M) * N.adolescent["VIII"]
  N.primiparous["X"]  <- (1-M) * N.adolescent["IX"]
  N.primiparous["XI"] <- (1-M) * N.adolescent["X"]
  
  # Selectivity curve:
  S <- function(x) return(1 / (1 + exp(-4 * exp(theta["log.S.scale"]) * (x-theta["S.location"]))))
  S.immature    <- S(mu.immature)
  S.adolescent  <- S(mu.adolescent)
  S.primiparous <- S(mu.primiparous)
  
  # Observation abundances:
  n.immature   <- S.immature   * N.immature
  n.adolescent <- S.adolescent * N.adolescent
  n.immature[N.immature == 0] <- 0
  n.adolescent[N.adolescent == 0] <- 0
  n.primiparous <- S.primiparous * N.primiparous
  n.primiparous[N.primiparous == 0] <- 0
  
  # Package output:
  res <- data.frame(S.immature        = S.immature,
                    S.adolescent      = S.adolescent,
                    S.primiparous     = S.primiparous,
                    N.immature        = N.immature, 
                    N.adolescent      = N.adolescent, 
                    N.primiparous     = N.primiparous,
                    n.immature        = n.immature, 
                    n.adolescent      = n.adolescent, 
                    n.primiparous     = n.primiparous,
                    mu.immature       = mu.immature, 
                    mu.adolescent     = mu.adolescent,
                    mu.primiparous    = mu.primiparous,
                    sigma.immature    = sigma.immature,
                    sigma.adolescent  = sigma.adolescent,
                    sigma.primiparous = sigma.primiparous)
  rownames(res) <- as.roman(4:(4+nrow(res)-1))
  
  return(res)
}

predict <- function(x, theta){
   # Population model outputs:
   r <- parse.theta(theta)
  
   # Calculate immature mixture density:
   d.immature <- NULL
   for (i in 1:nrow(r)){
      p.immature <- pnorm(x + 0.5, r$mu.immature[i], r$sigma.immature[i]) - pnorm(x - 0.5, r$mu.immature[i], r$sigma.immature[i])
      d.immature <- cbind(d.immature, r$n.immature[i] * p.immature)
   }
   colnames(d.immature) <- rownames(r)
   d.immature[is.na(d.immature)] <- 0
  
   # Calculate adolescent mixture density:
   d.adolescent <- NULL
   for (i in 1:nrow(r)){
      p.adolescent <- pnorm(x + 0.5, r$mu.adolescent[i], r$sigma.adolescent[i]) - pnorm(x - 0.5, r$mu.adolescent[i], r$sigma.adolescent[i])
      d.adolescent <- cbind(d.adolescent, r$n.adolescent[i] * p.adolescent)
   }
   colnames(d.adolescent) <- rownames(r)
   d.adolescent[is.na(d.adolescent)] <- 0
   
   # Calculate primiparous mixture density:
   d.primiparous <- NULL
   for (i in 1:nrow(r)){
     p.primiparous <- pnorm(x + 0.5, r$mu.primiparous[i], r$sigma.primiparous[i]) - pnorm(x - 0.5, r$mu.primiparous[i], r$sigma.primiparous[i])
     d.primiparous <- cbind(d.primiparous, r$n.primiparous[i] * p.primiparous)
   }
   colnames(d.primiparous) <- rownames(r)
   d.primiparous[is.na(d.primiparous)] <- 0
   
   mu.immature    <- apply(d.immature, 1, sum)
   mu.adolescent  <- apply(d.adolescent, 1, sum)
   mu.primiparous <- apply(d.primiparous, 1, sum)
   
   res <- data.frame(mu.immature = mu.immature, 
                     mu.adolescent = mu.adolescent, 
                     mu.primiparous = mu.primiparous)
   
   return(res)
}

# Selectivity curve:
S <- function(x, theta) return(1 / (1 + exp(-4 * exp(theta["log.S.scale"]) * (x-theta["S.location"]))))

# Immature and adolescent:
loglike <- function(theta, x, fi, fa, fp, fixed){
   # Append fixed parameters:
   if (!missing(fixed)) theta <- c(theta, fixed) 
  
   # Population model outputs:
   r <- parse.theta(theta)
  
   # Calculate immature mixture density:
   d.immature <- NULL
   for (i in 1:nrow(r)){
      p.immature <- pnorm(x + 0.5, r$mu.immature[i], r$sigma.immature[i]) - pnorm(x - 0.5, r$mu.immature[i], r$sigma.immature[i])
      d.immature <- cbind(d.immature, r$n.immature[i] * p.immature)
   }
   colnames(d.immature) <- rownames(r)
   d.immature[is.na(d.immature)] <- 0
   
   # Calculate adolescent mixture density:
   d.adolescent <- NULL
   for (i in 1:nrow(r)){
     p.adolescent <- pnorm(x + 0.5, r$mu.adolescent[i], r$sigma.adolescent[i]) - pnorm(x - 0.5, r$mu.adolescent[i], r$sigma.adolescent[i])
     d.adolescent <- cbind(d.adolescent, r$n.adolescent[i] * p.adolescent)
   }
   colnames(d.adolescent) <- rownames(r)
   d.adolescent[is.na(d.adolescent)] <- 0
   
   # Calculate primiparous mixture density:
   d.primiparous <- NULL
   for (i in 1:nrow(r)){
     p.primiparous <- pnorm(x + 0.5, r$mu.primiparous[i], r$sigma.primiparous[i]) - pnorm(x - 0.5, r$mu.primiparous[i], r$sigma.primiparous[i])
     d.primiparous <- cbind(d.primiparous, r$n.primiparous[i] * p.primiparous)
   }
   colnames(d.primiparous) <- rownames(r)
   d.primiparous[is.na(d.primiparous)] <- 0
   
   # Predicted means:
   sigma.obs <- exp(theta[["log.sigma.obs"]])
   mu.immature    <- apply(d.immature, 1, sum)
   mu.adolescent  <- apply(d.adolescent, 1, sum)
   mu.primiparous <- apply(d.primiparous, 1, sum)
   
   # Likelihoods:
   ll.immature    <- dnorm(fi[x <= 80], mu.immature[x <= 80], sigma.obs, log = TRUE)  
   ll.adolescent  <- dnorm(fa, mu.adolescent, sigma.obs, log = TRUE)
   ll.primiparous <- dnorm(fa, mu.primiparous, sigma.obs, log = TRUE)
   
   return(-(sum(ll.immature) + sum(ll.adolescent) + sum(ll.primiparous)))
}

# Read average female density data from the SCS:
x <- read.csv("summer 2024/female SF averages by maturity 1997-2023.csv")

fi <- x$immature
fa <- x$adolescent
fp <- x$primiparous
x <- x$carapace.width
fi[x >= 75] <- 0
fp[x < 22]  <- 0

clg()

theta <- c(N.immature0 = 1809,
           log.S.scale = -2,
           S.location = 23,
           logit.M = -4,
           mu0 = 10.1,  
           log.sigma0 = -0.7109,
           a = 1.28,           
           b = 1.66,
           a.mature = 1.1,
           b.mature = 4,
           log.epsilon = -1.3,
           logit.p.split.VII = 1.3,
           logit.t.split.VII = 0.2,
           logit.p.split.VIII = -2,
           logit.t.split.VIII = -2.5,
           log.sigma.obs = 3.2)

xx <- seq(0, 100, len = 1000)
p <- predict(xx, theta)
r <- parse.theta(theta)

plot(x, fi, type= "l", col = "green3", lwd = 2, lty = "dashed", xlim = c(0, 80), ylim = c(0, 250), xaxs = "i", yaxs = "i", xlab = "", ylab = "")
grid()
lines(xx, p$mu.immature, col = "green3", lwd = 2)
lines(x, fa, col = "yellow3", lwd = 2, lty = "dashed")
lines(xx, p$mu.adolescent, col = "yellow3", lwd = 2)
lines(x, fp, col = "red2", lwd = 2, lty = "dashed")
lines(xx, p$mu.primiparous, col = "red2", lwd = 2)
vline(r$mu.immature, lty = "dashed", col = "green3")
vline(r$mu.adolescent, lty = "dashed", col = "yellow3")
vline(r$mu.primiparous, lty = "dashed", col = "red2")
mtext("Carapace width (mm)", 1, 2.5, font = 2, cex = 1.25)
mtext("Abundance", 2, 2.5, font = 2, cex = 1.25)
box(col = "grey50")


log.S.scale <- 2.5
S.location  <- 18
S.offset    <- 15
S.theta <- (exp(log.S.scale)^2) / S.location
S.k     <- S.location / S.theta

plot(xx, pgamma(xx-S.offset, shape = S.k, scale = S.theta), type = "l", lwd = 2, col = "blue")
grid()
vline(S.offset + S.location)

parse.theta(theta)

loglike(theta, x = x, fi = fi, fa = fa, fp = fp)
vars <- c("N.immature0", "S.location", "log.sigma.obs")
fixed <- theta[setdiff(names(theta), vars)]
theta <- theta[setdiff(names(theta), names(fixed))]

loglike(theta, x = x, fi = fi, fa = fa, fixed = fixed)
theta <- optim(theta, fn = loglike, x = x, fi = fi, fa = fa, fixed = fixed, control = list(maxit = 2000, trace = 3))$par
parse.theta(c(fixed, theta))
theta <- c(theta, fixed)

vars <- c("N.immature0", "S.location", "log.sigma.obs", "logit.M")
fixed <- theta[setdiff(names(theta), vars)]
theta <- theta[setdiff(names(theta), names(fixed))]
loglike(theta, x = x, fi = fi, fa = fa, fixed = fixed)
theta <- optim(theta, fn = loglike, x = x, fi = fi, fa = fa, fixed = fixed, control = list(maxit = 2000, trace = 3))$par
theta <- c(theta, fixed)



plot(x, S(x, c(fixed, theta)), type = "l")

p <- predict(x, c(fixed, theta))

plot(x, fi, type= "l", col = "green3", lwd = 2, lty = "dashed", xlim = c(0, 80), ylim = c(0, 250), xaxs = "i", yaxs = "i", xlab = "", ylab = "")
grid()
lines(x, p$mu.immature, col = "green3", lwd = 2)
lines(x, fa, col = "yellow3", lwd = 2, lty = "dashed")
lines(x, p$mu.adolescent, col = "yellow3", lwd = 2)
lines(x, fp, col = "red3", lwd = 2, lty = "dashed")
lines(x, p$mu.primiparous, col = "red3", lwd = 2)
mtext("Carapace width (mm)", 1, 2.5, font = 2, cex = 1.25)
mtext("Abundance", 2, 2.5, font = 2, cex = 1.25)
box(col = "grey50")

theta <- c(theta, fixed)

vars <- c("N.immature0", "S.location", "log.S.scale", "log.sigma.obs", "log.epsilon")
fixed <- theta[setdiff(names(theta), vars)]
theta <- theta[setdiff(names(theta), names(fixed))]

loglike(theta, x = x, fi = fi, fa = fa, fixed = fixed)

theta <- optim(theta, fn = loglike, x = x, fi = fi, fa = fa, fixed = fixed, control = list(maxit = 2000, trace = 3))$par

plot(x, S(x, c(fixed, theta)), type = "l")

p <- predict(x, c(fixed, theta))

plot(x, fi, type= "l", col = "red")
lines(x, p$mu.immature, col = "blue")

plot(x, fa, type= "l", col = "red")
lines(x, p$mu.adolescent, col = "blue")


