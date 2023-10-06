# Fit Hiatt mixture density to immatures:
llmix <- function(theta, x, fixed){
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

theta <- c(logit.p = c(1,1,1,1,1), mu0 = 10, intercept = 2, slope = 1.27, log.sigma = -2)
fixed <- theta[-grep("logit", names(theta))]
theta <- theta[setdiff(names(theta), names(fixed))]

llmix(theta, mi, fixed)
theta <- optim(theta, llmix, x = mi, fixed = fixed, control = list(trace = 3))$par

theta <- c(theta, fixed["log.sigma"])
fixed <- fixed[-grep("log.sigma", names(fixed))]
theta <- optim(theta, llmix, x = mi, fixed = fixed, control = list(trace = 3))$par

theta <- c(theta, fixed)
theta <- optim(theta, llmix, x = mi, control = list(trace = 3, maxit = 2000))$par


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

# Data example:
png(file = "C:/Users/SuretteTJ/Desktop/gulf-population-modelling/studies/females size-based model/figures/female instar example.png",
    units = "in", res = 300, height = 6, width = 7.5)

gbarplot(mi, xlim = c(0, 65), ylim = c(0, 225), xaxs = "i", col = "grey80", border = "grey70", lwd = 0.5, grid = TRUE)
x0 <- seq(0, 100, len = 1000)
ll <- rep(0, length(x0))
for (i in 1:k){
   pp <- p[i] * dnorm(x0, mu[i], sigma[i])
   ll <- ll + pp
   #lines(x0, sum(mi) * pp, lwd = 1.5, col = "grey20")
   text(mu[i], sum(mi) * pp[which.max(pp)], names(mu)[i], pos = 3, font = 2, cex = 1.1)
   
}
mtext("Carapace width (mm)", 1, 2.25, cex = 1.25, font = 2)
mtext("Frequency", 2, 2.25, cex = 1.25, font = 2)
box(col = "grey50")
dev.off()

# Plot mixture fit:
png(file = "C:/Users/SuretteTJ/Desktop/gulf-population-modelling/studies/females size-based model/figures/female instar mixture example.png",
    units = "in", res = 300, height = 6, width = 7.5)

gbarplot(mi, xlim = c(0, 65), ylim = c(0, 225), xaxs = "i", col = "grey80", border = "grey70", lwd = 0.5, grid = TRUE)
x0 <- seq(0, 100, len = 1000)
ll <- rep(0, length(x0))
for (i in 1:k){
   pp <- p[i] * dnorm(x0, mu[i], sigma[i])
   ll <- ll + pp
   lines(x0, sum(mi) * pp, lwd = 1.5, col = "grey20")
   text(mu[i], sum(mi) * pp[which.max(pp)], names(mu)[i], pos = 3, font = 2, cex = 1.1)
}
#lines(x0, sum(mi) * ll, lwd = 2)
mtext("Carapace width (mm)", 1, 2.25, cex = 1.25, font = 2)
mtext("Frequency", 2, 2.25, cex = 1.25, font = 2)
box(col = "grey50")
dev.off()

# Simulate:
n <- 100000
np <- rmultinom(1, size = n, prob = p)[,1]
tmp <- NULL
y <- NULL
for (i in 1:length(np)){
   y <- c(y, rnorm(np[i], mu[i], sigma[i]))
   tmp <- c(tmp, rep(i, each = np[i]))
}
r <- data.frame(carapace.width = round(y,1), instar = tmp)
gbarplot(table(round(r$carapace.width)))


# Define growth transition matrices:
Gi <- growth.matrix(as.numeric(fvars), theta = c(intercept = 2.276, slope = 0.24, log.sigma = -2.5))   
Gp <- growth.matrix(as.numeric(fvars), theta = c(intercept = 7.7, slope = 0.126, log.sigma = -2))
Gm <- growth.matrix(as.numeric(fvars), theta = c(intercept = 5, slope = 0.11, log.sigma = -2))

# Define instar-stage transition probabilities:
p.I.to.I <- c(1, 1, 1, 0.43, 0.03, 0); names(p.I.to.I) <- names(mu)
p.I.to.P <- 1 - p.I.to.I; names(p.I.to.P) <- names(mu)
p.P.to.M <- c(1, 1, 1, 1, 1, 1); names(p.P.to.M) <- names(mu)

# Calculate discrete probabilities from fitted mixture:
pi <- rep(0, length(fvars))
names(pi) <- fvars
pp <- pi
for (i in 1:k){
   for (j in 1:length(pp)-1){
      pi[j] <- pi[j] + p.I.to.I[i] * p[i] * (pnorm(as.numeric(fvars[j])+0.5, mu[i], sigma[i]) - pnorm(as.numeric(fvars[j])-0.5, mu[i], sigma[i]))
      pp[j] <- pp[j] + p.I.to.P[i] * p[i] * (pnorm(as.numeric(fvars[j])+0.5, mu[i], sigma[i]) - pnorm(as.numeric(fvars[j])-0.5, mu[i], sigma[i]))
   }
}
pi <- (pi %*% Gi)[1, names(pi)]
pp <- (pp %*% Gp)[1, names(pp)]
pm <- (pp %*% Gm)[1, names(pp)]

png(file = "C:/Users/SuretteTJ/Desktop/gulf-population-modelling/studies/females size-based model/figures/female instar mixture.png",
    units = "in", res = 300, height = 6, width = 7.5)

plot(as.numeric(fvars), mi, type = "n", 
     xlim = c(0, 90), ylim = c(0, 250), xaxs = "i", yaxs = "i", xlab = "", ylab = "")
grid()
lines(as.numeric(fvars), mi, lwd = 2, col = cols[1])
lines(as.numeric(fvars), mp, lwd = 2, col = cols[2])
lines(as.numeric(fvars), mm, lwd = 2, col = cols[3])

lines(as.numeric(names(pi)), (1 / sum(pi)) * sum(mi) * pi, lwd = 2, lty = "dashed", col = cols[1])
lines(as.numeric(names(pp)), (1 / sum(pp)) * sum(mp) * pp, lwd = 2, lty = "dashed", col = cols[2])
lines(as.numeric(names(pm)), (1 / sum(pm)) * sum(mm) * pm, lwd = 2, lty = "dashed", col = cols[3])

for (i in 1:k){
   ppp <- p[i] * dnorm(x0, mu[i], sigma[i])
   text(mu[i], sum(mi) * ppp[which.max(ppp)], names(mu)[i], pos = 3, font = 2, cex = 1.1, col = cols[1])
}


legend("topright", 
       legend = c("Immature", "Pubescent", "Primiparous"),
       lwd = 2, col = cols, cex = 1.0)

mtext("Carapace width (mm)", 1, 2.25, cex = 1.25, font = 2)
mtext("Frequency", 2, 2.25, cex = 1.25, font = 2)

box()

dev.off()
