mu <- 28.4   # Mean size of instar to be split.
sigma <- sqrt(8.786695)  # Standard deviation of instar to be split.

x <- seq(0, 100, len = 1000)

p <- 0.87943074 #0.9724  # Proportion that remains immature.
t <- 0.289323455 #0.1687  # Instar separation parameter 0 =< t < 1.

#mu <- mu.a + t * sqrt((1-p)/p) * sigma
mu.a <- mu - t * sqrt((1-p)/p) * sigma   # Mean size of first instar component (immature).
mu.b <- (mu - p * mu.a)/(1-p)            # Mean size of second instar component (adolescent)..

# Calculate instar variances, assuming that they are equal:
V <- sigma^2 - p * (1-p) * (mu.b-mu.a)^2
sigma.a <- sigma.b <- sqrt(V)

plot(x, dnorm(x, mu, sigma), type = "l", yaxs = "i", 
     xlim = c(18, 37), ylim = c(0, 0.25), lwd = 2, col = "blue",
     xlab = "", ylab = "")
grid()
#vline(mu, lty = "dashed", lwd = 2, col = "blue")

lines(x, p * dnorm(x, mu.a, sigma.a), lty = "dashed", col = "green3", lwd = 1) 
lines(x, (1-p) * dnorm(x, mu.b, sigma.b), lty = "dashed", col = "yellow3", lwd = 1) 
abline(v=c(mu.a, mu, mu.b), lty = "dashed", lwd = 1)

lines(x, p * dnorm(x, mu.a, sigma.a) + (1-p) * dnorm(x, mu.b, sigma.b), lty = "solid", col = "red3", lwd = 2)

p * sigma.a^2 + (1-p) * sigma.b^2 + p * (1-p) * (mu.b-mu.a)^2

mtext(round(mu.a, 1), 3, at = mu.a, font = 2, cex = 0.75)
mtext(round(mu.b, 1), 3, at = mu.b, font = 2, cex = 0.75)
mtext("Carapace width (mm)", 1, 2.5, font = 2, cex = 1.25)
mtext("Density", 2, 2.5, font = 2, cex = 1.25)

legend("topleft", 
       legend = c("Original instar", "Split instar", "Split immature", "Split adolescent"), 
       lwd = c(2,2,1),
       lty = c("solid", "solid", "dashed", "dashed"),
       col = c("blue", "red3", "green3", "yellow3"))

box(col = "grey50")

#plot(x, (1-p) * dnorm(x, mu.b, sigma.b) / (p * dnorm(x, mu.a, sigma.a) + (1-p) * dnorm(x, mu.b, sigma.b)), type = "l")