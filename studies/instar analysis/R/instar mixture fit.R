library(gulf.data)
library(gulf.graphics)

x <- read.scsbio(1998:2022)

# Prepare data:
xi <- x$carapace.width[which(!is.na(x$carapace.width) & (x$sex == 2) & (x$carapace.width <= 70) & !is.mature(x) & !is.pubescent.scsbio(x))]
xp <- x$carapace.width[which(!is.na(x$carapace.width) & (x$carapace.width >= 28) & (x$carapace.width <= 80) & (x$sex == 2) & !is.mature(x) & is.pubescent.scsbio(x))]
xm <- x$carapace.width[which(!is.na(x$carapace.width) & (x$carapace.width >= 36) & (x$carapace.width <= 95) & (x$sex == 2) & is.primiparous.scsbio(x))]  

ti <- table(round(log(xi), 2))
tp <- table(round(log(xp), 2))
tm <- table(round(log(xm), 2))

theta <- c(mu4 = 2.2957724, 
           log.sigma = -2.5701350, 
           log.sigma.pubescent = -1,
           log.sigma.mature = -1,
           logit.p.immature4 = -1.9025520, 
           logit.p.immature5 = 0.3302821, 
           logit.p.immature6 = 2.1653939, 
           logit.p.immature7 = 3.3263875, 
           logit.p.immature8 = 2.9821234, 
           logit.p.pubescent8 = 3.3735471, 
           logit.p.pubescent9 = 2.9042755,
           logit.p.mature9    = 8.2518782, 
           logit.p.mature10    = 8.2026561,
           a = 0.8667117, b = 0.7225728,            
           a.pubescent = 0.7069271, b.pubescent = 1.3312277,          
           a.mature = 0.934210, b.mature = 0.4103975)
         
fixed <- c("a", "b", "mu4", "log.sigma", "log.sigma.pubescent", "log.sigma.mature", "a.pubescent", "b.pubescent", "a.mature", "b.mature")
fixed <- theta[fixed]
theta <- theta[setdiff(names(theta), names(fixed))]

loglike.instar(theta, xi = ti, fixed = fixed)
loglike.instar(theta, xi = ti, xp = tp, fixed = fixed)
loglike.instar(theta, xi = ti, xp = tp, xm = tm, fixed = fixed)

#theta <- optim(theta, loglike, xi = ti, fixed = fixed, control = list(maxit = 5000, trace = 3))$par
#theta <- optim(theta, loglike, xi = ti, fixed = fixed, control = list(maxit = 5000, trace = 3))$par

theta <- optim(theta, loglike.instar, xi = ti, xp = tp, xm = tm, fixed = fixed, control = list(maxit = 1000, trace = 3))$par
theta <- optim(theta, loglike.instar, xi = ti, xp = tp, xm = tm, fixed = fixed, control = list(maxit = 5000, trace = 3))$par
theta <- optim(theta, loglike.instar, xi = ti, xp = tp, xm = tm, fixed = fixed, control = list(maxit = 5000, trace = 3))$par
theta <- optim(theta, loglike.instar, xi = ti, xp = tp, xm = tm, fixed = fixed, control = list(maxit = 5000, trace = 3))$par
theta <- optim(theta, loglike.instar, xi = ti, xp = tp, xm = tm, fixed = fixed, control = list(maxit = 5000, trace = 3))$par

theta <- c(theta, fixed)
#fixed <- c("a.pubescent", "b.pubescent")
#fixed <- theta[fixed]
#theta <- theta[setdiff(names(theta), names(fixed))]

theta <- optim(theta, loglike.instar, xi = ti, xp = tp, xm = tm, control = list(maxit = 1000, trace = 3))$par
for (i in 1:40) theta <- optim(theta, loglike.instar, xi = ti, xp = tp, xm = tm, control = list(maxit = 5000, trace = 3))$par

r <- instar.stats(theta, species = 2526, sex = 2)
rm <- instar.stats(theta, species = 2526, sex = 2, log = FALSE)

clg()
png(file = "instar/Instar mixture fit 1998-2022.png", res = 500, units = "in", height = 11, width = 8.5)
cols <- c("palegreen3", "skyblue4")
xlim <- c(log(10), log(85))
ylim <- c(0, 2000, 400)
t <- seq(0, 5, len = 1000)
m <- kronecker(matrix(c(1,2,3)), matrix(1, ncol = 10, nrow = 10))
m <- rbind(0, cbind(0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
res <- 2
gbarplot(ti, xlim = xlim, ylim = ylim[1:2], xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", col = "grey90", border = "grey80")
vline(log(seq(0, 100, by = 5)), lty = "dashed", col = "grey65")
hline(seq(ylim[1], ylim[2]-ylim[3], by = ylim[3]), lty = "dashed", col = "grey65")
#vline(mu.immature, lty = "dashed", lwd = 2, col = "blue")
v <- rep(0, length(t))
for (i in 1:length(r$immature$mu)){
   d <- 0.01 * sum(ti) * r$immature$p[i] * dnorm(t, r$immature$mu[i],  r$immature$sigma[i])
   v <- v + d
   lines(t, d, lwd = 2, col = cols[1])
   if (r$immature$p[i] > 0.01) text(t[which.max(d)], max(d), names(r$immature$mu)[i], pos = 3, font = 2, cex = 1.0)
} 
lines(t, v, lwd = 2, col = cols[2])
mtext("Immature", 4, 1.0, font = 2)
axis(2, seq(ylim[1], ylim[2], by = ylim[3]))
box(col = "grey50")

gbarplot(tp, xlim = xlim, ylim = ylim[1:2], xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", col = "grey90", border = "grey80")
vline(log(seq(0, 100, by = 5)), lty = "dashed", col = "grey65")
hline(seq(ylim[1], ylim[2]-ylim[3], by = ylim[3]), lty = "dashed", col = "grey65")
#vline(mu.pubescent, lty = "dashed", lwd = 2, col = "blue")
v <- rep(0, length(t))
for (i in 1:length(mu.pubescent)){
   d <- 0.01 * sum(tp) * r$pubescent$p[i] * dnorm(t, r$pubescent$mu[i], r$pubescent$sigma[i])  
   v <- v + d
   lines(t, d, lwd = 2, col = cols[1])
   if (r$pubescent$p[i] > 0.01) text(t[which.max(d)], max(d), names(r$pubescent$mu)[i], pos = 3, font = 2, cex = 1.0)
} 
lines(t, v, lwd = 2, col = cols[2])
mtext("Pubescent", 4, 1.0, font = 2)
mtext(expression(Density~(n~per~km^{2})), 2, 2.5, cex = 1.25)
axis(2, seq(ylim[1], ylim[2]-ylim[3], by = ylim[3]))
box(col = "grey50")

gbarplot(tm, xlim = xlim, ylim = ylim[1:2], xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", col = "grey90", border = "grey80")
vline(log(seq(0, 100, by = 5)), lty = "dashed", col = "grey65")
hline(seq(ylim[1], ylim[2]-ylim[3], by = ylim[3]), lty = "dashed", col = "grey65")
#vline(mu.mature, lty = "dashed", lwd = 2, col = "blue")
v <- rep(0, length(t))
for (i in 1:length(mu.mature)){
   d <- 0.01 * sum(tm) * r$mature$p[i] * dnorm(t, r$mature$mu[i], r$mature$sigma[i])
   v <- v + d
   lines(t, d, lwd = 2, col = cols[1])
   if (r$mature$p[i] > 0.01) text(t[which.max(d)], max(d), names(r$mature$mu)[i], pos = 3, font = 2, cex = 1.0)
} 
lines(t, v, lwd = 2, col = cols[2])
axis(1, at = log(seq(0, 100, by = 5)), labels = seq(0, 100, by = 5), lty = "dashed", col = "grey65")
mtext("Carapace width (mm)", 1, 3.0, cex = 1.25)
mtext("Primiparous", 4, 1.0, font = 2)
axis(2, seq(ylim[1], ylim[2]-ylim[3], by = ylim[3]))
box(col = "grey50")
dev.off()

ix <- which(!is.na(x$carapace.width) & (x$sex == 2) & (x$carapace.width <= 70) & !is.mature(x) & !is.pubescent.scsbio(x))
ip <- which(!is.na(x$carapace.width) & (x$carapace.width >= 28) & (x$carapace.width <= 80) & (x$sex == 2) & !is.mature(x) & is.pubescent.scsbio(x))
im <- which(!is.na(x$carapace.width) & (x$carapace.width >= 36) & (x$carapace.width <= 95) & (x$sex == 2) & is.primiparous.scsbio(x))  

clg()
plot(range(ix), c(0, 80), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", xaxt = "n")
points(ix, xi, col = "brown2", cex = 0.2,)
points(ip, xp, col = "palegreen3", cex = 0.2,)
points(im, xm, col = "skyblue2", cex = 0.2,)
hline(rm$immature$mu, col = "grey20", lty = "dashed", lwd = 2)
box(col = "grey50")

clg()
plot(range(ix), log(c(10, 80)), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i", xaxt = "n")
points(ix, log(xi), col = "brown2", cex = 0.2,)
points(ip, log(xp), col = "palegreen3", cex = 0.2,)
points(im, log(xm), col = "skyblue2", cex = 0.2,)
hline(r$immature$mu, col = "grey20", lty = "dashed", lwd = 2)
hline(r$pubescent$mu, col = "palegreen4", lty = "dashed", lwd = 2)
hline(r$mature$mu, col = "skyblue4", lty = "dashed", lwd = 2)
box(col = "grey50")
