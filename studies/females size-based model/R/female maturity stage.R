library(gulf.data)

years <- 1997:2023

# Load data:
b <- read.scsbio(years, valid = 2, survey = "regular", sex = 2, species = 2526)
b$year <- year(b)

# Identify maturity stages:
b$maturity <- maturity(b)
b$stage <- ""
b$stage[which((b$maturity == "immature") & (b$gonad.colour == 1))] <- "immature"
b$stage[which((b$maturity == "immature") & (b$gonad.colour %in% 2:3))] <- "pubescent"
b$stage[which((b$maturity == "mature") & is.new.shell(b))]  <- "primiparous"
b$stage[which((b$maturity == "mature") & !is.new.shell(b))] <- "multiparous"

# Calculate size frequencies:
fi <- freq(b[b$stage == "immature", ],    by = c("date", "tow.id"))
fp <- freq(b[b$stage == "pubescent", ],   by = c("date", "tow.id"))
fm <- freq(b[b$stage == "primiparous", ], by = c("date", "tow.id"))

# Load tow data and attach size-frequencies:
s <- read.scsset(years, survey = "regular", valid = 1)
si <- sp <- sm <- s
import(si, fill = 0) <- fi
import(sp, fill = 0) <- fp
import(sm, fill = 0) <- fm

# Buffer frequency variables:
fvars <- function(x) return(names(x)[gsub("[0-9]", "", names(x)) == ""])
fvars <- unique(c(fvars(si), fvars(sp), fvars(sm)))
si[setdiff(fvars, names(si))] <- 0
sp[setdiff(fvars, names(sp))] <- 0
sm[setdiff(fvars, names(sm))] <- 0

# Standardize by swept area:
si[fvars] <- si[fvars] / repvec(si$swept.area, ncol = length(fvars))
sp[fvars] <- sp[fvars] / repvec(sp$swept.area, ncol = length(fvars))
sm[fvars] <- sm[fvars] / repvec(sm$swept.area, ncol = length(fvars))

# Summary frequency tables:
ti <- aggregate(si[fvars], by = list(year = year(si)), mean)
tp <- aggregate(sp[fvars], by = list(year = year(sp)), mean)
tm <- aggregate(sm[fvars], by = list(year = year(sm)), mean)

clg()
cols <- c("skyblue3", "green3", "brown2")
theta <- c(intercept = 0.276, transition = 38.2, slope = c(0.15, 0.116), window = 1.6, sigma = 0.135)
#G <- growth.matrix(as.numeric(fvars), theta = theta)  # In 'growth.R'
for (i in 1:(length(years)-2)){
   plot(c(0, 90), c(0, 300), type = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "")
   grid()
   lines(as.numeric(fvars), 1000000 * ti[i, fvars], lwd = 2, col = cols[1])
   lines(as.numeric(fvars), 1000000 * tp[i+1, fvars], lwd = 2, col = cols[2])
   lines(as.numeric(fvars), 1000000 * tm[i+2, fvars], lwd = 2, col = cols[3])
   
   # Growth projections:
   #p <- (t(matrix(as.numeric(tp[i+1, fvars]))) %*% as.matrix(G))[1, ]
   #lines(as.numeric(names(p)), 1000000 * p, lwd = 2, lty = "dashed", col = cols[3])
   
   
   mtext(paste0(years[i:(i+2)], collapse = "-"))
   
   legend("topright", legend = years[i:(i+2)], col = cols, lwd = 2)
   
   box(col = "grey50")
}

# Global patterns:
tiff(file = paste0("results/figures/scs female maturity stages.tiff"),
     compression = "lzw", units = "in", res = 300, height = 6.5, width = 8)
mi <- 1000000 * apply(ti[fvars], 2, mean)
mp <- 1000000 * apply(tp[fvars], 2, mean)
mm <- 1000000 * apply(tm[fvars], 2, mean)
plot(as.numeric(fvars), mi, type = "l", lwd = 2, col = cols[1], 
     xlim = c(0, 90), ylim = c(0, 250), xaxs = "i", yaxs = "i", xlab = "", ylab = "")
lines(as.numeric(fvars), mp, lwd = 2, col = cols[2])
lines(as.numeric(fvars), mm, lwd = 2, col = cols[3])
grid()
mtext("Carapace width (mm)", 1, 2.5, cex = 1.25, font = 1)
mtext(expression(paste("Density (# / km"^"2", ")",sep="")), 2, 2.25, cex = 1.25, font = 1)

box(col = "grey50")

mu <- c(10.5, 15, 21, 28.5, 37.5, 49.5, 58)
names(mu) <- c("IV", "V", "VI", "VII", "VIII", "IX", "X")
vline(mu[1:5], lower = 0, upper = mi[as.character(round(mu[1:5]))], lty = "dashed")
vline(mu[6], lower = 0, upper = mp[as.character(round(mu[6]))], lty = "dashed")
vline(mu[7], lower = 0, upper = mm[as.character(round(mu[7]))], lty = "dashed")

text(mu[1:5], mi[as.character(round(mu[1:5]))]-1, names(mu)[1:5], font = 2, pos = 3)
text(mu[6], mp[as.character(round(mu[6]))]-1, names(mu)[6], font = 2, pos = 3)
text(mu[5]+4.5, mp[as.character(round(mu[5]+4.5))]-1, "VIII ?", font = 2, pos = 3, col = cols[2])
text(mu[7], mm[as.character(round(mu[7]))]-1, names(mu)[7], font = 2, pos = 3)
text(mu[6]+1, mm[as.character(round(mu[6]+1))]-1, "IX ?", font = 2, pos = 3, col = cols[3])

dev.off()

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

logistic <- function(x, xp, window) return(1 / (1 + exp(-(x-xp)/xp)))

theta <- c(intercept = 0.276, transition = 38.2, slope = c(0.32, 0.126), window = 1.6, sigma = 0.135) # Snow crab.

# Immature-to-immature growth transition matrix:
Gi <- growth.matrix(as.numeric(fvars), theta = c(intercept = 2.276, slope = 0.24, log.sigma = -2.5))
Gp <- growth.matrix(as.numeric(fvars), theta = c(intercept = 7.7, slope = 0.126, log.sigma = -2))
plot(as.numeric(fvars), mi, type = "n", 
     xlim = c(0, 90), ylim = c(0, 250), xaxs = "i", yaxs = "i", xlab = "", ylab = "")
grid()
lines(as.numeric(fvars), mi, lwd = 2, col = cols[1])
lines(as.numeric(colnames(Gi)), as.numeric(mi) %*% Gi, lwd = 2, col = cols[1], lty = "dashed")
lines(as.numeric(fvars), mp, lwd = 2, col = cols[2])
lines(as.numeric(colnames(Gp)), as.numeric(mi) %*% Gp, lwd = 2, col = cols[2], lty = "dashed")
vline(mu, col = "brown3")

legend("topright", 
       legend = c("Immature observed", "Immature to pubescent projection", "Immature to pubescent projection", "Pubescent observed"),
       lwd = 2, col = c(cols[1], cols[1], cols[2], cols[2]),
       lty = c("solid", "dashed", "dashed", "solid"))

mtext("Carapace width (mm)", 1, 2.5, cex = 1.25, font = 1)
mtext(expression(paste("Density (# / km"^"2", ")",sep="")), 2, 2.25, cex = 1.25, font = 1)

box(col = "grey50")

# Immature-to-pubescent growth transition matrix:
Gp <- growth.matrix(as.numeric(fvars), theta = c(intercept = 7.7, slope = 0.126, log.sigma = -2))
gbarplot(mp, ylim = c(0, 200), xlim = c(0, 90), col = "grey90", border = "grey70")
lines(as.numeric(colnames(Gp)), as.numeric(mi) %*% Gp, lwd = 2, col = cols[1], lty = "dashed")
lines(as.numeric(colnames(Gi)), as.numeric(mi) %*% Gi, lwd = 2, col = "blue")
lines(as.numeric(fvars), mi, lwd = 2, col = cols[1])

lines(as.numeric(fvars), mp, lwd = 2, col = cols[2])
vline(mu, col = "brown3")

# Plot growth models:
plot(c(0, 90), c(0, 40), type = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "")
lines(as.numeric(rownames(Gi)), apply(Gi, 1, function(x) sum(x * as.numeric(colnames(Gi)))) - as.numeric(rownames(Gi))) 
lines(as.numeric(rownames(Gp)), apply(Gp, 1, function(x) sum(x * as.numeric(colnames(Gp)))) - as.numeric(rownames(Gp))) 



