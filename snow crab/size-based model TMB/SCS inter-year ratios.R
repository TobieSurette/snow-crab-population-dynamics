# Read SCS data:
load("results/tables/size-frequencies males scs 2001-2022.rdata")
fi <- f[,,c("immature")] + f[,,c("skip")]
fm <- f[,,c("recruit")] + f[,,c("residual")]
fi <- fi[,as.character(1:150)]
fm <- fm[,as.character(1:150)]

# Define growth model:
growth.matrix <- function(x, theta, ymax, ...){
   ux <- sort(unique(x))
   xmax <- max(ux)
   dx <- min(diff(ux))
   x0 <- seq(min(ux), max(ux), by = dx)
   if (length(x0) > 1000) stop("Growth matrix dimensions exceed 1000x1000.")
   
   # Define growth means and standard errors:
   m <- growth(theta = theta)(x0)
   s <- growth(theta = theta, error = TRUE)(x0)
   
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

# Calculate growth matrix:
theta <- c(intercept = 0.276, transition = 38.2, slope = c(0.32, 0.100), window = 1.6, sigma = 0.135)
G <- growth.matrix(1:150, theta)

# Calculate growth matrix:
fig <- fi %*% G 
fig <- fig[, as.character(1:150)]

# Calculate residuals:
tab <- NULL
tab2 <- NULL
for (i in 1:(nrow(f)-1)){
   tab <- rbind(tab, log(fi[i+1,] / fig[i,]))
   tab2 <- rbind(tab2, (fi[i+1,] - fig[i,]) / fi[i+1,])
}
rownames(tab) <- rownames(f)[-1] 
colnames(tab) <- 1:150
rownames(tab2) <- rownames(f)[-1] 
colnames(tab2) <- 1:150
tab[!is.finite(tab)] <- NA
tab[tab > 2] <- 2
tab[tab < -2] <- -2
tab2[!is.finite(tab2)] <- NA
tab2[tab2 < -2] <- -2


image(as.numeric(rownames(tab)), 1:150, tab, 
      col = hcl.colors(100, "YlOrRd", rev = TRUE), xlab = "", ylab = "")
hline(95, lty = "dashed")

mu <- apply(tab, 2, mean, na.rm= TRUE)
sigma <- apply(tab, 2, sd, na.rm= TRUE)
t <- (tab - repvec(mu, nrow = nrow(tab))) / repvec(sigma, nrow = nrow(tab))
image(as.numeric(rownames(tab)), 1:150, t, 
      col = c(colorRampPalette(c("blue", "white"))(50), colorRampPalette(c("white", "red"))(50)), 
      ylim = c(0, 130), xlab = "", ylab = "")
hline(95, lty = "dashed", lwd = 2)

# Relative difference between observed and predicted size-frequencies:
tiff(file = paste0("R/population model/figures/SCS year x year ratios relative.tiff"), 
     compression = "lzw", units = "in", res = 300, height = 7.5, width = 7.5)
image(as.numeric(rownames(tab2)), 1:150, tab2, 
      col = c(colorRampPalette(c("blue", "white"))(100), colorRampPalette(c("white", "red"))(50)), 
      breaks = seq(-2, 1, len = 151), 
      xlab = "", ylab = "", xaxt = "n", yaxt = "n", ylim = c(0, 120))
hline(95, lty = "dashed", lwd = 2)
axis(2, at = seq(0, 140, by = 10))
axis(1, at = seq(min(as.numeric(rownames(tab2))), max(as.numeric(rownames(tab2))), by = 4))
axis(1, at = seq(min(as.numeric(rownames(tab2))) + 2, max(as.numeric(rownames(tab2))), by = 4))
mtext("Year", 1, 2.5, cex = 1.25)
mtext("Carapace width (mm)", 2, 2.5, cex = 1.25)
box(col = "grey60")
dev.off()

# Average relative difference between observed and predicted size-frequencies:
tiff(file = paste0("R/population model/figures/SCS year x year ratios average.tiff"), 
     compression = "lzw", units = "in", res = 300, height = 7.5, width = 7.5)
gbarplot(100 * apply(tab2, 2, mean, na.rm = TRUE), grid = TRUE, 
         xaxs = "i", yaxs = "i", xlim = c(0, 120), ylim = 100 * c(-2, 1))
vline(95, lty = "dashed", lwd = 2)
mtext("Relative difference (%)", 2, 2.5, cex = 1.25)
mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
box(col = "grey60")
dev.off()

plot(c(30, 120), c(-1, 1.5), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")
grid()
for (i in 1:(nrow(f)-1)){
   lines(as.numeric(colnames(f)), log(f[i+1,] / fg[i,]), col = rainbow(nrow(f)-1)[i])
}
box(col = "grey60")
hline(0, col = "red", lwd = 2)
vline(95, col = "red", lwd = 2, lty = "dashed")


r <- log(f[-1, ] / fg[-nrow(fg),])
n <- f[-1, ] + fg[-nrow(fg),]
ix <- is.finite(r) & !is.na(r)
data <- data.frame(year = as.factor(repvec(as.numeric(rownames(r)), ncol = ncol(r))[ix]),
                   cw = repvec(as.numeric(colnames(r)), nrow = nrow(r))[ix],
                   ratio = r[ix],
                   n = n[ix])
model <- mgcv::gam(ratio ~ s(cw) + year, weight = n, data = data)
res <- as.data.frame(predict(model, newdata = data.frame(year = sort(unique(data$year)), cw = 90), se.fit = TRUE))
res$year <- as.numeric(as.character(sort(unique(data$year))))

plot(model)
grid()
hline(0, col = "red", lwd = 2)
vline(95, col = "red", lwd = 2, lty = "dashed")


gbarplot(res$fit- mean(res$fit), res$year, xaxt = "n", grid = TRUE)
error.bar(res$year, lower = res$fit- mean(res$fit) - 1.96 * res$se.fit, upper = res$fit- mean(res$fit) + 1.96 * res$se.fit)
axis(1, at = seq(2000, 2022, by = 2))
mtext("Year", 1, 2.5, cex = 1.25)
mtext("log(y+1 / y)", 2, 2.5, cex = 1.25)
box(col = "grey60")

# Read SCS data:
load("results/tables/size-frequencies males scs 2001-2022.rdata")
f <- f[,as.character(1:150),]
f[,,"immature"] <- f[,,"immature"] + f[,,"skip"]
f[,,"recruit"] <- f[,,"recruit"] + f[,,"residual"]
f <- f[,,c("immature", "recruit")]
dimnames(f) <- c(dimnames(f)[1:2], list(group = c("immature", "mature")))

fg <- (f[,,"immature"] %*% G)[-dim(f)[1], as.character(1:150)] + f[-dim(f)[1],,"mature"]
f <- apply(f, 1:2, sum)
f <- f[-1, ]
r <- log(f / fg)
n <- f + fg
ix <- is.finite(r) & !is.na(r)
data <- data.frame(year = as.factor(repvec(as.numeric(rownames(r)), ncol = ncol(r))[ix]),
                   cw = repvec(as.numeric(colnames(r)), nrow = nrow(r))[ix],
                   ratio = r[ix],
                   n = n[ix])
model <- mgcv::gam(ratio ~ s(cw) + year, weight = n, data = data)
res <- as.data.frame(predict(model, newdata = data.frame(year = sort(unique(data$year)), cw = 90), se.fit = TRUE))
res$year <- as.numeric(as.character(sort(unique(data$year))))
plot(model, xlim = c(20, 125), ylim = c(-1.25, 2), ylab = "", xlab = "", xaxs = "i")
grid()
hline(0, col = "red", lwd = 2)
vline(95, col = "red", lwd = 2, lty = "dashed")
mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
mtext("Log(f/f_pred)", 2, 2.5, cex = 1.25)
box(col = "grey60")

file = paste0("results/figures/size-frequencies/SCS inter-year ratios male GAM - english.tiff")
tiff(file = file, compression = "lzw", units = "in", res = 300, height = 7.5, width = 7.5)
gbarplot(res$fit- mean(res$fit), res$year, xaxt = "n", grid = TRUE)
error.bar(res$year, lower = res$fit- mean(res$fit) - 1.96 * res$se.fit, upper = res$fit- mean(res$fit) + 1.96 * res$se.fit)
axis(1, at = seq(2000, 2022, by = 2))
mtext("Year", 1, 2.5, cex = 1.25)
mtext("log(y+1 / y)", 2, 2.5, cex = 1.25)
box(col = "grey60")
dev.off()

plot(c(0, 150), c(-1, 1), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")
grid()
for (i in 1:nrow(f)){
   lines(as.numeric(colnames(f)), log(f[i,] / fg[i,]), col = rainbow(nrow(f)-1)[i])
}
box(col = "grey60")




