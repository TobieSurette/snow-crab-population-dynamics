library(gulf.data)
library(gulf.graphics)

# Read data:
mu.imm <- read.csv("R/female population model/results/tables/mu.imm.csv")
mu.mat <- read.csv("R/female population model/results/tables/mu.mat.csv")

years <- 1997:(1996+nrow(mu.mat))

clg()

png(file = "R/female population model/results/figures/population abundance.png",
    res = 600, units = "in", width = 6, height = 8)

m <- rbind(0, cbind(0, kronecker(1:2, matrix(1, nrow = 7, ncol = 7)), 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))

plot(c(min(years)-1, max(years)+1), c(0, 2500), 
     type = "n", xlab = "", ylab = "", xaxt = "n", 
     xaxs = "i", yaxs = "i")
grid()

w <- 0.6
for (i in 1:nrow(mu.imm)){
   xx <- c(years[i]-w/2, years[i]+w/2, years[i]+w/2, years[i]-w/2)
   yy <- c(mu.imm$val25.0pc[i], mu.imm$val25.0pc[i], mu.imm$val75.0pc[i], mu.imm$val75.0pc[i])
   if (years[i] > 2025) col = fade("red", 0.4) else col <- "grey80" 
   polygon(xx, yy, col = col, border = "grey30", lwd = 0.5)
   lines(c(years[i], years[i]), c(mu.imm$val2.5pc[i], mu.imm$val25.0pc[i]), col = "grey20", lwd = 0.75)
   lines(c(years[i], years[i]), c(mu.imm$val75.0pc[i], mu.imm$val97.5pc[i]), col = "grey20", lwd = 0.75)
   lines(c(years[i]-w/4, years[i]+w/4), c(mu.imm$val97.5pc[i], mu.imm$val97.5pc[i]), col = "grey20", lwd = 0.75)
   lines(c(years[i]-w/4, years[i]+w/4), c(mu.imm$val2.5pc[i], mu.imm$val2.5pc[i]), col = "grey20", lwd = 0.75)
   lines(c(years[i]-w/2, years[i]+w/2), c(mu.imm$median[i], mu.imm$median[i]), col = "grey20", lwd = 1.5)
}

hline(mean(mu.imm$mean), lty = "dashed", col = fade("red"), lwd = 2)
axis(4, at = mean(mu.imm$mean), labels = "", col = "red")
mtext(round(mean(mu.imm$mean)), 4, 1.0, at = mean(mu.imm$mean), font = 2, cex = 0.85, col = "red")

mtext("Abundance (millions)", 2, 2.5, at = 0, font = 2, cex = 1.25)
mtext("Immature & Pubescents", 3, -2, font = 2, cex = 1.00)

box(col = "grey50")

plot(c(min(years)-1, max(years)+1), c(0, 2500), 
     type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
     xaxs = "i", yaxs = "i")
grid()

w <- 0.6
for (i in 1:nrow(mu.mat)){
   xx <- c(years[i]-w/2, years[i]+w/2, years[i]+w/2, years[i]-w/2)
   yy <- c(mu.mat$val25.0pc[i], mu.mat$val25.0pc[i], mu.mat$val75.0pc[i], mu.mat$val75.0pc[i])
   if (years[i] > 2025) col = fade("red", 0.4) else col <- "grey80" 
   polygon(xx, yy, col = col, border = "grey30", lwd = 0.5)
   lines(c(years[i], years[i]), c(mu.mat$val2.5pc[i], mu.mat$val25.0pc[i]), col = "grey20", lwd = 0.75)
   lines(c(years[i], years[i]), c(mu.mat$val75.0pc[i], mu.mat$val97.5pc[i]), col = "grey20", lwd = 0.75)
   lines(c(years[i]-w/4, years[i]+w/4), c(mu.mat$val97.5pc[i], mu.mat$val97.5pc[i]), col = "grey20", lwd = 0.75)
   lines(c(years[i]-w/4, years[i]+w/4), c(mu.mat$val2.5pc[i], mu.mat$val2.5pc[i]), col = "grey20", lwd = 0.75)
   lines(c(years[i]-w/2, years[i]+w/2), c(mu.mat$median[i], mu.mat$median[i]), col = "grey20", lwd = 1.5)
}

hline(mean(mu.mat$mean), lty = "dashed", col = fade("red"), lwd = 2)

axis(4, at = mean(mu.mat$mean), labels = "", col = "red")
mtext(round(mean(mu.mat$mean)), 4, 1.0, at = mean(mu.mat$mean), font = 2, cex = 0.85, col = "red")
mtext("New & Old-shelled matures", 3, -2, font = 2, cex = 1.00)

axis(1, cex.axis = 1.25)
axis(2, at = 500 * (0:4))

mtext("Year", 1, 3.0, font = 2, cex = 1.25)

box(col = "grey50")

dev.off()

#mu.R <- read.csv("R/female population model/results/tables/mu.R.csv")
#mu.mat <- read.csv("R/female population model/results/tables/mu.mat.csv")
#N.kriging <- read.csv("R/female population model/results/tables/N.kriging.csv")

#years <- min(years):(min(years) + nrow(mu.R) - 1)
#plot(years, mu.R$mean, ylim = c(0, 400), type = "l", col = "green3")
#lines(1997:2025, N.kriging[, 1], col = "blue")
#lines(years+6, 0.35*mu.mat$mean, lwd = 2)


