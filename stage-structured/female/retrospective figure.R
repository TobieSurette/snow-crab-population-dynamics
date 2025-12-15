library(gulf.data)
library(gulf.graphics)

# Read data:
mu <- read.csv("R/female population model/results/tables/mu.csv")

# Read retrospective fits:
mu.r <- list()
files <- dir(path = "R/female population model/", pattern = "*mu-[0-9].csv", full.names = TRUE)
for (i in 1:length(files)) mu.r[[i]] <- read.csv(files[i])

years <- 1997:2029

clg()

png(file = "R/female population model/retrospective analysis.png",
    res = 600, units = "in", width = 9, height = 9)

m <- kronecker(matrix(2:5, ncol = 2), matrix(1, nrow = 7, ncol = 7))
m <- rbind(1, m)
m <- rbind(0, cbind(0, 0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))

plot(c(0, 1), c(0, 1), type = "n", xaxt = "n", yaxt = "n", axes = FALSE)
cols <- rainbow(length(mu.r))
legend("center",
       legend = 2025-(1:6),
       cex = 1.25,
       col = cols,
       lwd = 2, box.lwd = 0.5, box.col = "grey50",
       horiz = TRUE)

for (j in 1:4){
   plot(c(min(years), 2030), c(0, 1400), type = "n", 
        xlab = "", ylab = "", yaxs = "i", xaxt = "n", yaxt = "n")
   grid()
   r <- mu[grep(paste0(",", j, "]"), mu$variable), ]
   lines(years[years <= 2025], r$mean[years <= 2025], lwd = 2, col = "grey20")
   lines(years[years >= 2025], r$mean[years >= 2025], lty = "dashed", lwd = 2, col = "grey20")
   ix <- which(years <= 2025)
   polygon(c(years[ix], rev(years[ix])), c(r$val97.5pc[ix], rev(r$val2.5pc[ix])), col = fade("grey20", 0.15), border = NA)
   ix <- which(years >= 2025)
   polygon(c(years[ix], rev(years[ix])), c(r$val97.5pc[ix], rev(r$val2.5pc[ix])), col = fade("grey20", 0.30), border = NA)
   

   for (i in 1:length(mu.r)){
      r <- mu.r[[i]][grep(paste0(",", j, "]"), mu.r[[i]]$variable), ]
      ix <- which(years >= (2020-i) & years <= (2025-i))
      lines(years[ix], r$mean[ix], col = fade(cols[i], 0.85), lwd = 2)
      ix <- which(years >= (2025-i) & years <= (2029-i))
      lines(years[ix], r$mean[ix], col = fade(cols[i], 0.85), lwd = 2, lty = "dashed")
   }
   
   #legend("topleft", 
   #       legend = c("Immature", "Pubescent", "New-shelled", "Old-shelled"),
   #       lwd = 2, col = c("orange", fade("red"), "green3", fade("blue")),
   #       box.lwd = 0.5, box.col = "grey")
   
   if (j %in% c(4)) mtext("Year", 1, 3.0, at = par("usr")[1], font = 2, cex = 1.25)
   if (j %in% c(1)) mtext("Abundance (millions)", 2, 2.5, at = 0, font = 2, cex = 1.25)
   if (j %in% c(2, 4)) axis(1)
   if (j == 1) axis(2)
   if (j == 2) axis(2, at = 200 * (0:6))
   
   if (j == 1) mtext("Immature", 3, -2, cex = 1.0, font = 2)
   if (j == 2) mtext("Pubescent", 3, -2, cex = 1.0, font = 2)
   if (j == 3) mtext("New-shelled mature", 3, -2, cex = 1.0, font = 2)
   if (j == 4) mtext("Old-shelled mature", 3, -2, cex = 1.0, font = 2)
   
   box(col = "grey50")
}

dev.off()