library(gulf.data)
library(gulf.graphics)

# Read data:
mu.kriging <- read.csv("R/female population model/results/tables/N.kriging.csv")
mu.N <- read.csv("R/female population model/results/tables/mu.N.csv")

years <- 1997:2029

png(file = "R/female population model/survey abundance.png",
    res = 500, units = "in", width = 7, height = 7)

plot(c(min(years), 2025), c(0, 800), type = "n", xlab = "", ylab = "", yaxs = "i")
grid()
lines(years[years <= 2025], mu.N$mean[grep(",1]", mu.N$variable)], lwd = 2, col = "orange")
lines(years[years <= 2025], mu.kriging[, 1], lty = "dashed", lwd = 2, col = "orange")

lines(years[years <= 2025], mu.N$mean[grep(",3]", mu.N$variable)], lwd = 2, col = fade("red"))
lines(years[years <= 2025], mu.kriging[, 3], lty = "dashed", lwd = 2, col = fade("red"))

lines(years[years <= 2025], mu.N$mean[grep(",2]", mu.N$variable)], lwd = 2, col = "green3")
lines(years[years <= 2025], mu.kriging[, 2], lty = "dashed", lwd = 2, col = "green3")

lines(years[years <= 2025], mu.N$mean[grep(",4]", mu.N$variable)], lwd = 2, col = fade("blue"))
lines(years[years <= 2025], mu.kriging[, 4], lty = "dashed", lwd = 2, col = fade("blue"))

legend("topleft", 
       legend = c("Immature", "Pubescent", "New-shelled", "Old-shelled"),
       lwd = 2, col = c("orange", fade("red"), "green3", fade("blue")),
       box.lwd = 0.5, box.col = "grey")

mtext("Year", 1, 2.5, font = 2, cex = 1.25)
mtext("Abundance (millions)", 2, 2.5, font = 2, cex = 1.25)

box(col = "grey50")

dev.off()