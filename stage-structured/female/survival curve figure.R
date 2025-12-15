library(gulf.data)
library(gulf.graphics)

r <- read.csv("R/female population model/results/tables/survival.csv")

r <- rbind(r[1,], r)
r$variable[1] <- "survival[0]"
r[1, -1] <- 0
r[1, c("mean", "median")] <- 1

png(file = "R/female population model/survival curve.png",
    res = 600, units = "in", width = 7, height = 7)
gbarplot(100 * r$mean, (1:8) - 0.5, xaxt = "n", grid = TRUE, col = fade("grey50"), ylim = c(0, 100))
axis(1, at = (1:8) - 0.5)
mtext("Elapsed time since terminal moult (years)", 1, 2.5, cex = 1.25, font = 2)
mtext("Survival fraction (%)", 2, 2.5, cex = 1.25, font = 2)

vline(3.284, lty = "dashed", col = fade("red"), lwd = 2)
polygon(c(3.168, 3.425, 3.425, 3.168), c(0, 0, 100, 100), col = fade("red", 0.40), border = NA)	
polygon(c(2.801, 3.635, 3.635, 2.801), c(0, 0, 100, 100), col = fade("red", 0.20), border = NA)	
axis(3, at = 3.361, label = "Life expectancy = 3.28 yrs", col = "red", padj = 0.5, col.axis = "red")

box(col = "grey50")

dev.off()


