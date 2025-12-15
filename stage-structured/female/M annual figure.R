library(gulf.data)
library(gulf.graphics)

# Read data:
M.mean <- read.csv("R/female population model/results/tables/M.mean.csv")

clg()

png(file = "R/female population model/M annual.png",
    res = 600, units = "in", width = 7, height = 5.5)

par(mar = c(c(10, 4, 2, 2) + 0.1))

plot(c(0.5, 10.5), c(0, 0.7), type = "n", xlab = "", ylab = "", xaxt = "n", xaxs = "i", yaxs = "i")
grid()

w <- 0.6
for (i in 1:nrow(M.mean)){
   xx <- c(i-w/2, i+w/2, i+w/2, i-w/2)
   yy <- c(M.mean$val25.0pc[i], M.mean$val25.0pc[i], M.mean$val75.0pc[i], M.mean$val75.0pc[i])
   if (i < 3)  col <- "grey80" 
   if (i == 3) col <- fade("orange")
   if (i > 3)  col <- fade("red", 0.40)
   polygon(xx, yy, col = col, border = "grey30", lwd = 0.5)
   lines(c(i, i), c(M.mean$val2.5pc[i], M.mean$val25.0pc[i]), col = "grey20")
   lines(c(i, i), c(M.mean$val75.0pc[i], M.mean$val97.5pc[i]), col = "grey20")
   lines(c(i-w/4, i+w/4), c(M.mean$val97.5pc[i], M.mean$val97.5pc[i]), col = "grey20")
   lines(c(i-w/4, i+w/4), c(M.mean$val2.5pc[i], M.mean$val2.5pc[i]), col = "grey20")
   lines(c(i-w/2, i+w/2), c(M.mean$mean[i], M.mean$mean[i]), col = "grey20", lwd = 2)
}

axis(1, at = 1:10, 
     labels = c("Imm → Imm", "Imm → Pub", "Pub → New",
                "New → Old 1yr", "Old 1yr → Old 2yr",
                "Old 2yr → Old 3yr", "Old 3yr → Old 4yr",
                "Old 4yr → Old 5yr", "Old 5yr → Old 6yr",
                "Old 6yr → Old 7yr"), las = 2)
mtext("Annual mortality", 2, 2.5, font = 2, cex = 1.25)
mtext("Stage transition", 1, 8.5, font = 2, cex = 1.25)

box(col = "grey50")

dev.off()

