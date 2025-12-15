library(gulf.data)
library(gulf.graphics)

# Read data:
S <- read.csv("R/female population model/results/tables/S.csv")

clg()

png(file = "R/female population model/selectivity.png",
    res = 600, units = "in", width = 7, height = 5.5)
par(mar = c(c(10, 4, 2, 2) + 0.1))

plot(c(0.5, 6.5), c(0, 1), type = "n", xlab = "", ylab = "", xaxt = "n", xaxs = "i", yaxs = "i")
grid()

w <- 0.6
for (i in 1:nrow(S)){
   xx <- c(i-w/2, i+w/2, i+w/2, i-w/2)
   yy <- c(S$val25.0pc[i], S$val25.0pc[i], S$val75.0pc[i], S$val75.0pc[i])
   if (i == 1)  col <- "grey80" 
   if (i == 2) col <- fade("orange")
   if (i > 2)  col <- fade("red", 0.40)
   polygon(xx, yy, col = col, border = "grey30", lwd = 0.5)
   lines(c(i, i), c(S$val2.5pc[i], S$val25.0pc[i]), col = "grey20")
   lines(c(i, i), c(S$val75.0pc[i], S$val97.5pc[i]), col = "grey20")
   lines(c(i-w/4, i+w/4), c(S$val97.5pc[i], S$val97.5pc[i]), col = "grey20")
   lines(c(i-w/4, i+w/4), c(S$val2.5pc[i], S$val2.5pc[i]), col = "grey20")
   lines(c(i-w/2, i+w/2), c(S$mean[i], S$mean[i]), col = "grey20", lwd = 2)
}

axis(1, at = 1:6, 
     labels = c("Immature", "Pubescent", 
                "Emy Serge D.\n(1997-1998)", "Den C. Martin\n(1999-2002)", 
                "Marco-Michel\n(2003-2012)", "Jean-Mathieu\n(2013-2018)"), las = 2)
mtext("Catchability", 2, 2.5, font = 2, cex = 1.25)
mtext("Stage / Survey vessel", 1, 7.0, font = 2, cex = 1.25)

box(col = "grey50")

dev.off()

