library(gulf.data)
library(gulf.graphics)

# Read data:
M <- read.csv("R/female population model/results/tables/M.csv")
M <- M[grep(",1]", M$variable), ]

years <- 1997:(1996+nrow(M))

clg()

png(file = "R/female population model/mortality imm - imm figure.png",
    res = 600, units = "in", width = 7, height = 5.5)

plot(c(min(years)-1, max(years)+1), c(0, 1), type = "n", xlab = "", ylab = "", xaxt = "n", xaxs = "i", yaxs = "i")
grid()

w <- 0.6
for (i in 1:nrow(M)){
   xx <- c(years[i]-w/2, years[i]+w/2, years[i]+w/2, years[i]-w/2)
   yy <- c(M$val25.0pc[i], M$val25.0pc[i], M$val75.0pc[i], M$val75.0pc[i])
   if (years[i] > 2025) col = fade("red", 0.4) else col <- "grey80" 
   polygon(xx, yy, col = col, border = "grey30", lwd = 0.5)
   lines(c(years[i], years[i]), c(M$val2.5pc[i], M$val25.0pc[i]), col = "grey20")
   lines(c(years[i], years[i]), c(M$val75.0pc[i], M$val97.5pc[i]), col = "grey20")
   lines(c(years[i]-w/4, years[i]+w/4), c(M$val97.5pc[i], M$val97.5pc[i]), col = "grey20")
   lines(c(years[i]-w/4, years[i]+w/4), c(M$val2.5pc[i], M$val2.5pc[i]), col = "grey20")
   lines(c(years[i]-w/2, years[i]+w/2), c(M$median[i], M$median[i]), col = "grey20", lwd = 2)
}

axis(1)
mtext("Year", 1, 3.0, font = 2, cex = 1.25)
mtext("Annual mortality", 2, 2.5, font = 2, cex = 1.25)

box(col = "grey50")

dev.off()

