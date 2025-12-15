library(gulf.data)
library(gulf.graphics)

# Read data:
p.mat <- read.csv("R/female population model/results/tables/p.mat.csv")

years <- 1997:(1996+nrow(p.mat))

clg()

png(file = "R/female population model/pubescent maturation.png",
    res = 600, units = "in", width = 7, height = 5.5)

plot(c(min(years)-1, max(years)+1), c(0, 1), type = "n", xlab = "", ylab = "", xaxt = "n", xaxs = "i", yaxs = "i")
grid()

w <- 0.6
for (i in 1:nrow(p.mat)){
   xx <- c(years[i]-w/2, years[i]+w/2, years[i]+w/2, years[i]-w/2)
   yy <- c(p.mat$val25.0pc[i], p.mat$val25.0pc[i], p.mat$val75.0pc[i], p.mat$val75.0pc[i])
   if (years[i] > 2025) col = fade("red", 0.4) else col <- "grey80" 
   polygon(xx, yy, col = col, border = "grey30", lwd = 0.5)
   lines(c(years[i], years[i]), c(p.mat$val2.5pc[i], p.mat$val25.0pc[i]), col = "grey20")
   lines(c(years[i], years[i]), c(p.mat$val75.0pc[i], p.mat$val97.5pc[i]), col = "grey20")
   lines(c(years[i]-w/4, years[i]+w/4), c(p.mat$val97.5pc[i], p.mat$val97.5pc[i]), col = "grey20")
   lines(c(years[i]-w/4, years[i]+w/4), c(p.mat$val2.5pc[i], p.mat$val2.5pc[i]), col = "grey20")
   lines(c(years[i]-w/2, years[i]+w/2), c(p.mat$median[i], p.mat$median[i]), col = "grey20", lwd = 2)
}

axis(1)
mtext("Year", 1, 3.0, font = 2, cex = 1.25)
mtext("Pubescent maturation rate", 2, 2.5, font = 2, cex = 1.25)

box(col = "grey50")

dev.off()

