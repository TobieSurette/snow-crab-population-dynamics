library(gulf.data)
library(gulf.graphics)

# Read data:
mu.R <- read.csv("R/female population model/results/tables/mu.R.csv")

years <- 1997:(1996+nrow(mu.R))

clg()

png(file = "R/female population model/recruitment abundance.png",
    res = 600, units = "in", width = 7, height = 5.5)

plot(c(min(years)-1, max(years)+1), c(0, 1000), type = "n", xlab = "", ylab = "", xaxt = "n", xaxs = "i", yaxs = "i")
grid()

w <- 0.6
for (i in 1:nrow(mu.R)){
   xx <- c(years[i]-w/2, years[i]+w/2, years[i]+w/2, years[i]-w/2)
   yy <- c(mu.R$val25.0pc[i], mu.R$val25.0pc[i], mu.R$val75.0pc[i], mu.R$val75.0pc[i])
   if (years[i] > 2025) col = fade("red", 0.4) else col <- "grey80" 
   polygon(xx, yy, col = col, border = "grey30", lwd = 0.5)
   lines(c(years[i], years[i]), c(mu.R$val2.5pc[i], mu.R$val25.0pc[i]), col = "grey20")
   lines(c(years[i], years[i]), c(mu.R$val75.0pc[i], mu.R$val97.5pc[i]), col = "grey20")
   lines(c(years[i]-w/4, years[i]+w/4), c(mu.R$val97.5pc[i], mu.R$val97.5pc[i]), col = "grey20")
   lines(c(years[i]-w/4, years[i]+w/4), c(mu.R$val2.5pc[i], mu.R$val2.5pc[i]), col = "grey20")
   lines(c(years[i]-w/2, years[i]+w/2), c(mu.R$median[i], mu.R$median[i]), col = "grey20", lwd = 2)
}

axis(1)
mtext("Year", 1, 3.0, font = 2, cex = 1.25)
mtext("Abundance (millions)", 2, 2.5, font = 2, cex = 1.25)

box(col = "grey50")

dev.off()

