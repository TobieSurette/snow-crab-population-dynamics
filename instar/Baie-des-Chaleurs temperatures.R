
clg()


png(file = paste0("instar/Baie-des-Chaleurs temperatures 2022-2023.png"), res = 500, units = "in", height = 8.5, width = 8.5)

par(mar = c(8,4,4,2))
a <- read.csv("instar/Baie-des-Chaleur 2022.csv")
a <- a[a$Temperature < 5, ]

ix <- seq(1, nrow(a), by = 50)

plot(a$Temperature[ix], ylim = c(-1.5, 4), type = "n", xlab = "", xaxt = "n", ylab = "", yaxs = "i")

a$month <- as.numeric(substr(a$Date, 4,5))

iy <- c(par("usr")[1], which(diff(a$month[ix]) != 0), par("usr")[2])

for (i in 1:(length(iy)-1)){
   rect(iy[i], par("usr")[3], iy[i+1], par("usr")[4], col = c("grey75", "grey95")[(i %% 2) + 1])
}

lines(1:length(ix), a$Temperature[ix])

vline(iy, col = "grey60")
hline(-1:5, col = "grey60", lty = "dashed")
mtext("Temperature (C)", 2, 2.5, cex = 1.25)
axis(1, labels = month.name[c(5:12, 1:5)], at = iy[-length(iy)] + diff(iy)/2, las = 2)


box(col = "grey60")
dev.off()

