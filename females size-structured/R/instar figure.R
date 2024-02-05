library(gulf.data)
library(gulf.graphics)
library(TMB)

# Prepare data:
years <- 1998:2022
s <- read.scsset(years, valid = 1, survey = "regular")
b <- read.scsbio(years)
b <- b[!is.na(match(b[key(s)], s[key(s)])), ]
s <- s[order(s$date, s$tow.number), ]
s$group <- match(s[key(s)], unique(s[key(s)]))-1
b$group <- s$group[match(b[key(s)], s[key(s)])]
b$year <- year(b)
b <- b[order(b$date, b$tow.number), ]
      
# Define maturity stages:
b$maturity <- ""
b$maturity[which(!is.na(b$carapace.width) & (b$sex == 2) & (b$carapace.width <= 70) & !is.mature(b) & !is.pubescent.scsbio(b))] <- "immature"
b$maturity[which(!is.na(b$carapace.width) & (b$carapace.width >= 28) & (b$carapace.width <= 80) & (b$sex == 2) & !is.mature(b) & is.pubescent.scsbio(b))] <- "pubescent"
b$maturity[which(!is.na(b$carapace.width) & (b$carapace.width >= 36) & (b$carapace.width <= 95) & (b$sex == 2) & is.primiparous.scsbio(b))] <- "mature"
b <- b[b$maturity != "", ]
b$ix <- 1:nrow(b)

clg()
png(file = paste0("instar/Instar plot ", min(years), "-", max(years), ".png"), res = 300, units = "in", height = 8.5, width = 11)
m <- cbind(matrix(1, nrow = 8, ncol = 8), matrix(2, nrow = 8, ncol = 1))
m <- rbind(0, cbind(0, m, 0), 0)
layout(m)
par(mar = c(0,0,0,0))

plot(c(0, max(b$ix)), log(c(10, 85)), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i")
z <- c(0, b$ix[which(diff(b$year) > 0)], par("usr")[2])
for (i in 1:(length(z)-1)){
   polygon(c(z[i], z[i+1], z[i+1], z[i]), c(par("usr")[c(3,3)], par("usr")[c(4,4)]), col = c("grey95", "grey85")[(i %% 2) + 1], border = "grey50")
}
cols <- c("brown2", "palegreen3", "blue")
mats <- c("immature", "pubescent", "mature")
points(b$ix, log(b$carapace.width), pch = 21, cex = 0.08, bg = cols[match(b$maturity, mats)], col = cols[match(b$maturity, mats)])
iy <- seq(1, length(years), by = 2)
axis(1, years[iy], at = (z[-length(z)] + diff(z) / 2)[iy], las = 2)
iy <- seq(2, length(years), by = 2)
axis(1, years[iy], at = (z[-length(z)] + diff(z) / 2)[iy], las = 2)

at <- seq(10, 80, by = 10)
axis(2, at = log(at), labels = at)
at <- seq(15, 85, by = 10)
axis(2, at = log(at), labels = at)

mtext("Carapace width (mm)", 2, 2.5, cex = 1.25)
mtext("Year", 1, 4.0, cex = 1.25)

box(col = "grey60")

ix <- which(b$maturity == "immature")
tab <- aggregate(b$carapace.width[ix], by = list(cw = round(b$carapace.width[ix]*5)/5), length)
plot(c(0, 1.1 * max(tab$x)), par("usr")[3:4], type = "n", xaxs = "i", yaxs = "i", xaxt = "n", yaxt = "n", bty = "n")
polygon(c(tab$x) + 50, log(c(tab$cw)), col = adjustcolor(cols[1], alpha.f = 0.5) , border = "grey70", lwd = 0.5)
text(max(tab$x), log(tab$cw[which.max(tab$x)]), "immature", pos = 3, srt = -90, font = 2, cex = 1.25)
  
ix <- which(b$maturity == "pubescent")
tab <- aggregate(b$carapace.width[ix], by = list(cw = round(b$carapace.width[ix]*5)/5), length)
polygon(c(tab$x) + 50, log(c(tab$cw)), col = adjustcolor(cols[2], alpha.f = 0.5), border = "grey70", lwd = 0.5)
text(max(tab$x), log(tab$cw[which.max(tab$x)]), "pubescent", pos = 3, srt = -90, font = 2, cex = 1.25)

ix <- which(b$maturity == "mature")
tab <- aggregate(b$carapace.width[ix], by = list(cw = round(b$carapace.width[ix]*5)/5), length)
polygon(c(tab$x) + 50, log(c(tab$cw)), col = adjustcolor(cols[3], alpha.f = 0.5), border = "grey70", lwd = 0.5)
text(max(tab$x), log(tab$cw[which.max(tab$x)]), "mature", pos = 3, srt = -90, font = 2, cex = 1.25)

dev.off()
