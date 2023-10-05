library(gulf.data)
library(gulf.graphics)
library(gulf.spatial)

instar <- 11
maturity <- "Mature"
if (tolower(maturity) == "immature")  bg = "brown3"
if (tolower(maturity) == "pubescent") bg = "palegreen3"
if (tolower(maturity) == "mature")    bg = "skyblue3"
  
m <- rbind(0, 0, cbind(0, kronecker(matrix(1:32, ncol = 4), matrix(1, ncol = 6, nrow = 6)), 0), 0)
clg()

png(file = paste0("instar/Instar ", as.roman(instar), " ", maturity, " maps.png"), res = 500, units = "in", height = 11, width = 8.5)
layout(m)
par(mar = c(0,0,0,0))

for (i in 1:32){
  print(i)
  ix <- which((b$year == years[i+1]) & (b$maturity == tolower(maturity)) & (b$instar == instar))
  ss <- s[year(s) == years[i+1], ]
  ss$total <- 0
  if (length(ix) > 0) import(ss, fill = 0) <- catch(b[ix, ]) 
  
  # Immature instar map:
  map.new()
  points(lon(ss), lat(ss), pch = 21, bg = bg, cex = 0.2 * sqrt(ss$total), lwd = 0.5, col = "grey30")
  map("coast")
  text(par("usr")[1] + 0.9 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), years[i+1], cex = 1.25 , font = 2)
  box(col = "grey60")
  
  if (i == 17) mtext(paste0(maturity, " instar ", as.roman(instar)), 3, 0.5, at = par("usr")[1], font = 2)
}
dev.off()
