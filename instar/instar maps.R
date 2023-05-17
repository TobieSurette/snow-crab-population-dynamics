

instar <- 9
maturity <- "pubescent"

m <- rbind(0, 0, cbind(0, kronecker(matrix(1:32, ncol = 4), matrix(1, ncol = 6, nrow = 6)), 0), 0)
clg()

png(file = paste0("instar/Instar ", as.roman(instar), " ", maturity, " maps.png"), res = 500, units = "in", height = 11, width = 8.5)
layout(m)
par(mar = c(0,0,0,0))

for (i in 1:32){
  print(i)
  ix <- which((b$year == years[i+1]) & (b$maturity == maturity) & (b$instar == instar))
  ss <- s[year(s) == (year+i-4), ]
  ss$total <- 0
  if (length(ix) > 0) import(ss, fill = 0) <- catch(b[ix, ]) 
  
  # Immature instar map:
  map.new()
  points(lon(ss), lat(ss), pch = 21, bg = "brown3", cex = 0.2 * sqrt(ss$total))
  map("coast")
  text(par("usr")[1] + 0.9 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), years[i+1], cex = 1.25 , font = 2)
  box(col = "grey60")
}
dev.off()
