

year <- 2019

m <- rbind(0, 0, cbind(0, kronecker(t(matrix(1:18, ncol = 6)), matrix(1, ncol = 6, nrow = 6)), 0), 0)
clg()

png(file = paste0("instar/Year class map ", year, "-", year+7, ".png"), res = 500, units = "in", height = 11, width = 8.5)

layout(m)
par(mar = c(0,0,0,0))

for (i in 4:11){
   ix <- which((b$year == (year+i-4)) & (b$maturity == "immature") & (b$instar == i))
   ss <- s[year(s) == (year+i-4), ]
   ss$total <- 0
   if (length(ix) > 0) import(ss, fill = 0) <- catch(b[ix, ]) 
   
   # Immature instar map:
   map.new()
   if (i < 10){
     points(lon(ss), lat(ss), pch = 21, bg = "brown3", cex = 0.4 * sqrt(ss$total))
     map("coast")
     text(par("usr")[1] + 0.9 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), (year+i-4), cex = 1.25 , font = 2)
     text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), paste0("Instar ", as.roman(i)), cex = 1.25 , font = 2)
     if (i == 4){
       mtext("Immature", 2, 0.5, col = "brown3", font = 2)
       mtext("Immature", 3, 0.5, col = "brown3", font = 2)
     } 
     box(col = "grey60")
   }

   if (i == 7) {map.new(); map.new()}
   if (i >= 8){
      map.new()
      if (i < 11){
         ix <- which((b$year == (year+i-4)) & (b$maturity == "pubescent") & (b$instar == i))
         ss <- s[year(s) == (year+i-4), ]
         ss$total <- 0
         if (length(ix) > 0) import(ss, fill = 0) <- catch(b[ix, ]) 
        
         # Pubescent instar map:
         points(lon(ss), lat(ss), pch = 21, bg = "palegreen3", cex = 0.4 * sqrt(ss$total))
         map("coast")
         text(par("usr")[1] + 0.9 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), (year+i-4), cex = 1.25 , font = 2)
         text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), paste0("Instar ", as.roman(i)), cex = 1.25 , font = 2)
         if (i == 8) mtext("Pubescent", 3, 0.5, font = 2, col = "palegreen3")
         box(col = "grey60")
      }
     
     if (i >= 9){
       ix <- which((b$year == (year+i-4)) & (b$maturity == "mature") & (b$instar == i))
       ss <- s[year(s) == (year+i-4), ]
       ss$total <- 0
       if (length(ix) > 0) import(ss, fill = 0) <- catch(b[ix, ]) 
       
       map.new()
       points(lon(ss), lat(ss), pch = 21, bg = "skyblue3", cex = 0.4 * sqrt(ss$total))
       map("coast")
       text(par("usr")[1] + 0.9 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), (year+i-4), cex = 1.25 , font = 2)
       text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]), par("usr")[3] + 0.9 * diff(par("usr")[3:4]), paste0("Instar ", as.roman(i)), cex = 1.25 , font = 2)
       if (i == 9) mtext("Primiparous", 3, 0.5, font = 2, col = "skyblue3")
       box(col = "grey60")
     }else{
       map.new()
     }
   }
}
dev.off()

