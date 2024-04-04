library(gulf.data)
library(gulf.graphics)
library(gulf.spatial)
library(mgcv)

years <- 1997:2023

# Load female set data:
maturity <- c("immature", "pubescent", "primiparous", "multiparous")
res <- NULL
for (i in 1:length(years)){
   print(years[i])
   # Load data:
   s <- read.scsset(year = years[i], valid = 1, survey = "regular")
   b <- read.scsbio(year = years[i], sex = 2, survey = "regular")
   b <- b[which((b$carapace.width >= 1) & (b$carapace.width <= 95)), ]
   
   b$maturity <- maturity(b)
   ix <- which((b$maturity == "mature") & is.primiparous(b))
   b$maturity[ix] <- "primiparous"
   ix <- which((b$maturity == "mature") & is.multiparous(b))
   b$maturity[ix] <- "multiparous"

   # Length-frequency from biolgical data:
   f <- freq(b, by = c("date", "tow.id", "maturity"))
   
   for (j in 1:length(maturity)){
      ss <- s
      import(ss, fill = 0) <- f[f$maturity == maturity[j], ]
      fvars <- names(ss)[gsub("[-0-9.]", "", names(ss)) == ""]
      ss[, fvars] <- 1000000 * ss[, fvars] / repvec(ss$swept.area, ncol = length(fvars))
      
      # Disaggregate data using 10x10 minute grid:
      if (years[i] <= 2010){
         ss$grid <- deg2grid(lon(ss), lat(ss))
         ss <- aggregate(ss[,fvars], by = ss[c("date", "grid")], mean)
      }
      
      # Calculate mean densities:
      res <- rbind(res, data.frame(year = years[i], 
                                   maturity = maturity[j], 
                                   carapace.width = as.numeric(fvars), 
                                   density = apply(ss[fvars], 2, mean)))
   }
}

# Square off results:
tmp <- res
res <- array(0, dim = c(year = length(years), carapace.width = 95, maturity = length(maturity)))
dimnames(res) <- list(year = years, carapace.width = 1:95, maturity = maturity)
for (i in 1:length(years)){
   for (j in 1:length(maturity)){
      ix <- which((tmp$year == years[i]) & (tmp$maturity == maturity[j]))
      res[as.character(years[i]),as.character(tmp[ix,"carapace.width"]), maturity[j]] <- tmp$density[ix]
   }
}

png(file = "U de M/figures/Average female size-frequency by maturity.png", units = "in", 
    res = 400, width = 8.5, height = 7)

plot(c(0, 90), c(0, 400), xaxs = "i", yaxs = "i", type = "n", xlab = "", ylab = "")
grid()
cols <- c("green2", "orange", "red2", "red4")
for (i in 1:length(maturity)){
   xx <- as.numeric(dimnames(res)$carapace.width)
   yy <- apply(res[,,i], 2, mean)
   polygon(c(xx, rev(xx)), c(rep(0, length(yy)), rev(yy)), border = cols[i], col = fade(cols[i], 0.4))
}
mtext("Carapace width (mm)", 1, 2.25, cex = 1.25, font = 2, col = "grey25")
mtext("Average density (# / km2)", 2, 2.25, cex = 1.25, font = 2, col = "grey25")

mu <- list(immature = c(10.0, 14.5, 20.5, 28.38, 37.5),
           adolescent = c(41.71,	49.1),
           mature = c(49.1, 58.1, 69.1))

vline(mu$immature, lty = "dashed", col = cols[1], lwd = 2)
vline(mu$adolescent, lty = "dashed", col = cols[2], lwd = 2)
vline(mu$mature, lty = "dashed", col = cols[3], lwd = 2)

mtext(round(mu$immature,1), 1, 0.1, at = mu$immature, cex = 0.75, font = 2)
mtext(as.roman(4:8), 3, 0.1, at = mu$immature, cex = 0.75, font = 2)
mtext(round(mu$adolescent,1), 1, 0.1, at = mu$adolescent, cex = 0.75, font = 2)
mtext(as.roman(8:9), 3, 0.1, at = mu$adolescent, cex = 0.75, font = 2)
mtext(round(mu$mature,1), 1, 0.1, at = mu$mature, cex = 0.75, font = 2)
mtext(as.roman(9:11), 3, 0.1, at = mu$mature, cex = 0.75, font = 2)

legend("topleft", 
       legend = c("Immature", "Adolescent", "New-mature", "Old-mature"),
       pch = 22,
       col = cols, 
       pt.bg = fade(cols), 
       pt.cex = 3,
       cex = 1.25)

box(col = "grey50")

dev.off()