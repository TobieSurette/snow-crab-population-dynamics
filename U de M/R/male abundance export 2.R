library(gulf.data)
library(gulf.graphics)
library(gulf.spatial)
library(mgcv)

years <- 1997:2023

# Adolescent proportion function (Surette & Allard, 2011)
adolescent <- function(x, xp = 40.1, s = 0.059)  return(1 / (1 + exp(-4*0.059*(x - 40.1))))

# Load female set data:
maturity <- c("immature", "new-shelled", "old-shelled")
res <- NULL
for (i in 1:length(years)){
  print(years[i])
  # Load data:
  s <- read.scsset(year = years[i], valid = 1, survey = "regular")
  b <- read.scsbio(year = years[i], sex = 1, survey = "regular")
  b <- b[which((b$carapace.width >= 1) & (b$carapace.width <= 145)), ]
  
  b$maturity <- morphometric.maturity(b, probability = TRUE) >= 0.5
  b$maturity <- ifelse(b$maturity, "mature", "immature")
  ix <- which((b$maturity == "mature") & is.new.shell(b))
  b$maturity[ix] <- "new-shelled"
  ix <- which((b$maturity == "mature") & !is.new.shell(b))
  b$maturity[ix] <- "old-shelled"

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

# Add adolescent partition:
ix <- which(res$maturity == "immature")
p <- adolescent(res$carapace.width[ix])
r <- res[ix, ]
res$density[ix] <- (1-p) * res$density[ix]
r$maturity <- "adolescent"
r$density <- p * r$density
res <- rbind(res, r)
maturity <- c(maturity[1], "adolescent", maturity[2:3]) 

# Square off results:
tmp <- res
res <- array(0, dim = c(year = length(years), carapace.width = 150, maturity = length(maturity)))
dimnames(res) <- list(year = years, carapace.width = 1:150, maturity = maturity)
for (i in 1:length(years)){
  for (j in 1:length(maturity)){
    ix <- which((tmp$year == years[i]) & (tmp$maturity == maturity[j]))
    res[as.character(years[i]),as.character(tmp$carapace.width[ix]), maturity[j]] <- tmp$density[ix]
  }
}

png(file = "U de M/figures/Average male size-frequency by maturity.png", units = "in", 
    res = 400, width = 8.5, height = 7)

plot(c(0, 135), c(0, 300), xaxs = "i", yaxs = "i", type = "n", xlab = "", ylab = "")
grid()
cols <- c("green2", "orange", "red2", "red4")
for (i in 1:length(maturity)){
  xx <- as.numeric(dimnames(res)$carapace.width)
  yy <- apply(res[,,i], 2, mean)
  polygon(c(xx, rev(xx)), c(rep(0, length(yy)), rev(yy)), border = cols[i], col = fade(cols[i], 0.4))
}
mtext("Carapace width (mm)", 1, 2.25, cex = 1.25, font = 2, col = "grey25")
mtext("Average density (# / km2)", 2, 2.25, cex = 1.25, font = 2, col = "grey25")

mu <- list(immature = c(10.0, 14.5, 20.5, 28.38, 37.5, 48),
           adolescent = c(40.0,	50.2, 65.0, 78.6, 94.4, 112),
           mature = c(73.4, 99.5, 116.2))

vline(mu$immature, lty = "dashed", col = cols[1], lwd = 2)
vline(mu$adolescent, lty = "dashed", col = cols[2], lwd = 2)

mtext(round(mu$immature,1), 1, 0.1, at = mu$immature, cex = 0.75, font = 2)
mtext(as.roman(4:9), 3, 0.1, at = mu$immature, cex = 0.75, font = 2)
mtext(round(mu$adolescent,1), 1, 0.1, at = mu$adolescent, cex = 0.75, font = 2)
mtext(as.roman(8:13), 3, 0.1, at = mu$adolescent, cex = 0.75, font = 2)


legend("topleft", 
       legend = c("Immature", "Adolescent", "New-mature", "Old-mature"),
       pch = 22,
       col = cols, 
       pt.bg = fade(cols), 
       pt.cex = 3,
       cex = 1.25)

box(col = "grey50")

dev.off()