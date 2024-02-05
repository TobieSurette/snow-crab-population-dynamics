library(gulf.data)
library(gulf.graphics)

years <- 2012:2014

clg()
png(file = paste0("Female size-frequencies ", years[1], "-", max(years), ".png"), res = 500, units = "in", height = 11.0, width = 8.5)
m <- kronecker(matrix(1:3), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0, cbind(0, m, 0), 0)
par(mar = c(0,0,0,0))
layout(m)
for (i in 1:length(years)){
  year <- years[i]
  
  # Load immature females (y-2)
  s <- read.scsset((year-2):(year), valid = 1, survey = "regular")
  b <- read.scsbio((year-2):year, sex = 2)
  
  ix <- which(!is.mature(b) & (b$gonad.colour != 3) & (year(b) == (year-2)))  # Immature females year - 2.
  iy <- which(!is.mature(b) & (b$gonad.colour == 3) & (year(b) == (year-1)) & (b$carapace.width >= 30))  # Pubescent females year - 1.
  iz <- which(is.primiparous(b) & (year(b) == year))  # Primiparous females year.
  
  si <- s[year(s) == year-2, ]
  sp <- s[year(s) == year-1, ]
  sm <- s[year(s) == year, ]
  
  # Function which returns size-frequency variables:
  fvars <- function(x) return(names(x)[gsub("[0-9]", "", names(x)) == ""])
  
  # Attach size-frequencies:
  if (nrow(si) > 0) import(si, fill = 0) = freq(b[ix, ], by = key(s), step = 0.5)
  if (nrow(sp) > 0) import(sp, fill = 0) = freq(b[iy, ], by = key(s), step = 0.5)
  if (nrow(sm) > 0) import(sm, fill = 0) = freq(b[iz, ], by = key(s), step = 0.5)
  
  plot(c(0, 80), c(0, 1000), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i")
  grid()
  if (nrow(si) > 0) lines(as.numeric(fvars(si)), 1000 * apply(si[fvars(si)], 2, mean, na.rm = TRUE), col = "skyblue3", lwd = 2)
  if (nrow(sp) > 0) lines(as.numeric(fvars(sp)), 1000 * apply(sp[fvars(sp)], 2, mean, na.rm = TRUE), col = "orange", lwd = 2)
  if (nrow(sm) > 0) lines(as.numeric(fvars(sm)), 1000 * apply(sm[fvars(sm)], 2, mean, na.rm = TRUE), col = "palegreen3", lwd = 2)
  
  legend("topright", 
         legend = c(paste("Immature", year-2),
                    paste("Pubescent", year-1),
                    paste("Primiparous", year)),
         lwd = 2,
         col = c("skyblue3", "orange", "palegreen3"))
  
  box(col = "grey50")
  
  if (i == length(years)){
     axis(1)
     mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
  } 
  if (i == 1) axis(2) else axis(2, at = seq(0, 800, by = 200))
  if (i == 2) mtext("Density", 2, 2.5, cex = 1.25)
}
dev.off()
       