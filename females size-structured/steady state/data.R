library(gulf.data)
library(gulf.graphics)
library(gulf.spatial)

years <- 1997:2023

# Load data:
b <- read.scsbio(years, valid = 2, survey = "regular", sex = 2, species = 2526)
b$year <- year(b)

# Identify maturity stages:
b$maturity <- maturity(b)

# Calculate size frequencies:
ix <- which((b$maturity == "immature") & (b$gonad.colour == 1) & (b$carapace.width > 5) & (b$carapace.width < 65))
fi <- freq(b[ix, ], by = c("date", "tow.id"))

ix <- which((b$maturity == "immature") & (b$gonad.colour %in% 2:3) & (b$carapace.width > 24) & (b$carapace.width <= 75))
fp <- freq(b[ix, ], by = c("date", "tow.id"))

ix <- which((b$maturity == "mature") & is.new.shell(b) & (b$carapace.width > 32) & (b$carapace.width < 90))
fn <- freq(b[ix, ], by = c("date", "tow.id"))

ix <- which((b$maturity == "mature") & !is.new.shell(b) & (b$carapace.width > 32) & (b$carapace.width < 90))
fo <- freq(b[ix, ], by = c("date", "tow.id"))

# Load tow data and attach size-frequencies:
si <- read.scsset(years, survey = "regular", valid = 1)
import(si, fill = 0) <- fi

sp <- read.scsset(years, survey = "regular", valid = 1)
import(sp, fill = 0) <- fp

sn <- read.scsset(years, survey = "regular", valid = 1)
import(sn, fill = 0) <- fn

so <- read.scsset(years, survey = "regular", valid = 1)
import(so, fill = 0) <- fo

# Buffer frequency variables:
fvars <- function(x) return(names(x)[gsub("[0-9]", "", names(x)) == ""])
fvars <- unique(c(fvars(si), fvars(sp), fvars(sn), fvars(so)))
si[setdiff(fvars, names(si))] <- 0
sp[setdiff(fvars, names(sp))] <- 0
sn[setdiff(fvars, names(sn))] <- 0
so[setdiff(fvars, names(so))] <- 0

# Standardize by swept area:
si[fvars] <- si[fvars] / repvec(si$swept.area, ncol = length(fvars))
sp[fvars] <- sp[fvars] / repvec(sp$swept.area, ncol = length(fvars))
sn[fvars] <- sn[fvars] / repvec(sn$swept.area, ncol = length(fvars))
so[fvars] <- so[fvars] / repvec(so$swept.area, ncol = length(fvars))

# Define sampling grids and survey years:
si$grid <- deg2grid(lon(si), lat(si))
sp$grid <- deg2grid(lon(si), lat(sp))
sn$grid <- deg2grid(lon(si), lat(sn))
so$grid <- deg2grid(lon(si), lat(so))
si$year <- year(si)
sp$year <- year(sp)
sn$year <- year(sn)
so$year <- year(so)

# Spatially disaggregate 10' x 10' grids used before 2012 surveys:
ix <- (year(si) < 2012) 
si <- rbind(aggregate(si[ix,fvars], by = si[ix, c("year", "grid")], mean), si[!ix, c("year", "grid", fvars)])
sp <- rbind(aggregate(sp[ix,fvars], by = sp[ix, c("year", "grid")], mean), sp[!ix, c("year", "grid", fvars)])
sn <- rbind(aggregate(sn[ix,fvars], by = sn[ix, c("year", "grid")], mean), sn[!ix, c("year", "grid", fvars)])
so <- rbind(aggregate(so[ix,fvars], by = so[ix, c("year", "grid")], mean), so[!ix, c("year", "grid", fvars)])

# Summary frequency tables:
ti <- aggregate(si[fvars], by = list(year = year(si)), mean)
tp <- aggregate(sp[fvars], by = list(year = year(si)), mean)
tn <- aggregate(sn[fvars], by = list(year = year(si)), mean)
to <- aggregate(so[fvars], by = list(year = year(si)), mean)

# Aggregated surgvey summaries:
mi <- 1000000 * apply(ti[fvars], 2, mean)
mp <- 1000000 * apply(tp[fvars], 2, mean)
mn <- 1000000 * apply(tn[fvars], 2, mean)
mo <- 1000000 * apply(to[fvars], 2, mean)

png(filename = "females size-structured/figures/female SF by maturity.png", 
    width = 7.5, height = 7.5, units = "in", res = 300)

# Abundance figure:
cols <- c("deepskyblue", "green3", "royalblue", "brown3")
plot(c(0, 90), c(0, 400), xaxs = "i", yaxs = "i", type = "n", xlab = "", ylab = "")
grid()
lines(as.numeric(fvars), mi, col = cols[1], lwd = 2)
lines(as.numeric(fvars), mp, col = cols[2], lwd = 2)
lines(as.numeric(fvars), mn, col = cols[3], lwd = 2)
lines(as.numeric(fvars), mo, col = cols[4], lwd = 2)

# Instar labels:
mu <- c(9.9, 14.5, 20.5, 29, 37, 49, 58)
names(mu) <- as.roman(4:10)
text(mu[1:5], approx(as.numeric(fvars), mi, mu[1:5])$y, names(mu)[1:5], pos = 3, font = 2)
text(mu[6], approx(as.numeric(fvars), mp, mu[6])$y, names(mu)[6], pos = 3, font = 2)
text(mu[7], approx(as.numeric(fvars), mn, mu[7])$y, names(mu)[7], pos = 3, font = 2)

mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
mtext(expression(paste("Average density (# / km"^"2",")", sep="")), 2, 2.25, cex = 1.25, font = 2)

legend("topleft", 
       legend = c("Immature", "Pubescent", "Primiparous", "Multiparous"),
       lwd = 2,
       col = cols,
       cex = 1.25)

box(col = "grey50")

dev.off()



plot(c(0, 90), c(0, 400), xaxs = "i", yaxs = "i", type = "n", xlab = "", ylab = "")
grid()
for (i in 1:nrow(ti)){
   tmp <- 1000000 * tn[i,fvars]
   tmp <- 4000 * tmp / sum(tmp)
   lines(as.numeric(fvars), tmp)
}
  
clg()
png(filename = "females size-structured/figures/female SF by maturity and year.png", 
    width = 8.5, height = 11, units = "in", res = 300)
m <- kronecker(matrix(1:28, ncol = 4), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0, cbind(0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))

for (i in 3:nrow(ti)){
  plot(c(10, 80), c(0, 400), xaxs = "i", yaxs = "i", type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n")
  grid()
  lines(as.numeric(fvars), 1000000 * ti[i-2,fvars], col = "deepskyblue", lwd = 2)
  lines(as.numeric(fvars), 1000000 * tp[i-1,fvars], col = "green3", lwd = 2)
  lines(as.numeric(fvars), 1000000 * tn[i,fvars], col = "royalblue", lwd = 2)
  text(par("usr")[1] + 0.5 * diff(par("usr")[1:2]), par("usr")[3] + 0.85 * diff(par("usr")[3:4]), paste0(years[i-2], "-", years[i]))
  
  if (i %in% c(9, 16, 23, nrow(ti))) axis(1)
  
  box(col = "grey50")
}

dev.off()




  


