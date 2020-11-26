library(gulf.data)
library(gulf.spatial)
library(gulf.graphics)

x <- read.nssset(2020, experiment = 1)

x$longitude.start <- -dmm2deg(x$longitude.start)
x$longitude.stop <- -dmm2deg(x$longitude.stop)
x$latitude.start <- dmm2deg(x$latitude.start)
x$latitude.stop <- dmm2deg(x$latitude.stop)

# Check for coordinate errors:
x[which(round(depth(x$longitude.start, x$latitude.start)) <= 0), ]
x[which(round(depth(x$longitude.stop, x$latitude.stop)) <= 0), ]

map.new(xlim = c(-65.5, -63), ylim = c(45.75, 47.25))
map("bathymetry")
map("coast")
points(longitude(x), latitude(x))
box()
axis(1)
axis(2)

# Read lobster length observations:
y <- read.nsslen(2020, species = 2550)
y$year <- year(y)
y$sex[which(y$year == 2020 & y$species == 2550 & is.na(y$egg.condition))] <- 1
y$sex[which(y$year == 2020 & y$species == 2550 & !is.na(y$egg.condition))] <- 2
y$gear <- x$gear[match(y$set.number, x$set.number)]

clg()
dev.new(height = 11, width = 8.5)
m <- kronecker(1:4, matrix(1, nrow = 5, ncol = 5))
m <- rbind(0, cbind(0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
for (gear in c(13, 19)){
   for (sex in 1:2){
      plot(c(0, 120), c(0, 140), xaxs = "i", yaxs = "i", xaxt = "n")
      grid()
      gbarplot(table(y$length[y$sex == sex & y$gear == gear]), width = 1, border = "grey50", add = TRUE)
      text(10, 120, c("Male", "Female")[sex])
      text(100, 120, gear(gear))
   }
}
axis(1)
mtext("Carapace length(mm)", 1, 2, cex = 1.25)

import(x, by = "set.number") <- freq(y, by = "set.number")


