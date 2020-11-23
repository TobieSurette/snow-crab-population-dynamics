library(gulf.data)
library(gulf.spatial)
library(gulf.graphics)

x <- read.nssset(2018, experiment = 1)

x$longitude.start <- -dmm2deg(x$longitude.start)
x$longitude.stop <- -dmm2deg(x$longitude.stop)
x$latitude.start <- dmm2deg(x$latitude.start)
x$latitude.stop <- dmm2deg(x$latitude.stop)

which(round(depth(longitude(x), latitude(x))) <= 0)

clg()
map.new(xlim = c(-65.5, -61), ylim = c(45.5, 47.25))
map("bathymetry")
map("coast")
points(longitude(x), latitude(x))
box()
axis(1)
axis(2)

# Read lobster length observations:
y <- read.nsslen(2018:2019, species = 2550)
y$sex[which(y$year == 2020 & y$species == 2550 & is.na(y$egg.condition))] <- 1
y$sex[which(y$year == 2020 & y$species == 2550 & !is.na(y$egg.condition))] <- 2
y$gear <- x$gear[match(y$set.number, x$set.number)]

clg()
dev.new(height = 11, width = 8.5)
m <- kronecker(1:2, matrix(1, nrow = 5, ncol = 5))
m <- rbind(0, cbind(0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
years <- sort(unique(year(y)))
   
plot(c(0, 120), c(0, 180), xaxs = "i", yaxs = "i", xaxt = "n")
grid()
gbarplot(table(y$length[y$sex == 1 & year(y) == years[1]]), width = 1, border = "grey50", add = TRUE)

plot(c(0, 120), c(0, 180), xaxs = "i", yaxs = "i", xaxt = "n")
grid()
gbarplot(table(y$length[y$sex == 1 & year(y) == years[2]]), width = 1, border = "grey50", add = TRUE)

z <- round(grow(y$length[y$sex == 1 & year(y) == years[1]], species = "lobster"))
t <- table(z)
lines(as.numeric(names(t)), t, lwd = 2, col = "blue")

axis(1)
mtext("Carapace length(mm)", 1, 2.5, cex = 1.25)



