# Program to calculate and map out variations in mean sizes:
years <- 1998:2022
language <- language("en")
b <- read.scsbio(years)

# Target variable:
b <- b[which(is.pubescent.scsbio(b)), ] # Pubescent females. 
map.file       <- "Pubescent female size anomaly maps.png"
mean.size.file <- paste0("Pubescent female sizes by year ", language, ".png")

clg()
png(file = paste0(map.file), res = 500, units = "in", height = 7.0, width = 10)
m <- kronecker(matrix(1:25, ncol = 5), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0, cbind(0, m, 0), 0)
layout(m)
par(mar = c(0, 0, 0, 0))
for (i in 1:length(years)){
   print(i)
   
   set <- read.scsset(years[i], valid = 1, survey = "regular")
   set$depth <- depth(lon(set), lat(set))
   #set$temperature <- temperature(years[i]-1, depth = "bottom", longitude = lon(set), latitude = lat(set))
   
   bio <- b[year(b) == years[i], ]
   
   res <- aggregate(list(mu = bio$carapace.width), by = bio[key(set)], mean, na.rm = TRUE)
   res$sigma <- aggregate(list(x = bio$carapace.width), by = bio[key(set)], sd, na.rm = TRUE)$x
   
   import(set) <- res
   import(set, fill = 0) <- catch(bio)
   set$sigma <- set$sigma / sqrt(set$total)
   
   r <- set$mu - mean(set$mu, na.rm = TRUE)
   #plot(set$temperature, r, xlim = c(-1, 7), ylim = c(-10, 10))
   #grid()
   #title(years)
   
   map.new()
   tmp <- drop(temperature(years[i]-1, depth = "bottom"))
   breaks <- seq(-1, 10, len = 101)
   cols <- c(colorRampPalette(c("blue", "white"))(which.min(abs(breaks-3))), colorRampPalette(c("white", "red"))(length(breaks)-which.min(abs(breaks-3))-1))
   image(as.numeric(dimnames(tmp)$longitude), as.numeric(dimnames(tmp)$latitude), tmp, 
         xlab = "", ylab = "", xaxt = "n", yaxt = "n", col = cols, zlim = c(0, 10), breaks = breaks, add = TRUE)
   
   #map("bathymetry")
   
   ix <- which(r >= 0)
   points(lon(set)[ix], lat(set)[ix], pch = 21, bg = "grey70", col = "grey40", cex = 0.25* sqrt(abs(r)[ix]), lwd = 0.5)
   ix <- which(r < 0)
   points(lon(set)[ix], lat(set)[ix], pch = 21, bg = "brown1", col = "brown3", cex = 0.25 * sqrt(abs(r)[ix]), lwd = 0.5)
   map("coast")
   
   # Display survey year:      
   text(par("usr")[1] + 0.8 * diff(par("usr")[1:2]),
        par("usr")[3] + 0.8 * diff(par("usr")[3:4]),
        years[i], cex = 1.25, pos = 3)
   
   if (i == 5) map.axis(1:2)
   if (i == 21) map.axis(3:4)
   
   # Draw box:
   box(col = "grey50")
}
dev.off()

# Mean sizes
for (i in 1:length(years)){
   print(years[i])
   
   # Load set data:
   s <- read.scsset(year = years[i], valid = 1, survey = "regular")
   s$grid <- deg2grid(lon(s), lat(s))
   
   # Target group:
   bio <- b[year(b) == years[i], ]
   
   # Length-frequency from biolgical data:
   f <- freq(bio, by = c("date", "tow.id"), step = 0.5)
   
   # Attatch catch data to set card:
   import(s, fill = 0) <- f
   fvars <- names(s)[gsub("[-0-9.]", "", names(s)) == ""]
   s[, fvars] <- 1000000 * s[, fvars] / repvec(s$swept.area, ncol = length(fvars))
   
   # Disaggregate spatial data for early years:
   if (years[i] <= 2012) s <- aggregate(s[fvars], by = s["grid"], mean)
   
   # Calculate size quantiles:
   f <- apply(s[fvars], 2, mean)
   m <- sum(as.numeric(names(f)) * (f / sum(f)))
   f <- cumsum(f) / sum(f)
   p <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.90, 0.95, 0.975)
   q <- c(m, approx(f, as.numeric(names(f)), p)$y)
   
   if (i == 1) res <- q else res <- rbind(res, q)
}
rownames(res) <- years
colnames(res) <- c("mean", paste0("p", 100 * p))
rf <- as.data.frame(res)
rf$year <- years

# Mean sizes:
clg()
png(file = mean.size.file, units = "in", res = 500, height = 7, width = 7)

plot(range(rf$year), c(30, 70), type = "n", xlab = "", ylab = "", cex.lab = 1.5, yaxs = "i", xaxt = "n", yaxt = "n")
grid()
w <- 0.4
for (i in 1:nrow(rf)){
   rect(rf$year[i] - w, rf$p25[i], rf$year[i] + w, rf$p75[i], col = "grey70", border = "grey30", lwd = 0.5)
   lines(rf$year[i] + c(-w,w), rep(rf$mean[i], 2), lwd = 1.5)
   lines(rep(rf$year[i], 2), c(rf$p2.5[i], rf$p25[i]), lwd = 1.0, lty = "dashed")
   lines(rep(rf$year[i], 2), c(rf$p75[i], rf$p97.5[i]), lwd = 1.0, lty = "dashed")
   lines(rf$year[i] + c(-w/2,w/2), rep(rf$p2.5[i], 2), lwd = 1.0)
   lines(rf$year[i] + c(-w/2,w/2), rep(rf$p97.5[i], 2), lwd = 1.0)
}
axis(2, at = seq(30, 70, by = 5))
box(col = "grey30")

if (language == "english"){
   mtext("Year", 1, 2.5, cex = 1.25)
   mtext("Carapace width (mm)", 2, 2.5, cex = 1.25)
   #text(par("usr")[1] + 0.9*diff(par("usr")[1:2]), par("usr")[4] - 2, "Female", cex = 1.5)
}else{
   mtext("Année", 1, 2.5, cex = 1.25)
   mtext("Largeur de carapace (mm)", 2, 2.5, cex = 1.25)
   #text(par("usr")[1] + 0.9*diff(par("usr")[1:2]), par("usr")[4] - 2, "Femelle", cex = 1.5)
}
hline(mean(rf$mean), lwd = 2, lty = "dashed", col = "red")
axis(1)
box(col = "grey50")

dev.off()


# Classify white gonads as either 



clg()
plot(log(b$carapace.width), col = c("grey40", "orange")[v+1], cex = b$cex, ylim = c(2, 4.5), yaxs = "i")
grid()

hist(b$carapace.width[v == 0], n = 500)

b$pubescent <- v

tmp <- aggregate(b$pubescent, by = list(cw = round(b$carapace.width)), function(x) return(sum(x, na.rm = TRUE)/sum(!is.na(x))))
gbarplot(logit(tmp[, 2]), tmp[, 1], xlim = c(0, 70), ylim = c(0, 1))
title(main = year)

sum(is.na(b$pubescent))
b$year <- year(b)
b$yearf <- as.factor(b$year)

y <- b[which(b$year == 2022 & !is.na(b$pubescent)), ]

m <- list()
m[[1]] <- glmmTMB(pubescent ~ carapace.width, data = b)
m[[2]] <- glmmTMB(pubescent ~ yearf + carapace.width, data = b)
m[[3]] <- glmmTMB(pubescent ~ (1|yearf) + carapace.width, data = b)
m[[4]] <- glmmTMB(pubescent ~ 0 + (carapace.width | yearf), data = b)
m[[5]] <- glmmTMB(pubescent ~ 0 + (carapace.width | tow.id), data = y)
m[[6]] <- glmmTMB(pubescent ~ 0 + (carapace.width | station), data = b)

years 

b$station <- paste0(b$year, "-", b$tow.id)

a <- ranef(m[[6]])$station
a$x50 <- -a[,1] / a[,2]

a$year <- as.numeric(unlist(lapply(strsplit(rownames(a), "-"), function(x) x[1])))
a$tow.id <- unlist(lapply(strsplit(rownames(a), "-"), function(x) x[2]))

s <- read.scsset(years)
s$year <- year(s)
ix <- match(a[c("year", "tow.id")], s[c("year", "tow.id")])
a$longitude <- lon(s[ix,])
a$latitude <- lat(s[ix,])

map.new()
ix <- a$x50 > mean(a$x50, na.rm = TRUE)
points(a$longitude[ix], a$latitude[ix], pch = 21, bg = "grey50", cex = 0.5 * sqrt(abs(a$x50[ix] - mean(a$x50, na.rm = TRUE))))
ix <- a$x50 < mean(a$x50, na.rm = TRUE)
points(a$longitude[ix], a$latitude[ix], pch = 21, bg = "brown2", cex = 0.5 * sqrt(abs(a$x50 - mean(a$x50, na.rm = TRUE))))

clg()
map.new()
for (i in 1:length(years)){
   
   aa <- a[a$year == years[i] & !is.na(a$x50) & a$x50 > 15 & a$x50 < 55, ]
   v <- aa$x50
   ix <- v > mean(v)
   points(aa$longitude[ix], aa$latitude[ix], pch = 21, bg = "grey50", cex = 0.8 * sqrt(abs(v[ix] - mean(v))))
   ix <- v < mean(v)
   points(aa$longitude[ix], aa$latitude[ix], pch = 21, bg = "brown2", cex = 0.8 * sqrt(abs(v[ix] - mean(v))))
   map("coast")
   title(main = years[i])
}

# Mean sizes of pubescent females:
ix <- b$pubescent == 1
res <- aggregate(b$carapace.width[ix], by = b[ix, "year", drop = FALSE], mean, na.rm = TRUE)

ix <- b$pubescent == 1
res <- aggregate(b$carapace.width[ix], by = b[ix, "station", drop = FALSE], mean, na.rm = TRUE)
res$year <- as.numeric(unlist(lapply(strsplit(res$station, "-"), function(x) x[1])))
res$tow.id <- unlist(lapply(strsplit(res$station, "-"), function(x) x[2]))
res$n <- aggregate(b$carapace.width[ix], by = b[ix, "station", drop = FALSE], length)$x
s <- read.scsset(years)
s$year <- year(s)
ix <- match(res[c("year", "tow.id")], s[c("year", "tow.id")])
res$longitude <- lon(s[ix,])
res$latitude <- lat(s[ix,])

clg()
for (i in 1:length(years)){
   map.new()
   r <- res[which(res$year == years[i] & res$n > 0), ]
   ix <- r$x > mean(r$x)
   points(r$longitude[ix], r$latitude[ix], pch = 21, bg = "grey50", cex = 0.8 * sqrt(abs(r$x[ix] - mean(r$x, na.rm = TRUE))))
   ix <- r$x < mean(r$x)
   points(r$longitude[ix], r$latitude[ix], pch = 21, bg = "brown2", cex = 0.8 * sqrt(abs(r$x[ix] - mean(r$x, na.rm = TRUE))))
   map("coast")
   title(main = years[i])
}



boxplot(a$x50 ~ a$year, ylim = c(20, 30))


