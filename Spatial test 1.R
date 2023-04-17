library(gulf.data)
library(gulf.graphics)
library(gulf.spatial)
library(glmmTMB)
library(splines)

# Load survey polygon:
p <- read.gulf.spatial("kriging polygons revised")$gulf
p$x <- p$longitude
p$y <- p$latitude
p <- as.polygon(list(p))

# Load data:
x <- read.scsset(year = 2021, survey = "regular", valid = 1)
b <- read.scsbio(2021)
import(x, fill = 0) <- catch(b, category = "COM")
names(x) <- gsub("COM", "n", names(x))
x$depth <- round(depth(lon(x), lat(x)),1)
x$log.depth <- log(x$depth)
x$fixed <- as.numeric(x$station.type == "fixed")

# Assemble data:
data <- cbind(x, deg2km(lon(x), lat(x)))
data$x <- round(data$x)
data$y <- round(data$y)
data$pos <- numFactor(data$x, data$y)
data$group <- factor(rep(1, nrow(data)))
data$fixed <- x$fixed
data$fixed <- x$fixed

# Fit model:
model <- glmmTMB(n ~ 1 + fixed + log.depth + I(log.depth^2) + exp(pos + 0 | group) + offset(log(swept.area)), 
                 data = data, family = nbinom2)

# Predictions:
xx <- seq(min(data$x), max(data$x), by = 2)
yy <- seq(min(data$y), max(data$y), by = 2)
newdata <- data.frame(pos = numFactor(expand.grid(x = xx, y = yy)))
newdata$x <- as.numeric(unlist(lapply(strsplit(as.character(newdata$pos), "[(),]"), function(x) x[2])))
newdata$y <- as.numeric(unlist(lapply(strsplit(as.character(newdata$pos), "[(),]"), function(x) x[3])))
newdata$group <- factor(rep(1, nrow(newdata)))
newdata <- cbind(newdata, km2deg(newdata$x, newdata$y))
newdata$depth <- round(depth(newdata$longitude, newdata$latitude), 1)
newdata$log.depth <- log(newdata$depth)
newdata$swept.area <- 1
newdata$n <- NA
ix <- which(in.polygon(p, newdata$longitude, newdata$latitude))

stop <- FALSE
while (!stop){
  index <- ix[1:min(100, length(ix))]
  print(index)
  if (length(index) > 0){
     newdata$n[index] <- predict(model, newdata[index, ], type = "response", allow.new.levels = TRUE)
     ix <- setdiff(ix, index)
     print(paste(length(ix), "remaining"))
  }else{
     stop <- TRUE
  }
}

n <- newdata$n
dim(n) <- c(length(xx), length(yy))

image(n, main = "Reconstruction", col = rev(terrain.colors(100)))

