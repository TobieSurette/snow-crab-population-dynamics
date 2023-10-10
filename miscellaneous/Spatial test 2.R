library(glmmTMB)

d <- data.frame(z = as.vector(volcano),
                x = as.vector(row(volcano)),
                y = as.vector(col(volcano)))

set.seed(1)
d$z <- d$z + rnorm(length(volcano), sd=15)
d <- d[sample(nrow(d), 100), ]


volcano.data <- array(NA, dim(volcano))
volcano.data[cbind(d$x, d$y)] <- d$z
image(volcano.data, main="Spatial data", useRaster=TRUE)

d$pos <- numFactor(d$x, d$y)
d$group <- factor(rep(1, nrow(d)))

f <- glmmTMB(z ~ 1 + exp(pos + 0 | group), data = d)
confint(f, "sigma")

newdata <- data.frame( pos=numFactor(expand.grid(x=1:3,y=1:3)) )
newdata$group <- factor(rep(1, nrow(newdata)))

predict(f, newdata, type = "response", allow.new.levels = TRUE)

predict_col <- function(i) {
    newdata <- data.frame(pos = numFactor(expand.grid(1:87,i)))
    newdata$group <- factor(rep(1,nrow(newdata)))
    predict(f, newdata = newdata, type = "response", allow.new.levels = TRUE)
}

pred <- sapply(1:61, predict_col)

image(pred, main = "Reconstruction")








