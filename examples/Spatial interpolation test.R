library(fields)
library(raster)
library(glmmTMB)
library(ggplot2)
library(geoR)
library(viridis)
library(tidyverse)
library(gridExtra)
library(NLMR)
library(DHARMa)

data(ca20)
dat <- data.frame(x = ca20$coords[,1], y = ca20$coords[,2], calcium = ca20$data, elevation = ca20$covariate[,1], region = factor(ca20$covariate[,2]))
dat$pos <- numFactor(scale(dat$x), scale(dat$y))
dat$ID <- factor(rep(1, nrow(dat)))

m_tmb <- glmmTMB(calcium ~ elevation + region + mat(pos + 0 | ID), dat) # takes some time to fit

summary(m_tmb)

sims <- simulateResiduals(m_tmb)
plot(sims)

# some R magic to extract and re-order the estimated correlation between pairs of locations
fit_cor <- matrix(as.numeric(attr(VarCorr(m_tmb)$cond$ID, "correlation")), nrow = 178, ncol = 178, byrow = FALSE, 
                  dimnames = attr(attr(VarCorr(m_tmb)$cond$ID, "correlation"),"dimnames"))

ff <- dimnames(fit_cor)[[1]]
ff <- gsub("pos","",ff)
fit_cor2 <- fit_cor[order(match(ff, dat$pos)), order(match(ff, dat$pos))]

# plot
plot(as.numeric(dd), fit_cor2[lower.tri(fit_cor2)],
     xlab = "Distance between pairs of location [m]",
     ylab = "Estimated correlation")

# the effect of elevation
newdat <- data.frame(elevation = seq(3, 7, length = 10), region = factor(1, levels = 1:3))
# turn this into a model matrix
mm <- model.matrix(~ elevation + region, newdat)
newdat$calcium <- mm %*% fixef(m_tmb)$cond + mean(c(0, fixef(m_tmb)$cond[3:4])) # predicted values removing region effects
pvar <- diag(mm %*% tcrossprod(vcov(m_tmb)$cond, mm))
newdat$lci <- newdat$calcium - 1.96 * sqrt(pvar)
newdat$uci <- newdat$calcium + 1.96 * sqrt(pvar)

gg1 <- ggplot(dat, aes(x = elevation, y = calcium)) +
   geom_point() +
   geom_line(data = newdat) +
   geom_ribbon(data = newdat, aes(ymin = lci, ymax = uci), alpha = 0.2)

# the effect of region
newdat <- data.frame(elevation = mean(dat$elevation), region = factor(1:3))
# turn this into a model matrix
mm <- model.matrix(~ elevation + region, newdat)
newdat$calcium <- mm %*% fixef(m_tmb)$cond # predicted values 
pvar <- diag(mm %*% tcrossprod(vcov(m_tmb)$cond, mm))
newdat$lci <- newdat$calcium - 1.96 * sqrt(pvar)
newdat$uci <- newdat$calcium + 1.96 * sqrt(pvar)

gg2 <- ggplot(dat, aes(x = region, y = calcium)) +
   geom_jitter() +
   geom_point(data = newdat, color = "red", size = 2) +
   geom_linerange(data = newdat, aes(ymin = lci, ymax = uci), color = "red")

# plot together
grid.arrange(gg1, gg2, ncol = 2)



# predict at any location
# derive a DEM
elev_m <- Tps(dat[,c("x","y")], dat$elevation)

newdat <- expand.grid(x = seq(4960, 5960, length.out = 50), y = seq(4830, 5710, length.out = 50))
newdat$ID <- factor(rep(1, nrow(newdat)))
newdat$elevation <- extract(elev, newdat[,1:2])
newdat$region <- factor(extract(region, newdat[,1:2]))
# remove NAs
newdat <- na.omit(newdat)
newdat$pos <- numFactor(((newdat$x - mean(dat$x)) / sd(dat$x)), ((newdat$y - mean(dat$y)) / sd(dat$y)))
# predict in slices of 100 predictions to speed up computation
pp <- rep(NA, 1927)
for(i in seq(1, 1927, by = 100)){
   if(i == 1901){
      pp[1901:1927] <- predict(m_tmb, newdat[1901:1927,], allow.new.levels = TRUE)
   }
   else{
      pp[i:(i+99)] <- predict(m_tmb, newdat[i:(i+99),], allow.new.levels = TRUE)
   }
   # print(i)
}
newdat$calcium <- pp
(gg_tmb <- ggplot(newdat,aes(x=x, y=y, fill = calcium)) +
      geom_raster() +
      scale_fill_viridis())





