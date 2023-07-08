library(gulf.data)
library(gulf.graphics)
library(RColorBrewer)
library(TMB)

# Compile TMB model:
clc()
compile("instar/instar.cpp")
dyn.load(dynlib("instar/instar"))

# Prepare data:
years <- 1990:2022
s <- read.scsset(years, valid = 1, survey = "regular")
s$grid <- grid.scs(lon(s), lat(s), correct = FALSE)

b <- read.scsbio(years)
b <- b[!is.na(match(b[key(s)], s[key(s)])), ]
b$year <- year(b)
b <- b[which(b$sex == 2), ]
ix <- match(b[key(s)], s[key(s)])

b$grid <- s$grid[ix]


# recommendation: Group by grid for testing.

# Define maturity stages:
b$maturity <- ""
b$maturity[which(!is.na(b$carapace.width) & (b$sex == 2) & (b$carapace.width >= 9) & (b$carapace.width <= 70) & !is.mature(b) & !is.pubescent.scsbio(b))] <- "immature"
b$maturity[which(!is.na(b$carapace.width) & (b$sex == 2) & (b$carapace.width >= 28) & (b$carapace.width <= 80) & !is.mature(b) & is.pubescent.scsbio(b))] <- "pubescent"
b$maturity[which(!is.na(b$carapace.width) & (b$sex == 2) & (b$carapace.width >= 36) & (b$carapace.width <= 95) & is.primiparous.scsbio(b))] <- "mature"
b <- b[b$maturity != "", ]

# Define sampling group:
b$group <- match(b[key(s)], unique(b[key(s)]))-1

# Sub-sample for testing:
#tmp <- aggregate(b$carapace.width, by = b["group"], length)
#b <- b[b$group %in% tmp[rev(order(tmp$x))[1:500], "group"], ]
#b$group <- match(b$group, sort(unique(b$group)))

tab <- NULL
for (i in 1:length(years)){
   if (years[i] < 1998) precision <- 1 else precision <- 0.1
   bb <- b[b$year == years[i]  , ]
   bb$cw <- round(bb$carapace.width, -log10(precision))
   tmp <- aggregate(list(n = bb$carapace.width), by = bb[c("year", "tow.id", "group", "maturity", "cw")], length)
   tmp$precision <- precision
   tab <- rbind(tab, tmp)
}

# Define data:
data <- list(x_imm = tab$cw[tab$maturity == "immature"], 
             x_pub = tab$cw[tab$maturity == "pubescent"],  
             x_mat = tab$cw[tab$maturity == "mature"], 
             f_imm = tab$n[tab$maturity == "immature"],
             f_pub = tab$n[tab$maturity == "pubescent"], 
             f_mat = tab$n[tab$maturity == "mature"],
             precision_imm = tab$precision[tab$maturity == "immature"],
             precision_pub = tab$precision[tab$maturity == "pubescent"], 
             precision_mat = tab$precision[tab$maturity == "mature"],
             group_imm = tab$group[tab$maturity == "immature"],
             group_pub = tab$group[tab$maturity == "pubescent"], 
             group_mat = tab$group[tab$maturity == "mature"],
             n_instar = 11,
             n_group = max(tab$group)+1)

# Define parameters:
parameters <- list(mu_imm_0 = 0.908,              # Size of immature first instar.
                   a_imm    = 0.891,              # Slope parameter for immature growth.
                   b_imm    = 0.635,              # Intercept parameter for immature growth.
                   a_pub    = 0.642,              # Slope parameter for pubescent growth.
                   b_pub    = 1.58,               # Intercept parameter for pubescent growth.
                   a_mat    = 1.03,               # Slope parameter for mature growth.
                   b_mat    = 0.0405,             # Intercept parameter for mature growth.   
                   log_sigma_delta_mu = -3,
                   delta_mu_imm = rep(0, data$n_instar * data$n_group),
                   delta_mu_pub = rep(0, data$n_instar * data$n_group),
                   delta_mu_mat = rep(0, data$n_instar * data$n_group),
                   log_sigma = -2.76,
                   logit_p_imm_global = rep(0,6),
                   logit_p_pub_global = rep(0,2),
                   logit_p_mat_global = rep(0,2),
                   delta_logit_p_imm = matrix(0.1*runif(data$n_group * 6)-0.05, nrow = data$n_group, ncol = 6),
                   delta_logit_p_pub = matrix(0.1*runif(data$n_group * 2)-0.05, nrow = data$n_group, ncol = 2),
                   delta_logit_p_mat = matrix(0.1*runif(data$n_group * 2)-0.05, nrow = data$n_group, ncol = 2),
                   log_sigma_delta_logit_p = 0.66) 

# Define functions:
free <- function(x, p){
   if (missing(p)) p <- names(x)
   p <- p[which(p %in% names(x))]
   if (length(p) > 0){
      for (i in 1:length(p)){
         d <- dim(x[[p[i]]])
         x[[p[i]]] <- as.factor(1:length(x[[p[i]]]))
         dim(x[[p[i]]]) <- d
      } 
   }
   return(x)
}

fix <- function(x, p){
   if (missing(p)) p <- names(x)
   p <- p[which(p %in% names(x))]
   if (length(p) > 0){
      for (i in 1:length(p)){
         d <- dim(x[[p[i]]])
         x[[p[i]]] <- as.factor(NA * x[[p[i]]])
         dim(x[[p[i]]]) <- d
      } 
   }
   return(x)
}

update <- function(x, p){
   if (!missing(p)){
       str <- unique(names(p))
       for (i in 1:length(str)){
          d <- dim(x[[str[i]]])
          x[[str[i]]] <- p[names(p) == str[i]]
          dim(x[[str[i]]]) <- d
       }
   }
  return(x)
}
map <- function(x){
   v <- lapply(x, function(x) as.factor(NA * x))
   for (i in 1:length(x)) dim(v[[i]]) <- dim(x[[i]])
   return(v)
}
  
map <- map(parameters)

# Define set of random variables:
random <- c("delta_mu_imm", "delta_mu_pub", "delta_mu_mat", "delta_logit_p_imm", "delta_logit_p_pub", "delta_logit_p_mat")

# Fit global proportions:
map <- free(map, c("logit_p_imm_global", "logit_p_pub_global", "logit_p_mat_global"))
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "instar")
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 500))$par
parameters <- update(parameters, obj$par)

# Fit instar standard error:
map$log_sigma <- as.factor(1)
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "instar")
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
parameters <- update(parameters, obj$par)

# Fit global growth parameters:
map <- free(map, c("mu_imm_0", "a_imm", "b_imm", "a_pub", "b_pub", "a_mat", "b_mat"))
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "instar")
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
parameters <- update(parameters, obj$par)

# Fit group proportions:
map <- free(map, c("delta_logit_p_imm", "delta_logit_p_pub", "delta_logit_p_mat", "log_sigma_delta_logit_p"))
obj <- MakeADFun(data = data, parameters = parameters, random = random, map = map, DLL = "instar")
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 500))$par
parameters <- update(parameters, obj$par)
rep <- sdreport(obj)
parameters <- update(parameters, summary(rep, "random")[,1])

# Keep global growth parameters constant, fit delta_mu random effects:
# Write 'fix' function:

save.image(file = paste0("instar/results ", min(years), "-", max(years), ".rdata"))

# Fix global growth parameters and allow for group-level instar size variations:
map <- fix(map, c("mu_imm_0", "a_imm", "b_imm", "a_pub", "b_pub", "a_mat", "b_mat"))
map <- free(map, c("delta_mu_imm", "delta_mu_pub", "delta_mu_mat", "log_sigma_delta_mu"))
parameters$log_sigma_delta_mu <- -3
obj <- MakeADFun(data = data, parameters = parameters, random = random, DLL = "instar")
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 500))$par
parameters <- update(parameters, obj$par)
rep <- sdreport(obj)
parameters <- update(parameters, summary(rep, "random")[,1])

# Instar classification probabilities:
r <- obj$report()

save.image(file = paste0("instar/results 2 ", min(years), "-", max(years), ".rdata"))

# Assign instar membership probabilities:
b$instar <- NA
for (i in 1:nrow(r$mu_imm)){
   if ((i %% 50) == 0) print(i)
   # Immatures:
   pi <- r$p_imm[i,]
   mi <- r$mu_imm[i,]
   si <- r$sigma_imm[i,]
   ix <- which((b$maturity == "immature") & (b$group == (i-1)))
   if (length(ix) > 0){
      for (j in 1:length(ix)){
         u <- pi * dnorm(log(b$carapace.width[ix[j]]), mi, si)
         u <- u / sum(u)
         b$instar[ix[j]] <- which(as.numeric(rmultinom(1, size = 1, prob = u / sum(u))) == 1)
      }
   }
   
   # Pubescents:
   pp <- r$p_pub[i,]
   mp <- r$mu_pub[i,]
   sp <- r$sigma_pub[i,]
   ix <- which((b$maturity == "pubescent") & (b$group == (i-1)))
   if (length(ix) > 0){
      for (j in 1:length(ix)){
         u <- pp * dnorm(log(b$carapace.width[ix[j]]), mp, sp)
         u <- u / sum(u)
         b$instar[ix[j]] <- which(as.numeric(rmultinom(1, size = 1, prob = u / sum(u))) == 1)
      }
   }
   
   # Matures:
   pm <- r$p_mat[i,]
   mm <- r$mu_mat[i,]
   sm <- r$sigma_mat[i,]
   ix <- which((b$maturity == "mature") & (b$group == (i-1)))
   if (length(ix) > 0){
      for (j in 1:length(ix)){
         u <- pm * dnorm(log(b$carapace.width[ix[j]]), mm, sm)
         u <- u / sum(u)
         b$instar[ix[j]] <- which(as.numeric(rmultinom(1, size = 1, prob = u / sum(u))) == 1)
      }
   }
}

clg()
png(file = paste0("instar/Immature classification ", min(years), "-", max(years), ".png"), res = 500, units = "in", height = 8.5, width = 11)
par(mar = c(5, 4, 4, 4) + 0.1) # c(bottom, left, top, right)
plot(c(0, nrow(b)), c(7, 55), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i")
ix <- which(b$maturity == "immature")
ix <- sample(ix)
instars <- sort(unique(b$instar[ix]))
#cols <- rainbow(length(instars))
cols <- brewer.pal(11, name = "Paired")
yy <- sort(unique(b$year))
at <- NULL
for (i in 1:length(yy)){
   xlim <- range(which(b$year == yy[i]))
   rect(xlim[1], par("usr")[3], xlim[2], par("usr")[4], col = c("grey80", "grey95")[((i%%2)==0)+1], border = "grey80")
   at <- c(at, mean(xlim))
}
points(ix, b$carapace.width[ix], col = cols[b$instar[ix]], cex = 0.1)
for (i in 1:length(instars)) hline(exp(r$mu_imm_global[instars[i]]), col = cols[instars[i]], lwd = 2)
axis(1, labels = yy, at = at, las = 2)
axis(2)
axis(4, labels = as.roman(1:length(r$mu_imm_global)), at = exp(r$mu_imm_global), font = 2, padj = -0.5)
mtext("Carapace width (mm)", 2, 2.5, cex = 1.25)
mtext("Year", 1, 3.5, cex = 1.25)
mtext("Instar", 4, 2.0, cex = 1.25)
box(col = "grey60")
dev.off()

clg()
png(file = paste0("instar/Pubescent classification ", min(years), "-", max(years), ".png"), res = 500, units = "in", height = 8.5, width = 11)
par(mar = c(5, 4, 4, 4) + 0.1) # c(bottom, left, top, right)
plot(c(0, nrow(b)), c(30, 70), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i")
ix <- which(b$maturity == "pubescent")
ix <- sample(ix)
instars <- sort(unique(b$instar[ix]))
#cols <- rainbow(length(instars))
cols <- brewer.pal(11, name = "Paired")
yy <- sort(unique(b$year))
at <- NULL
for (i in 1:length(yy)){
   xlim <- range(which(b$year == yy[i]))
   rect(xlim[1], par("usr")[3], xlim[2], par("usr")[4], col = c("grey80", "grey95")[((i%%2)==0)+1], border = "grey80")
   at <- c(at, mean(xlim))
}
points(ix, b$carapace.width[ix], col = cols[b$instar[ix]], cex = 0.1)
for (i in 1:length(instars)) hline(exp(r$mu_pub_global[instars[i]]), col = cols[instars[i]], lwd = 2)
axis(1, labels = yy, at = at, las = 2)
axis(2)
axis(4, labels = as.roman(1:length(r$mu_pub_global)), at = exp(r$mu_pub_global), font = 2, padj = -0.5)
mtext("Carapace width (mm)", 2, 2.5, cex = 1.25)
mtext("Year", 1, 3.5, cex = 1.25)
mtext("Instar", 4, 2.0, cex = 1.25)
box(col = "grey60")
dev.off()

clg()
png(file = paste0("instar/Mature classification ", min(years), "-", max(years), ".png"), res = 500, units = "in", height = 8.5, width = 11)
par(mar = c(5, 4, 4, 4) + 0.1) # c(bottom, left, top, right)
plot(c(0, nrow(b)), c(40, 80), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", xaxs = "i", yaxs = "i")
ix <- which(b$maturity == "mature")
ix <- sample(ix)
instars <- sort(unique(b$instar[ix]))
#cols <- rainbow(length(instars))
cols <- brewer.pal(11, name = "Paired")
yy <- sort(unique(b$year))
at <- NULL
for (i in 1:length(yy)){
   xlim <- range(which(b$year == yy[i]))
   rect(xlim[1], par("usr")[3], xlim[2], par("usr")[4], col = c("grey80", "grey95")[((i%%2)==0)+1], border = "grey80")
   at <- c(at, mean(xlim))
}
points(ix, b$carapace.width[ix], col = cols[b$instar[ix]], cex = 0.1)
for (i in 1:length(instars)) hline(exp(r$mu_mat_global[instars[i]]), col = cols[instars[i]], lwd = 2)
axis(1, labels = yy, at = at, las = 2)
axis(2)
axis(4, labels = as.roman(1:length(r$mu_mat_global)), at = exp(r$mu_mat_global), font = 2, padj = -0.5)
mtext("Carapace width (mm)", 2, 2.5, cex = 1.25)
mtext("Year", 1, 3.5, cex = 1.25)
mtext("Instar", 4, 2.0, cex = 1.25)
box(col = "grey60")
dev.off()

# Generate instar abundance summary:
tab <- aggregate(list(n_imm = b$maturity == "immature"), by = b[c("date", "tow.id", "group")], sum)
tab$n_pub <- aggregate(list(n = b$maturity == "pubescent"), by = b[c("date", "tow.id", "group")], sum)$n
tab$n_mat <- aggregate(list(n = b$maturity == "mature"), by = b[c("date", "tow.id", "group")], sum)$n
n_imm <- as.data.frame(r$p_imm * repvec(tab$n_imm, ncol = ncol(r$p_imm)))
names(n_imm) <- paste0("n_instar_", 1:ncol(r$p_imm), "_imm")
n_pub <- as.data.frame(r$p_pub * repvec(tab$n_pub, ncol = ncol(r$p_pub)))
names(n_pub) <- paste0("n_instar_", 1:ncol(r$p_pub), "_pub")
n_mat <- as.data.frame(r$p_mat * repvec(tab$n_mat, ncol = ncol(r$p_mat)))
names(n_mat) <- paste0("n_instar_", 1:ncol(r$p_mat), "_mat")
tab <- cbind(tab, n_imm, n_pub, n_mat)
fvars <- names(tab)[grep("instar", names(tab))]
fvars <- fvars[which(apply(tab[fvars], 2, sum) > 0)]
tab <- tab[c("date", "tow.id", "group", "n_imm", "n_pub", "n_mat", fvars)]
ix <- match(tab[key(s)], s[key(s)])
tab$grid <- s$grid[ix]
tab$longitude <- lon(s)[ix]
tab$latitude  <- lat(s)[ix]
tab$swept.area <- s$swept.area[ix]
fvars <- names(tab)[grep("instar", names(tab))]
tab[fvars] <- 1000000 * tab[fvars] / repvec(tab$swept.area, ncol = length(fvars))
tab <- tab[!is.na(tab$grid), ]
tab$year <- year(tab)
tab <- aggregate(tab[fvars], by = tab[c("date", "year", "tow.id", "group", "grid")], mean)
   

plot(c(30,150), c(3, 11), type = "n")
for (i in 1990:2022) points(tab$depth[tab$year == i], log(tab$n_instar_8_imm[tab$year == i]), pch = 21, bg = rainbow(31)[i-1989])


exp(r$mu_imm_global + 0.5*r$sigma_imm[1,1]^2)
sqrt((exp(r$sigma_imm[1,1]^2)-1) * exp(2*r$mu_imm_global + r$sigma_imm[1,1]^2))

ii <- exp(r$mu_pub_global + 0.5*r$sigma_pub[1,1]^2)
names(ii) <- 1:length(ii)
jj <- sqrt((exp(r$sigma_pub[1,1]^2)-1) * exp(2*r$mu_pub_global + r$sigma_pub[1,1]^2))
names(jj) <- 1:length(jj)

ii <- exp(r$mu_mat_global + 0.5*r$sigma_mat[1,1]^2)
names(ii) <- 1:length(ii)
jj <- sqrt((exp(r$sigma_mat[1,1]^2)-1) * exp(2*r$mu_mat_global + r$sigma_mat[1,1]^2))
names(jj) <- 1:length(jj)

# 
delta_mu_imm <- exp(r$mu_imm) - repvec(exp(r$mu_imm_global), nrow = nrow(r$mu_imm))

# Generate instar size summary:
colnames(r$mu_imm) <- paste0("mu_instar_", 1:ncol(r$mu_imm), "_imm")
colnames(r$mu_pub) <- paste0("mu_instar_", 1:ncol(r$mu_pub), "_pub")
colnames(r$mu_mat) <- paste0("mu_instar_", 1:ncol(r$mu_mat), "_mat")
mu <- cbind(unique(b[c("date", "tow.id", "group")]), r$mu_imm, r$mu_pub, r$mu_mat) 
ix <- match(mu[key(s)], s[key(s)])
mu$grid <- s$grid[ix]
mu$longitude <- lon(s)[ix]
mu$latitude  <- lat(s)[ix]
mu$year <- year(mu)
fvars <- names(mu)[grep("instar", names(mu))]
ref <- aggregate(mu[fvars], mu["year"], function(x) as.numeric(names(table(x))[which.max(table(x))]))
rr <- aggregate(mu[fvars], mu["year"], mean)
rr[fvars] - ref[fvars]

gbarplot(exp(rr$mu_instar_7_imm) - exp(ref$mu_instar_7_imm), rr$year)
gbarplot(exp(rr$mu_instar_8_imm) - exp(ref$mu_instar_8_imm), rr$year)

gbarplot(exp(rr$mu_instar_8_pub) - exp(ref$mu_instar_8_pub), rr$year)
gbarplot(exp(rr$mu_instar_9_pub) - exp(ref$mu_instar_9_pub), rr$year)

gbarplot(exp(rr$mu_instar_10_mat) - exp(ref$mu_instar_10_mat), rr$year)

for (i in 2:length(years)){
   plot(c(2.5, 12.5), c(2.5, 12.5), type = "n", xaxs = "i", yaxs = "i")
   grid()
   ix <- which(tab$year == years[i-1])
   iy <- which(tab$year == years[i])
   
   xx <- tab[ix, ]
   yy <- tab[iy, ]
   
   iz <- match(xx$grid, yy$grid)
   
   # points(log(xx$n_instar_7_imm), log(yy$n_instar_8_pub[iz]) - log(xx$n_instar_7_imm), pch = 21, bg = rainbow(length(years))[i], cex = 0.75)
   #points(log(xx$n_instar_9_pub), log(yy$n_instar_10_mat[iz]), pch = 21, bg = rainbow(length(years))[i], cex = 1)
   
   #cex <- 0.01 * sqrt(xx$n_instar_7_imm + yy$n_instar_8_imm[iz] + yy$n_instar_8_pub[iz])
   #points(log(xx$n_instar_7_imm), log(yy$n_instar_8_imm[iz] + yy$n_instar_8_pub[iz]), pch = 21, bg = rainbow(length(years))[i], cex = cex)

   points(log(xx$n_instar_7_imm), log(yy$n_instar_8_imm[iz]), pch = 22, bg = rainbow(length(years))[i], cex = 2)
          
   fun <- function(x) mean(x[!is.na(x) & is.finite(x)])
   
   #points(fun(log(xx$n_instar_7_imm)), fun(log(yy$n_instar_8_imm[iz] + yy$n_instar_8_pub[iz])), pch = 22, bg = rainbow(length(years))[i], cex = 4)
}
box(col = "grey60")
legend("topleft", legend = years, pch = 22, pt.cex = 3, pt.bg = rainbow(length(years)))


clg()
for (i in 2:length(years)){
   xx <- tab[which(tab$year == years[i-1]), ]
   yy <- tab[which(tab$year == years[i]), ]
   
   iz <- match(xx$grid, yy$grid)
   
   map.new()
   grid()
   points(xx$longitude, xx$latitude, pch = 21, bg = adjustcolor("green", alpha.f = 0.3), cex = 0.02 * sqrt(xx$n_instar_7_imm))
   ix <- which(xx$n_instar_7_imm == 0)
   points(xx$longitude[ix], xx$latitude[ix], pch = "x", col = "green")
   points(yy$longitude, yy$latitude, pch = 21, bg = adjustcolor("red", alpha.f = 0.3), cex = 0.02 * sqrt(yy$n_instar_8_imm))
   points(yy$longitude, yy$latitude, pch = 21, bg = adjustcolor("blue", alpha.f = 0.3), cex = 0.02 * sqrt(yy$n_instar_8_pub))
   map("coast")
   title(main = paste0(years[i-1], "-", years[i]))
   box(col = "grey60")
   
   legend("topright", legend = years[c(i-1,i)], pch = 21, pt.cex = 3, 
          pt.bg = adjustcolor(c("green", "red"), alpha.f = 0.3), cex = 1.5)
   
   #points(fun(log(xx$n_instar_7_imm)), fun(log(yy$n_instar_8_imm[iz] + yy$n_instar_8_pub[iz])), pch = 22, bg = rainbow(length(years))[i], cex = 4)
}
box(col = "grey60")

plot(c(2.5, 12.5), c(0, 1), type = "n", xaxs = "i", yaxs = "i")
grid()
p <- p2 <- m <- NULL
for (i in 2:length(years)){
   xx <- tab[which(tab$year == years[i-1]), ]
   yy <- tab[which(tab$year == years[i]), ]

   iz <- match(xx$grid, yy$grid)
   
   points(log(xx$n_instar_7_imm), yy$n_instar_8_imm[iz] / (yy$n_instar_8_imm[iz] + yy$n_instar_8_pub[iz]), 
          pch = 21, bg = rainbow(length(years))[i], cex = 0.02 * sqrt(yy$n_instar_8_imm[iz] + yy$n_instar_8_pub[iz]))

   p <- c(p, mean(yy$n_instar_8_imm[iz], na.rm = TRUE) / (mean(yy$n_instar_8_imm[iz], na.rm = TRUE) + mean(yy$n_instar_8_pub[iz], na.rm = TRUE) ))
   p2 <- c(p2, mean(yy$n_instar_8_imm, na.rm = TRUE) / (mean(yy$n_instar_8_imm, na.rm = TRUE) + mean(yy$n_instar_8_pub, na.rm = TRUE) ))

   m <- c(m, (mean(yy$n_instar_8_imm[iz], na.rm = TRUE) + mean(yy$n_instar_8_pub[iz], na.rm = TRUE)) /  mean(xx$n_instar_7_imm, na.rm = TRUE))
}
box(col = "grey60")
legend("topleft", legend = years, pch = 22, pt.cex = 3, pt.bg = rainbow(length(years)))

gbarplot(p, years[-1], ylim = c(0.0, 1.0))
mtext("Proportion of instar VIIIs that are immature", 2, 2.5, cex = 1.25)
points(years[-1], p2)

# Mortality exploration:
gbarplot(m, years[-1])
mtext("Proportion of instar VIIIs that are immature", 2, 2.5, cex = 1.25)
points(years[-1], p2)

# Immature instars:
layout(rbind(0, cbind(0, kronecker(1:5, matrix(1, ncol = 8, nrow = 8)), 0), 0, 0, 0))
par(mar = c(0,0,0,0))
for (i in 4:8){
   gbarplot(aggregate(tab[, paste0("n_instar_", i, "_imm")], by = tab["year"], mean)[,2], years, grid = TRUE, xaxt = "n")
   box(col = "grey60")
   mtext(paste0("Instar ", as.roman(i)), 4, 1.0)
}
axis(1)

# Immature VII & Pubescent instars:
layout(rbind(0, cbind(0, kronecker(1:6, matrix(1, ncol = 8, nrow = 8)), 0), 0, 0, 0))
par(mar = c(0,0,0,0))
gbarplot(aggregate(tab[, paste0("n_instar_", 7, "_imm")], by = tab["year"], mean)[,2], years, grid = TRUE, xaxt = "n")
gbarplot(aggregate(tab[, paste0("n_instar_", 8, "_imm")], by = tab["year"], mean)[,2], years, grid = TRUE, xaxt = "n")
gbarplot(aggregate(tab[, paste0("n_instar_", 8, "_pub")], by = tab["year"], mean)[,2], years, grid = TRUE, xaxt = "n")
gbarplot(aggregate(tab[, paste0("n_instar_", 9, "_pub")], by = tab["year"], mean)[,2], years, grid = TRUE, xaxt = "n")
gbarplot(aggregate(tab[, paste0("n_instar_", 9, "_mat")], by = tab["year"], mean)[,2], years, grid = TRUE, xaxt = "n")
gbarplot(aggregate(tab[, paste0("n_instar_", 10, "_mat")], by = tab["year"], mean)[,2], years, grid = TRUE, xaxt = "n")
axis(1)

# Calculate mean sizes:
m <- NULL
for (i in 2:length(years)){
   xx <- tab[which(tab$year == years[i-1]), ]
   yy <- tab[which(tab$year == years[i]), ]
   
   iz <- match(xx$grid, yy$grid)
   
   points(log(xx$n_instar_7_imm), yy$n_instar_8_imm[iz] / (yy$n_instar_8_imm[iz] + yy$n_instar_8_pub[iz]), 
          pch = 21, bg = rainbow(length(years))[i], cex = 0.02 * sqrt(yy$n_instar_8_imm[iz] + yy$n_instar_8_pub[iz]))
   
   p <- c(p, mean(yy$n_instar_8_imm[iz], na.rm = TRUE) / (mean(yy$n_instar_8_imm[iz], na.rm = TRUE) + mean(yy$n_instar_8_pub[iz], na.rm = TRUE) ))
   p2 <- c(p2, mean(yy$n_instar_8_imm, na.rm = TRUE) / (mean(yy$n_instar_8_imm, na.rm = TRUE) + mean(yy$n_instar_8_pub, na.rm = TRUE) ))
   
   m <- c(m, (mean(yy$n_instar_8_imm[iz], na.rm = TRUE) + mean(yy$n_instar_8_pub[iz], na.rm = TRUE)) /  mean(xx$n_instar_7_imm, na.rm = TRUE))
}

clg()
for (i in 6000:7000){
   p     <- r$p_imm[i, ]
   mu    <- r$mu_imm[i, ]
   sigma <- r$sigma_imm[i, ] 
   t <- seq(0, 4.5, len = 1000)
   d <- rep(0, length(t))
   for (j in 1:length(mu)){
      d <- d + p[j] * dnorm(t, mu[j], sigma[j])
   }
   #plot(t, d, type = "l")
   ix <- which(data$group_imm == i-1)
   if (length(ix) > 10){
      rep(data$x_imm[ix], each = data$f_imm[ix])
      
      gbarplot(table(round(log(data$x_imm[ix])*25)/25), xlim = c(2, 4), col = "grey70")
      lines(t, 0.025 * sum(data$f_imm[ix])*d, col = "red", lwd = 2)
      
      points(log(data$x_imm[ix]), sum(data$f_imm[ix]) * r$d_imm[ix])
      
      vline(r$mu_imm_global, col = "blue", lty = "dashed", lwd = 2)
   }
}

clg()
for (i in 6600:7000){
   p     <- r$p_pub[i, ]
   mu    <- r$mu_pub[i, ]
   sigma <- r$sigma_pub[i, ] 
   t <- seq(0, 4.5, len = 1000)
   d <- rep(0, length(t))
   for (j in 1:length(mu)){
      d <- d + p[j] * dnorm(t, mu[j], sigma[j])
   }
   #plot(t, d, type = "l")
   ix <- which(data$group_pub == i-1)
   if (length(ix) > 10){
      rep(data$x_pub[ix], each = data$f_pub[ix])
      
      gbarplot(table(round(log(data$x_pub[ix])*50)/50), xlim = c(2, 4), col = "grey70")
      lines(t, 0.050 * sum(data$f_pub[ix])*d, col = "red", lwd = 2)
      
      points(log(data$x_pub[ix]), sum(data$f_pub[ix]) * r$d_pub[ix])
      
      vline(r$mu_pub_global, col = "blue", lty = "dashed", lwd = 2)
   }
}


res <- aggregate(list(n_imm = data$f_imm), by = data["group_imm"], sum)
r$n_imm <- r$p_imm[res$group_imm, ] * repvec(res$n_imm, ncol = ncol(r$p_imm))

plot(c(0,max(data$group_imm)), c(10, 80))
points(exp(r$mu_imm[,4]), col = "blue", cex = 0.5)
points(exp(r$mu_imm[,5]), col = "blue", cex = 0.5)
points(exp(r$mu_imm[,6]), col = "blue", cex = 0.5)
points(exp(r$mu_imm[,7]), col = "blue", cex = 0.5)
points(exp(r$mu_imm[,8]), col = "blue", cex = 0.5)
points(exp(r$mu_imm[,9]), col = "blue", cex = 0.5)
points(exp(r$mu_pub[,8]), col = "purple", cex = 0.5)
points(exp(r$mu_pub[,9]), col = "purple", cex = 0.5)
points(exp(r$mu_pub[,10]), col = "purple", cex = 0.5)
points(exp(r$mu_mat[,9]), col = "palegreen3", cex = 0.5)
points(exp(r$mu_mat[,10]), col = "palegreen3", cex = 0.5)
points(exp(r$mu_mat[,11]), col = "palegreen3", cex = 0.5)
hline(exp(r$mu_imm_global))




plot(sample(data$x_imm), cex = 0.6)
hline(exp(r$mu_imm_global), lwd = 2, col = "red")
points(exp(r$mu_imm[,5]), col = "blue")
points(exp(r$mu_imm[,6]), col = "red")
points(exp(r$mu_imm[,7]), col = "blue")
points(exp(r$mu_imm[,8]), col = "green")


plot(sample(data$x_pub), cex = 0.6)
hline(exp(r$mu_pub_global), lwd = 2, col = "red")
points(exp(r$mu_pub[,5]), col = "blue")
points(exp(r$mu_pub[,6]), col = "red")
points(exp(r$mu_imm[,7]), col = "blue")
points(exp(r$mu_imm[,8]), col = "green")

plot(sample(data$x_mat), cex = 0.6)
hline(exp(r$mu_mat_global), lwd = 2, col = "red")


# Plot instar classification results:
plot(r$mu_imm[, 8])
plot(r$mu_pub[, 8])

p <- apply(r$p_pub, 2, mean)
sigma <- apply(r$sigma_pub, 2, mean)
mu <- apply(r$mu_pub, 2, mean)

x <- seq(0, 5, len = 1000)
d <- rep(0, length(x))
for (i in 1:length(mu)) d <- d + p[i] * dnorm(x, mu[i], sigma[i])
gbarplot(table(round(log(b$carapace.width[b$maturity == "pubescent"]), 2)))
lines(x, 100*d)

