library(gulf.data)
library(gulf.graphics)
library(TMB)
source("U:/TMB/TMB utilities.R")

setwd("C:/Users/SuretteTJ/Desktop/github/snow-crab-stock-assessment-2022/R/population model")
clc(); compile("model_basis.cpp")
dyn.load(dynlib("model_basis"))

# Prepare data:
load("C:/Users/SuretteTJ/Desktop/github/snow-crab-stock-assessment-2022/results/tables/size-frequencies males scs 2001-2022.rdata")
f <- f[, as.character(1:150), ]
f[,as.character(1:10),"immature"] <- 0
f[,as.character(121:150),"immature"] <- 0
f[,as.character(1:33),"skip"] <- 0
f[,as.character(121:150),"skip"] <- 0
f[,as.character(1:40),"recruit"] <- 0
f[,as.character(1:40),"residual"] <- 0
data <- list(x = as.numeric(dimnames(f)[[2]]),
             f_imm = f[,,"immature"] + f[,,"skip"],
             f_mat = f[,,"recruit"] + f[,,"residual"])

# Define instar basis functions:
mu_instars <- c(14.5, 20.5, 28.1, 37.8, 51.1, 68.0, 88.0, 106, 127)
names(mu_instars) <- 5:13
sigma_instars <- 0.08 * mu_instars
basis <- matrix(NA, nrow = length(mu_instars), ncol = length(data$x))
for (i in 1:length(mu_instars)){
   basis[i, ] <- sqrt(2*pi) * sigma_instars[i] * (pnorm(data$x + 0.5, mu_instars[i], sigma_instars[i]) - pnorm(data$x - 0.5, mu_instars[i], sigma_instars[i]))  
}
plot(c(0, 140), c(0, 1.1), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")
grid()
for (i in 1:nrow(basis)){
   lines(data$x, basis[i,], lty = "dashed", col = rainbow(length(mu_instars))[i], lwd = 2)
   text(mu_instars[i], 1.035, as.roman(as.numeric(names(mu_instars[i]))))
}
#lines(data$x, apply(basis, 2, sum), lwd = 2, col = "grey30")
mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)      
mtext("Value", 2, 2.5, cex = 1.25)  
years <- as.numeric(rownames(f))
box(col = "grey60")
data$basis <- basis

logit <- function(x) return(log(x/(1-x)))

# Define initial parameters: 
parameters <- list(intercept_growth = 0.276,             # Growth intercept parameter.
                   xp_growth = 38.2,                     # Growth pivot point parameter.
                   log_window_growth = 1.6,              # Growth transition window parameter.
                   slope_growth = c(0.32, 0.100),        # Growth slope parameters.
                   log_sigma_growth = log(0.135),        # Growth error parameter.
                   log_shape_selectivity = log(5),       # Shape parameter for cumulative gamma distribution function.
                   log_scale_selectivity = log(10),      # Scale parameter for cumulative gamma distribution function.
                   xp_maturity = c(55, 105),             # Maturation pivot point parameters.
                   log_window_maturity = c(-1.5, -0.1),  # Maturation transition window parameters.
                   logit_p_mix_maturity = 0.1,           # Mixing proportion between early and late maturation (logit-scale).
                   xp_skp = 70,                          # Skip-moulting pivot point parameter.
                   log_window_skp = -1.6,                # Skip-moulting transition window parameter.
                   logit_p_max_skp = 0,                  # Maximum skip-moulting probability (logit-scale).
                   xp_F = 98,                                       # Exploitation rate switching point.
                   log_window_F = -1,                               # Exploitation switching point window.
                   log_mu_F = logit(0.6),                           # Global exploitation rate.      
                   log_dev_F = 0.0 * rnorm(dim(data$f_imm)[1]),     # Annual deviations exploitation rate.
                   log_mu_M_imm = logit(0.25),                      # Global natural mortality rate.
                   log_dev_M_imm = rep(0, dim(data$f_imm)[1]),      # Annual deviations natural mortality rate.
                   log_mu_M_mat = logit(0.40),                      # Global natural mortality rate.
                   log_dev_M_mat = rep(0, dim(data$f_imm)[1]),      # Annual deviations natural mortality rate.
                   size_recruitment_mu = 14,                        # Mean recruitment size.
                   log_sigma_size_recruitment = -2,                 # Error for recruitment annual deviations.
                   size_recruitment = rep(0,length(years)),         # Annual deviations from recruitment instar mean size.
                   log_sigma_recruitment_instar = log(0.8),         # Error for recruitment instar.
                   log_scale_recruitment_mu = log(30000),           # Recruitment scale mean.
                   log_sigma_scale_recruitment = 10,                # Recruitment scale error.
                   log_scale_recruitment = rep(0,length(years)),    # Recruitment scale deviations.
                   log_N_imm_reg_beta = rep(1, nrow(data$basis)),   # Basis coefficients for immature crab for year 0.    
                   log_N_imm_skp_beta = rep(1, nrow(data$basis)-4), # Basis coefficients for immature skip-moulting crab for year 0.
                   log_N_mat_new_beta = rep(1, nrow(data$basis)-4), # Basis coefficients for new mature crab for year 0.
                   log_N_mat_old_beta = rep(1, nrow(data$basis)-4)) # Basis coefficients for old mature crab for year 0. 

# Set parameters:
parameters$log_sigma_dev_F = 0                            # Error for annual exploitation rate random effect.
parameters$log_sigma_dev_M_imm = 0                        # Error for annual mortality rate random effect.
parameters$log_sigma_dev_M_mat = 0                        # Error for annual mortality rate random effect.

# Error parameters:
parameters$log_r_imm <- 0.1
parameters$log_r_mat <- 0.1

# Adjust to better initial parameters values:
parameters$log_mu_F <- 0

parameters$log_scale_recruitment <- c(-0.8,-0.7,-0.2,0,0.2,-0.2,-0.3,0.2,0.1,-0.2,-0.2,-0.5,0.5,-0.2,0.3,0.2,-1.1,1.5,0.9,0.4,-1.2,-1.2)
parameters$log_scale_recruitment_mu <- 10.4
parameters$log_sigma_scale_recruitment <- -0.17

tolist <- function(x){
   str <- unique(names(x))
   r <- list()
   for (i in 1:length(str)){
      r[[i]] <- as.numeric(x[str[i] == names(x)])
   }
   names(r) <- str
   return(r)
}
update <- function(x, theta){
   theta <- tolist(theta)
   for (i in 1:length(theta)){
      x[[names(theta[i])]] <- theta[[i]]
   }
   return(x)
}

# parameters$size_recruitment <- 0.5*rnorm(length(years))

obj <- MakeADFun(data = data, parameters = parameters, DLL = "model_basis")
r <- obj$report()
image(r$p_recruitment)

map <- lapply(parameters, function(x) as.factor(NA * x))

# Fit dispersion parameters:
map$log_r_imm <- as.factor(1)
map$log_r_mat <- as.factor(1)
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
parameters <- update(parameters, theta)

# Add initial immature:
map$log_N_imm_reg_beta <- as.factor(1:length(map$log_N_imm_reg_beta))
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 2500))$par
parameters <- update(parameters, theta)

# Add initial skip:
map$log_N_imm_skp_beta <- as.factor(1:length(map$log_N_imm_skp_beta))
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 2500))$par
parameters <- update(parameters, theta)

# Add initial new mature:
map$log_N_mat_new_beta <- as.factor(1:length(map$log_N_mat_new_beta))
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 2500))$par
parameters <- update(parameters, theta)

# Add initial old mature:
map$log_N_mat_old_beta <- as.factor(1:length(map$log_N_mat_old_beta))
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 3000))$par
parameters <- update(parameters, theta)

# Add recruitment:
map$size_recruitment_mu <- as.factor(1)                                     # Mean recruitment size.
map$log_scale_recruitment <- as.factor(1:length(map$log_scale_recruitment)) # Recruitment scale deviations.
map$log_scale_recruitment_mu <- as.factor(1)                                # Recruitment scale mean.
map$log_sigma_scale_recruitment <- as.factor(1)                             # Recruitment scale error.
#map$log_sigma_recruitment <- as.factor(1)
#size_recruitment_mu = 14,                        
#log_sigma_size_recruitment = -2,                 # Error for recruitment annual deviations.
#size_recruitment = rep(0,length(years)),         # Annual deviations from recruitment instar mean size.
#log_sigma_recruitment_instar = log(0.8),         # Error for recruitment instar.
#log_scale_recruitment_mu = log(30000),           
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 2500))$par
parameters <- update(parameters, theta)

# Add maturity parameters:
#map$log_window_maturity <- as.factor(1:length(map$log_window_maturity))
map$logit_p_mix_maturity <- as.factor(1)
map$xp_maturity <- as.factor(1:length(map$xp_maturity))
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 5000))$par
parameters <- update(parameters, theta)

# Add selectivity parameters:
map$log_shape_selectivity <- as.factor(1)
map$log_scale_selectivity <- as.factor(1)
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 3000))$par
parameters <- update(parameters, theta)

# Add skip-moulters parameters:
map$xp_skp <- as.factor(1)
map$logit_p_max_skp <- as.factor(1)
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 5000))$par
parameters <- update(parameters, theta)

# Add fishing mortality:
map$log_mu_F <- as.factor(1)
map$log_dev_F <- as.factor(1:length(map$log_dev_F))
map$log_sigma_dev_F <- as.factor(1)

map$log_window_F <- as.factor(1)
map$xp_F <- as.factor(1)
for (i in 1:5){
   obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
   theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 5000))$par
   parameters <- update(parameters, theta)
}

# Add mortality parameters:
map$log_mu_M_imm <- as.factor(1)
map$log_dev_M_imm <- as.factor(1:length(map$log_dev_M_imm))
map$log_mu_M_mat <- as.factor(1)
map$log_dev_M_mat <- as.factor(1:length(map$log_dev_M_mat))
map$log_sigma_dev_M_imm <- as.factor(1)
map$log_sigma_dev_M_mat <- as.factor(1)
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 5000))$par
parameters <- update(parameters, theta)

# Add growth parameters:
map$intercept_growth <- as.factor(1)
map$xp_growth
map$log_window_growth
map$slope_growth <- as.factor(1:length(map$slope_growth))
map$log_sigma_growth <- as.factor(1)
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 5000))$par
parameters <- update(parameters, theta)

# Add recruitment sizes:
#map$log_sigma_recruitment <- as.factor(1)
#size_recruitment_mu = 14,                        
map$log_sigma_size_recruitment <- as.factor(1)                    # Error for recruitment annual deviations.
map$size_recruitment <- as.factor(1:length(map$size_recruitment)) # Annual deviations from recruitment instar mean size.
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 2500))$par
parameters <- update(parameters, theta)

r <- obj$report()
dimnames(r$n_imm) <- list(year = (min(years)-1):max(years), cw = colnames(f))
dimnames(r$n_mat) <- list(year = (min(years)-1):max(years), cw = colnames(f))

random <- summary(obj)

gbarplot(exp(parameters$log_scale_recruitment_mu + parameters$log_scale_recruitment))
gbarplot(exp(parameters$log_N_imm_reg_beta))
gbarplot(exp(parameters$log_N_imm_skp_beta))
gbarplot(exp(parameters$log_N_mat_new_beta))
gbarplot(exp(parameters$log_N_mat_old_beta))
gbarplot(r$N_imm_reg[1,])
gbarplot(r$N_imm_skp[1,])
gbarplot(r$N_mat_new[1,])
gbarplot(r$N_mat_old[1,])

gbarplot(r$N_imm_reg[1,] + r$N_imm_skp[1,], data$x)
lines(data$x, r$n_imm[1,])
gbarplot(r$N_mat_new[1,] + r$N_mat_old[1,], data$x)
lines(data$x, r$n_mat[1,])

plot(r$p_selectivity, type = "l", lwd = 2)
plot(r$p_maturity, type = "l", lwd = 2)
plot(r$p_skp, type = "l", lwd = 2)
gbarplot(1/(1+ exp(-exp(-parameters$log_mu_F - parameters$log_dev_F))))
plot(r$mu_growth, type = "l", lwd = 2)
gbarplot(1/(1+ exp(-parameters$log_mu_M_imm - parameters$log_dev_M_imm)))
gbarplot(1/(1+ exp(-parameters$log_mu_M_mat - parameters$log_dev_M_mat)))

#F <- exp(parameters$log_mu_F + parameters$log_dev_F[i]) / (1+exp(-exp(parameters$log_window_F) * (data$x - parameters$xp_F)))
gbarplot(r$F[1,], data$x)

clg()
#windows()
r <- obj$report()
dimnames(r$n_imm) <- list(year = (min(years)-1):max(years), cw = colnames(f))
dimnames(r$n_mat) <- list(year = (min(years)-1):max(years), cw = colnames(f))

m <- kronecker(matrix(1:21, ncol = 3), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0, cbind(0, 0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))

for (i in 2:length(years)){
   year <- as.character(years[i])
   gbarplot(cbind(r$n_imm[year, ], r$n_mat[year, ]) , as.numeric(colnames(f)), 
            legend = FALSE, grid = TRUE, xaxs = "i", xlim = c(0, 140), yaxs = "i", # ylim = c(0, 800),  
            xaxt = "n", yaxt = "n")
   lines(data$x, data$f_imm[year, ] + data$f_mat[year, ], lwd = 2, col = "blue")
   lines(data$x, data$f_imm[year, ], lwd = 2, col = "red")
   
   if (i %in% 2:8) axis(2)
   if (i %in% (c(7, 14, 21)+1)) axis(1)
   text(par("usr")[1] + 0.8 * diff(par("usr")[1:2]),
        par("usr")[3] + 0.8 * diff(par("usr")[3:4]), year)
   
   vline(c(9.4, 14, 20.3, 28.5, 38), lty = "dashed")
   box(col = "grey60")
}

random <- c("log_dev_F", "log_dev_M_imm", "log_dev_M_mat", "log_scale_recruitment", "size_recruitment")

obj <- MakeADFun(data = data, parameters = parameters, random = random, map = map, DLL = "model_basis")
theta <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
parameters <- update(parameters, theta)

rep <- sdreport(obj)
theta <- c(theta, summary(rep, "random")[,1])
parameters <- update(parameters, theta)
# parameters <- update.parameters(parameters, summary(rep, "fixed"), summary(rep, "random"))



