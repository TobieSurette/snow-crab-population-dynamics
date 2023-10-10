library(gulf.data)
library(gulf.graphics)
library(TMB)
source("U:/TMB/TMB utilities.R")

setwd("C:/Users/SuretteTJ/Desktop/github/snow-crab-stock-assessment-2022/R/population model")
clc(); compile("basis.cpp")
dyn.load(dynlib("basis"))

# Prepare data:
load("C:/Users/SuretteTJ/Desktop/github/snow-crab-stock-assessment-2022/results/tables/size-frequencies males scs 2001-2022.rdata")
f <- f[, as.character(1:150), ]
f[,as.character(1:10),"immature"] <- 0
f[,as.character(121:150),"immature"] <- 0
f[,as.character(1:33),"skip"] <- 0
f[,as.character(121:150),"skip"] <- 0
f[,as.character(1:40),"recruit"] <- 0
f[,as.character(1:40),"residual"] <- 0
data <- list(x = as.numeric(dimnames(f)[[2]]))

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

# Add basis to data:
data$basis <- basis

# Initialize parameters:
parameters <- list(log_N_imm_reg_beta = rnorm(nrow(data$basis)),       
                   log_N_imm_skp_beta = rnorm(nrow(data$basis)),
                   log_N_mat_new_beta = rnorm(nrow(data$basis)),
                   log_N_mat_old_beta = rnorm(nrow(data$basis)))


# Create TMB model object:
obj <- MakeADFun(data = data, parameters = parameters, DLL = "basis")

print(parameters$log_N_imm_reg_beta)
plot(obj$report()$N_mat_old_0, type = "l")

# Plot skip-moulting function:
plot(data$x, obj$report()$p_skip, type = "l", lwd = 2, col = "blue")
gbarplot(data$f, data$x, grid = TRUE)
lines(data$x, obj$report()$f_skip, lwd = 2, col = "red")

# Plot maturation function:
plot(data$x, obj$report()$p_maturity, type = "l", lwd = 2, col = "blue")
gbarplot(data$f, data$x, grid = TRUE)
lines(data$x, obj$report()$f_mat, lwd = 2, col = "red")

