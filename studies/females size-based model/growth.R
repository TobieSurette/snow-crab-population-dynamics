# Snow crab - Framework Assessment 2021

# Objectives:
# Refine survey catch standardization:
# - Change swept area assessment method to Hierarchical Gaussian Process model.
# - Retroactively re-estimate touchdown times using Minilog + Acoutic Depth profiles.
# - Project trawl track using vessel winch profiles.
# - Retroactively estimate passive phase swept area with error.
# - Summarize winch and vessel differences through time.
# - Show depth-temperature relationships through time.
# - Incorporate season time as a predictor (through temperature?)
#
# Time series standardization:
# - Global population model.
# - Prediction using population model.
# - Spatial investigations.

clc(); compile("growth.cpp")
dyn.load(dynlib("growth"))

b <- read.scsbio(1998:2021, survey = "regular", sex = 1)
b <- b[!is.mature(b), ]
b <- b[which((b$carapace.width >= exp(2)) & (b$carapace.width < exp(5))), ]

t <- table(round(log(b$carapace.width),2))
f <- seq(2, 5, by = 0.01)
names(f) <- f
f <- f * 0
f[names(t)] <- as.numeric(t)

# Define data:
data <- list(dx = NA, 
             x = round(as.numeric(names(f)), 2),
             f = as.numeric(f))
data$dx <- unique(round(diff(data$x),2))

# Specify initial parameter values:
parameters = list(intercept_growth = 0.276,
                  xp_growth = 38.2,
                  log_window_growth = 1.6,
                  slope_growth = c(0.32, 0.126),
                  log_sigma_growth = log(0.135))


map <- lapply(parameters, function(x) factor(rep(NA, length(x))))

# Estimate initial abundance parameters:
obj <- MakeADFun(data = data, parameters = parameters, DLL = "growth")

gbarplot(data$f, data$x)
lines(data$x, obj$report()$mu, col = "red", lwd = 2)

obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1500))$par
parameters <- update.parameters(parameters, obj, map = map)


plot(data$x, obj$report()$mu_growth, lwd = 2, type = "l")
lines(data$x, obj$report()$mu_growth - obj$report()$sigma_growth, lwd = 2, lty= "dashed", col =  "red")
lines(data$x, obj$report()$mu_growth + obj$report()$sigma_growth, lwd = 2, lty= "dashed", col =  "red")

