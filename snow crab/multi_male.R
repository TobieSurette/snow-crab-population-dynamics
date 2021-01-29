library(TMB)
library(gulf.data)
library(gulf.spatial)
library(gulf.graphics)
library(gulf.stats)

# Model parameters:
years  <- 2006:2020
sex <- 1
step <- 0.5

# Set derived quantities:
n_year <- length(years)
if (sex == 1){
   n_instar <- 10
   xlim <- c(0, 140)
   ylim <- c(0, 200)
}else{
   n_instar <- 10
   xlim <- c(0, 100)
   ylim <- c(0, 400)
}
instars <- as.character(as.roman(4:(4+n_instar-1)))

setwd("snow crab")
#source("snow crab/male/instar.year.data.R")
if (sex == 1) load("males.2006-2020.rdata")
if (sex == 2) load("females.2006-2020.rdata")

# Work computer fix:
if (Sys.getenv("RSTUDIO_USER_IDENTITY") == "SuretteTJ") Sys.setenv(BINPREF = "C:/Rtools/mingw_64/bin/")

compile("multi_male.cpp")
dyn.load(dynlib("multi_male"))

# Define initial parameters:
parameters <- list(mu0                 = 10,                             # First instar mean size.
                   log_sigma0          = log(0.8),                       # Log-scale standard error for first instar.
                   log_hiatt_slope     = log(c(0.350, 0.055)),           # Hiatt slope parameters.
                   log_hiatt_intercept = log(c(0.689, 10.000)),          # Hiatt intercept parameters.
                   log_growth_error    = log(c(0.01, 0.25)),             # Growth increment error inflation parameters.
                   mu_year_instar      = rep(0, n_instar * n_year),      # Log-scale instar mean year interaction (n_instar x n_year).
                   log_sigma_mu_year_instar   = -0.7,                    # Instar mean year interaction error term.
                   delta_mat = -1.5,                                     # Size offset between mature and immature instar sizes.
                   log_n_imm_year_0    = rep(4, n_instar-1),             # First year immature instar abundances (n_instar-1).
                   log_n_imm_instar_0  = rep(4, n_year),                 # First instar recruitment for all years (n_year).
                   log_sigma_n_imm_instar_0 = -1,                        # Log-scale first instar annual recruitment error parameter.
                   log_n_skp_instar_0  = rep(0, n_instar-5),             # First year skip abundances (n_instar-5).
                   log_n_mat_instar_0 = rep(0, (n_instar-5)*6),          # First year mature group abundances (n_instar-5)x6.
                   selectivity_x50     = 25,                             # Size-at-50% trawl selectivity.
                   log_selectivity_slope = -1,                           # Log-scale trawl selectivity slope.
                   log_year_effect = rep(0, n_year),                     # Abundance year effect (n_year).
                   log_sigma_year_effect = -2,                           # Log-scale year effect error parameter.
                   logit_p_skp = c(rep(-8, 4), rep(-6, n_instar-5), -8), # Logit-scale skip-moulting probabilities (n_instar).
                   logit_p_mat = c(rep(-8, 4), rep(0, n_instar-6)),      # Logit-scale moult-to-maturity probabilities (n_instar-2).
                   logit_p_mat_year = rep(0, (n_instar-2) * (n_year-1)), # Logit-scale moult-to-maturity instar x year interaction (n_instar-2 x n_year-1).
                   log_sigma_p_mat_year = -1,                            # Moult-to-maturity instar x year interaction error term.
                   logit_M_imm = -2,                                      # Logit-scale immature mortality.
                   logit_M_mat = c(-1.10, -1.0),                          # Logit-scale mature mortality.
                   selectivity_x50_fishing = 95,                         # Size-at-50% selectivity for fishing.
                   log_selectivity_slope_fishing = -1.5,                 # Log-scale trawl selectivity slope for fishing.
                   logit_fishing_effect_skp = -1,
                   logit_fishing_effect_rec = -1,
                   logit_fishing_effect_res = 0)

#parameters$log_hiatt_slope <- log(c(0.38619475, 0.05))

#parameters$log_hiatt_intercept <- log(c(1.28, 5))

#parameters$log_growth_error <- c(-4.390540, -1.6)
#parameters$delta_mat <- -1.5

#parameters$log_year_effect <- rep(0, n_year)

parameters$logit_M_imm = -2                                     # Logit-scale immature mortality.
parameters$logit_M_mat = c(-1.10, -1.0)                         # Logit-scale mature mortality.
parameters$selectivity_x50_fishing = 95                         # Size-at-50% selectivity for fishing.
parameters$log_selectivity_slope_fishing = 1.5                 # Log-scale trawl selectivity slope for fishing.
parameters$logit_fishing_effect_skp = -1
parameters$logit_fishing_effect_rec = 2
parameters$logit_fishing_effect_res = 2


# Define random variables in model:
random <- c("mu_year_instar", "log_n_imm_instar_0", "logit_p_mat_year", "log_year_effect")
data.vars <- names(data)[-grep("(rec)|(res)|(skp)", names(data))]

# Initialize parameter mapping:
map <- lapply(parameters, function(x) factor(rep(NA, length(x))))

# Estimate initial abundance parameters:
map <- update.map(map, free = c("log_n_imm_instar_0", "log_n_imm_year_0", "log_n_skp_instar_0", "log_n_mat_instar_0", "log_sigma_n_imm_instar_0"))
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
parameters <- update.parameters(parameters, obj, map = map)

# Add moult to maturity parameters:
map <- update.map(map, free = c("logit_p_mat", "logit_p_mat_year", "log_sigma_p_mat_year"))
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 5000))$par
parameters <- update.parameters(parameters, obj, map = map)

# Add year effect parameters:
map <- update.map(map, free = c("log_sigma_year_effect"))
parameters$log_year_effect[length(parameters$log_year_effect)] <- 0
map$log_year_effect <- factor(c(1:(length(parameters$log_year_effect)-1 ), NA))
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 500))$par
parameters <- update.parameters(parameters, obj, map = map)

# Add fishing selectivity parameters:
map <- update.map(map, free = c("selectivity_x50_fishing", "log_selectivity_slope_fishing", "logit_fishing_effect_skp", "logit_fishing_effect_rec", "logit_fishing_effect_res"))
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1500))$par
parameters <- update.parameters(parameters, obj, map = map)

# Add skip-moulting probability parameters:
map$logit_p_skp <- as.factor(c(rep(NA, 5), 1:(length(parameters$logit_p_skp)-5)))
map <- update.map(map, free = c("logit_p_skp"))
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
parameters <- update.parameters(parameters, obj, map = map)

parameters$log_selectivity_slope_fishing <- -1.5
parameters$logit_fishing_effect_skp <- 1
parameters$logit_fishing_effect_rec <- 1
parameters$logit_fishing_effect_res <- 0

# Add mortality parameters:
map <- update.map(map, free = c("logit_M_mat","logit_M_imm"))
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1500))$par
parameters <- update.parameters(parameters, obj, map = map)

# Add selectivity parameters:
map <- update.map(map, free = c("selectivity_x50", "log_selectivity_slope"))
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1500))$par
parameters <- update.parameters(parameters, obj, map = map)

parameters$log_selectivity_slope_fishing <- 1.5

# Add delta_mat:
map <- update.map(map, free = c("delta_mat"))
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 300))$par
parameters <- update.parameters(parameters, obj, map = map)

# Add annual growth parameters:
map <- update.map(map, free = c("mu_year_instar"))
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 300))$par
parameters <- update.parameters(parameters, obj, map = map)

# Add some growth parameters:
map <- update.map(map, free = c("log_hiatt_slope", "log_hiatt_intercept", "log_sigma0"))
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 500))$par
parameters <- update.parameters(parameters, obj, map = map)

# Add instar error parameter:
map <- update.map(map, free = c("log_growth_error"))
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 200))$par
parameters <- update.parameters(parameters, obj, map = map)

# Add annual growth parameters:
map <- update.map(map, free = "log_sigma_mu_year_instar")
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 1000))$par
parameters <- update.parameters(parameters, obj, map = map)

# Add annual growth parameters:
map <- update.map(map, free = "mu0") 
obj <- MakeADFun(data[data.vars], parameters, DLL = "multi_male",  random = random, map = map)
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 2000))$par
parameters <- update.parameters(parameters, obj, map = map)


parameters$mu_year_instar[abs(parameters$mu_year_instar) > 5] <- 0
parameters$log_sigma_mu_year_instar[1] <- 1

