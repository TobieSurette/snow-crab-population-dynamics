library(gulf.data)
library(TMB)

# Compile TMB model:
clc()
compile("instar/instar.cpp")
dyn.load(dynlib("instar/instar"))

# Prepare data:
years <- 2015:2022
s <- read.scsset(years, valid = 1, survey = "regular")
b <- read.scsbio(years)
b <- b[!is.na(match(b[key(s)], s[key(s)])), ]
s$group <- match(s[key(s)], unique(s[key(s)]))-1
b$group <- s$group[match(b[key(s)], s[key(s)])]
b$year <- year(b)

# recommendation: Group by grid for testing.

# Define maturity stages:
b$maturity <- ""
b$maturity[which(!is.na(b$carapace.width) & (b$sex == 2) & (b$carapace.width <= 70) & !is.mature(b) & !is.pubescent.scsbio(b))] <- "immature"
b$maturity[which(!is.na(b$carapace.width) & (b$carapace.width >= 28) & (b$carapace.width <= 80) & (b$sex == 2) & !is.mature(b) & is.pubescent.scsbio(b))] <- "pubescent"
b$maturity[which(!is.na(b$carapace.width) & (b$carapace.width >= 36) & (b$carapace.width <= 95) & (b$sex == 2) & is.primiparous.scsbio(b))] <- "mature"
b <- b[b$maturity != "", ]

tmp <- aggregate(b$carapace.width, by = b["group"], length)

b <- b[b$group %in% tmp[rev(order(tmp$x))[1:100], "group"], ]
b$group <- match(b$group, sort(unique(b$group)))

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
parameters <- list(mu_imm_0 = 0.7930157,             # Size of immature first instar.
                   a_imm = 0.890443829,              # Slope parameter for immature growth.
                   b_imm = 0.645895530,              # Intercept parameter for immature growth.
                   a_pub = 0.663637366,              # Slope parameter for pubescent growth.
                   b_pub = 1.515007633,              # Intercept parameter for pubescent growth.
                   a_mat = 1.040998406,              # Slope parameter for mature growth.
                   b_mat = 0.004727189,              # Intercept parameter for mature growth.   
                   log_sigma_delta_mu_imm = 0,
                   log_sigma_delta_mu_pub = 0,
                   log_sigma_delta_mu_mat = 0,
                   delta_mu_imm = rep(0, data$n_instar * data$n_group),
                   delta_mu_pub = rep(0, data$n_instar * data$n_group),
                   delta_mu_mat = rep(0, data$n_instar * data$n_group),
                   log_sigma = -2.3685149,
                   logit_p_imm_global = rep(0,5),
                   logit_p_pub_global = rep(0,2),
                   logit_p_mat_global = rep(0,2),
                   delta_logit_p_imm = matrix(2*runif(data$n_group * 5)-1, nrow = data$n_group, ncol = 5),
                   delta_logit_p_pub = matrix(2*runif(data$n_group * 2)-1, nrow = data$n_group, ncol = 2),
                   delta_logit_p_mat = matrix(2*runif(data$n_group * 2)-1, nrow = data$n_group, ncol = 2),
                   log_sigma_delta_logit_p_imm = -1,
                   log_sigma_delta_logit_p_pub = -1,
                   log_sigma_delta_logit_p_mat = -1) 

# Define functions:
free <- function(x, p){
   if (missing(p)){
      d <- dim(x) 
      v <- as.factor(1:length(x))
      dim(v) <- d
      return(v)
   }
   p <- p[p %in% names(x)]
   for (i in 1:length(p)) x[[p[i]]] <- free(x[[p[i]]])
   return(x)
}

update <- function(x, p, r){
   if (!missing(p)){
       str <- unique(names(p))
       for (i in 1:length(str)){
          d <- dim(x[[str[i]]])
          x[[str[i]]] <- p[names(p) == str[i]]
          dim(x[[str[i]]]) <- d
       }
   }
   if (!missing(r)){
       str <- unique(names(r))
       for (i in 1:length(str)){
          d <- dim(x[[str[i]]])
          x[[str[i]]] <- r[names(r) == str[i]]
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
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 2000))$par
parameters <- update(parameters, obj$par)

# Fit group proportions:
map <- free(map, c("delta_logit_p_imm", "log_sigma_delta_logit_p_imm", 
                   "delta_logit_p_pub", "log_sigma_delta_logit_p_pub",
                   "delta_logit_p_mat", "log_sigma_delta_logit_p_mat"))
obj <- MakeADFun(data = data, parameters = parameters, map = map, DLL = "instar")
obj$par <- optim(obj$par, obj$fn, control = list(trace = 3, maxit = 10000))$par
parameters <- update(parameters, obj$par)

rep <- sdreport(obj)
summary(rep, "random") 
parameters <- update(parameters, obj$par)

# Fit second:
log_sigma 

