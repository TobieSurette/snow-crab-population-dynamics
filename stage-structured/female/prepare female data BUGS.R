library(gulf.data)
library(gulf.stats)
library(gulf.graphics)
library(gulf.spatial)

# Survey year:
years <- 1997:2025

# Load kriging polygons:
p <- read.gulf.spatial("kriging polygons revised")["gulf"]

# Read three years of data (for variogram averaging):
s <- read.scsset(year = (min(years)-2):max(years), survey = "regular", valid = 1) # Tow data.
b <- read.scsbio(year = (min(years)-2):max(years), survey = "regular", species = 2526, sex = 2)  # Biological data.]

# Define female categories:
b$category <- ""
b$category[which(!is.mature(b) & (b$gonad.colour %in% 1:2))]   <- "immature"
b$category[which(!is.mature(b) & (b$gonad.colour %in% 3:4))]   <- "adolescent"
b$category[which(is.mature(b) & (b$shell.condition %in% 1:2))] <- "primiparous"
b$category[which(is.mature(b) & (b$shell.condition %in% 3:5))] <- "multiparous"
b <- b[which(b$category != ""), ]

# Import catch data:
categories <- c("immature", "adolescent", "primiparous", "multiparous")
for (i in 1:length(categories)){
   import(s, fill = 0) <- catch(b[b$category == categories[i], ]) # Merge catches.
   names(s) <- gsub("total", categories[i], names(s))
}
s[categories] <- s[categories] / repvec(s$swept.area, ncol = length(categories)) # Standardize by swept area.

# Perform kriging:
res <- NULL
for (i in 1:length(years)){
   # Perform kriging with external drift:
   ix <- (year(s) >= (years[i]-2)) & (year(s) <= years[i])
   m <- ked(s[ix, ], variables = categories, variogram.average = 3, lag = 3, max.distance = 75)
   
   # Calculate abundance:
   tmp <- summary(m, polygon = p)
   tmp$year <- years[i]
   tmp <- tmp[c("year", setdiff(names(tmp), "year"))]
   res <- rbind(res, tmp)
}

# BUGS export:
mu <- matrix(NA, nrow = length(years), ncol = 4)
colnames(mu) <- categories
rownames(mu) <- years
sigma <- mu
for (i in 1:length(categories)){
   mu[, categories[i]]    <- res$mean[res$variable == categories[i]]
   sigma[, categories[i]] <- res$sd[res$variable == categories[i]]
}

# Output 'mu':
cat("N = structure(\n")
cat(paste0("    .Data = c(", paste(round(mu[1, ],1), collapse = ","), ",\n"))
for (i in 2:(nrow(mu)-1)){
   cat(paste0("              ", paste(round(mu[i, ],1), collapse = ","), ",\n"))
}
cat(paste0("              ", paste(round(mu[nrow(mu), ],1), collapse = ","), "), \n"))
cat(paste0("       .Dim = c(", nrow(mu), ",", 4, "),\n"))

cat("sigma = structure(\n")
cat(paste0("    .Data = c(", paste(round(sigma[1, ],1), collapse = ","), ",\n"))
for (i in 2:(nrow(sigma)-1)){
   cat(paste0("              ", paste(round(sigma[i, ],1), collapse = ","), ",\n"))
}
cat(paste0("              ", paste(round(sigma[nrow(sigma), ],1), collapse = ","), "), \n"))
cat(paste0("       .Dim = c(", nrow(sigma), ",", 4, "),\n"))

# BUGS export:
mu <- matrix(NA, nrow = length(years) + length(categories) - 1, ncol = 4)
colnames(mu) <- rev(categories)
sigma <- mu
for (i in 1:length(categories)){
   mu[(5-i):(nrow(mu) - (i-1)), categories[i]]    <- res$mean[res$variable == categories[i]]
   sigma[(5-i):(nrow(sigma) - (i-1)), categories[i]] <- res$sd[res$variable == categories[i]]
}

# Output 'mu':
cat("mu.r = structure(\n")
cat(paste0("    .Data = c(", paste(round(mu[1, ],1), collapse = ","), ",\n"))
for (i in 2:(nrow(mu)-1)){
   cat(paste0("           ", paste(round(mu[i, ],1), collapse = ","), ",\n"))
}
cat(paste0("              ", paste(round(mu[nrow(mu), ],1), collapse = ","), "), \n"))
cat(paste0("       .Dim = c(", nrow(mu), ",", 4, "),\n"))
	
# Output 'sigma':
cat("sigma.r = structure(\n")
cat(paste0("    .Data = c(", paste(round(sigma[1, ],1), collapse = ","), ",\n"))
for (i in 2:(nrow(sigma)-1)){
   cat(paste0("           ", paste(round(sigma[i, ],1), collapse = ","), ",\n"))
}
cat(paste0("              ", paste(round(sigma[nrow(sigma), ],1), collapse = ","), "), \n"))
cat(paste0("       .Dim = c(", nrow(sigma), ",", 4, "),\n"))

