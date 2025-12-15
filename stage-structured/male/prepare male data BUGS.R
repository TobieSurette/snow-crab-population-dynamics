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
b <- read.scsbio(year = (min(years)-2):max(years), survey = "regular", species = 2526, sex = 1)  # Biological data.]

categories <- c("MIGE56L69SC12", "MIGE56L69SC345", "MMGE56L69SC12", "MMGE56L69SC345", 
                "MIGE69L83",     "MIGE69L83SC345", "MMGE69L83SC12", "MMGE69L83SC345", 
                "MIGE83L95",     "MIGE83L95SC345", "MMGE83L95SC12", "MMGE83L95SC345",       
                "MIGE95",        "MIGE95SC345",    "MMGE95SC12",    "MMGE95SC345")    

# Define male categories:
b$category <- ""
for (i in 1:length(categories)){
   print(i)
   b$category[which(is.category(b, categories[i]))] <- categories[i]
}
b <- b[which(b$category != ""), ]

# Import catch data:
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
mu <- matrix(NA, nrow = length(years), ncol = length(categories))
colnames(mu) <- categories
rownames(mu) <- years
sigma <- mu
for (i in 1:length(categories)){
   mu[, categories[i]]    <- res$mean[res$variable == categories[i]]
   sigma[, categories[i]] <- res$sd[res$variable == categories[i]]
}

clc()
for (j in 1:4){
   # Output 'mu':
   cat(paste0("N", 5-j, " = structure(\n"))
   cat(paste0("    .Data = c(", paste(round(mu[1, (1:4) + 4*(j-1)],1), collapse = ","), ",\n"))
   for (i in 2:(nrow(mu)-1)){
      cat(paste0("              ", paste(round(mu[i, (1:4) + 4*(j-1)],1), collapse = ","), ",\n"))
   }
   cat(paste0("              ", paste(round(mu[nrow(mu), (1:4) + 4*(j-1)],1), collapse = ","), "), \n"))
   cat(paste0("       .Dim = c(", nrow(mu), ",", 4, ")),\n\n"))
}

for (j in 1:4){
   cat(paste0("N", 5-j, ".sigma = structure(\n"))
   cat(paste0("    .Data = c(", paste(round(sigma[1, (1:4) + 4*(j-1)],1), collapse = ","), ",\n"))
   for (i in 2:(nrow(sigma)-1)){
      cat(paste0("              ", paste(round(sigma[i, (1:4) + 4*(j-1)],1), collapse = ","), ",\n"))
   }
   cat(paste0("              ", paste(round(sigma[nrow(sigma), (1:4) + 4*(j-1)],1), collapse = ","), "), \n"))
   cat(paste0("       .Dim = c(", nrow(sigma), ",", 4, ")),\n\n"))
}

# Calculate commercial biomass:
import(s, fill = 0) <- catch(b, category = "COM", weight = TRUE, as.hard.shelled = TRUE, units = "t") # Merge catches.
s["COM"] <- 1000000 * s["COM"] / s$swept.area   # Convert to tonnes per km2.

# Perform kriging:
res <- NULL
for (i in 1:length(years)){
  # Perform kriging with external drift:
  ix <- (year(s) >= (years[i]-2)) & (year(s) <= years[i])
  m <- ked(s[ix, ], variables = "COM", variogram.average = 3, lag = 3, max.distance = 75)
  
  # Calculate abundance:
  tmp <- summary(m, polygon = p)
  tmp$year <- years[i]
  tmp <- tmp[c("year", setdiff(names(tmp), "year"))]
  res <- rbind(res, tmp)
}


