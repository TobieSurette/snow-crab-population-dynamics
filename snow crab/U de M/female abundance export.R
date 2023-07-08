library(gulf.data)
library(gulf.stats)
library(gulf.graphics)
library(gulf.spatial)

# Survey year:
years <- 2000:2022
output <- "results/tables/"

# Load kriging polygons:
p <- read.gulf.spatial("kriging polygons revised")["gulf"]

# Read three years of data (for variogram averaging):
s <- read.scsset(year = (min(years)-2):max(years), survey = "regular", valid = 1) # Tow data.
b <- read.scsbio(year = (min(years)-2):max(years), survey = "regular")            # Biological data.
b$tow.id <- tow.id(b)

# Partition female categories:
b$FI <- FALSE
b$FI[which(is.category(b, "FI") & (b$gonad.colour != 3))] <- TRUE # Immature.
b$FA <- FALSE
b$FA[which(is.category(b, "FI") & (b$gonad.colour == 3))] <- TRUE # Adolescent.
b$FP <- FALSE
b$FP[which(is.category(b, "FP"))] <- TRUE # Primiparous.
b$FM <- FALSE
b$FM[which(is.category(b, "FMULT"))] <- TRUE # Multiparous.

categories <- c("FI", "FA", "FP", "FM")

# Import catch data:
import(s, fill = 0) <- aggregate(b[categories], by = b[key(s)], sum) # Merge catches.
s[categories] <- s[categories] / repvec(s$swept.area, ncol = length(categories)) # Standardize by swept area.

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

# write.csv(res, file = paste0(categories, " ", min(years), "-", max(years), ".csv"))

# Add description:
res <- cbind(res["variable"],
             data.frame(description = category(res$variable, language = "english")),
             res[setdiff(names(res), "variable")])

# Output:

# English:
#write.table(res, file = paste0(output, "english/", "abundance.", year, ".csv"), row.names = FALSE, sep = ",")

