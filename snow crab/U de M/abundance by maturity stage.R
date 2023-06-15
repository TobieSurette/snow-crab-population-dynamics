library(gulf.data)
library(gulf.stats)
library(gulf.graphics)
library(gulf.spatial)

# Survey year:
years <- 2000:2022
output <- "results/tables/"
categories <- c("MILE38",       # Immature.
                "MIG38L95",     # Adolescent sub-legal.
                "MIGE95",       # Adolescent legal.
                "MML95SC12",    # Mature new-shelled sub-legal.
                "MMGE95SC12",   # Mature new-shelled legal.
                "MML95SC345",    # Mature old-shelled sub-legal.
                "MMGE95SC345")   # Mature old-shelled legal.

# Load kriging polygons:
p <- read.gulf.spatial("kriging polygons revised")["gulf"]

# Read three years of data (for variogram averaging):
s <- read.scsset(year = (min(years)-2):max(years), survey = "regular", valid = 1) # Tow data.
b <- read.scsbio(year = (min(years)-2):max(years), survey = "regular")            # Biological data.
b$tow.id <- tow.id(b)

# Import catch data:
import(s, fill = 0) <- catch(b, category = categories) # Merge catches.
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

# French:
names <- names(res)
names <- gsub("polygon", "polygone", names)
names <- gsub("mean", "moyenne", names)
names <- gsub("sample", "echantillon", names)
names <- gsub("sd", "ecart.type", names)
names <- gsub("mad", "deviation.absolue", names)
names <- gsub("lowerCI", "int.conf.bas", names)
names <- gsub("upperCI", "int.conf.haut", names)
names(res) <- names
res$description <- category(res$variable, language = "french")
#write.table(res, file = paste0(output, "français/", "abondance.", year, ".csv"), row.names = FALSE, sep = ",")
