library(gulf.data)

years <- 1997:2023

# Load data:
b <- read.scsbio(years, valid = 2, survey = "regular", sex = 2, species = 2526)
b$year <- year(b)

# Identify maturity stages:
b$maturity <- maturity(b)
b$stage <- ""
b$stage[which((b$maturity == "immature") & (b$gonad.colour == 1))] <- "immature"
b$stage[which((b$maturity == "immature") & (b$gonad.colour %in% 2:3))] <- "pubescent"
b$stage[which((b$maturity == "mature") & is.new.shell(b))]  <- "primiparous"
b$stage[which((b$maturity == "mature") & !is.new.shell(b))] <- "multiparous"

# Calculate size frequencies:
fi <- freq(b[b$stage == "immature", ],    by = c("date", "tow.id"))
fp <- freq(b[b$stage == "pubescent", ],   by = c("date", "tow.id"))
fm <- freq(b[b$stage == "primiparous", ], by = c("date", "tow.id"))

# Load tow data and attach size-frequencies:
s <- read.scsset(years, survey = "regular", valid = 1)
si <- sp <- sm <- s
import(si, fill = 0) <- fi
import(sp, fill = 0) <- fp
import(sm, fill = 0) <- fm

# Buffer frequency variables:
fvars <- function(x) return(names(x)[gsub("[0-9]", "", names(x)) == ""])
fvars <- unique(c(fvars(si), fvars(sp), fvars(sm)))
si[setdiff(fvars, names(si))] <- 0
sp[setdiff(fvars, names(sp))] <- 0
sm[setdiff(fvars, names(sm))] <- 0

# Standardize by swept area:
si[fvars] <- si[fvars] / repvec(si$swept.area, ncol = length(fvars))
sp[fvars] <- sp[fvars] / repvec(sp$swept.area, ncol = length(fvars))
sm[fvars] <- sm[fvars] / repvec(sm$swept.area, ncol = length(fvars))

# Summary frequency tables:
ti <- aggregate(si[fvars], by = list(year = year(si)), mean)
tp <- aggregate(sp[fvars], by = list(year = year(sp)), mean)
tm <- aggregate(sm[fvars], by = list(year = year(sm)), mean)

# Aggregated surgvey summaries:
mi <- 1000000 * apply(ti[fvars], 2, mean)
mp <- 1000000 * apply(tp[fvars], 2, mean)
mm <- 1000000 * apply(tm[fvars], 2, mean)

