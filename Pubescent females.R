library(gulf.data)
library(gulf.graphics)

year <- 2020

# Load immature females (y-2)
s <- read.scsset((year-2):(year), valid = 1, survey = "regular")
b <- read.scsbio((year-2):year, sex = 2)

ix <- which(!is.mature(b) & (b$gonad.colour != 3) & (year(b) == (year-2)))  # Immature females year - 2.
iy <- which(!is.mature(b) & (b$gonad.colour == 3) & (year(b) == (year-1)) & (b$carapace.width >= 30))  # Pubescent females year - 1.
iz <- which(is.primiparous(b) & (year(b) == year))  # Primiparous females year.

si <- s[year(s) == year-2, ]
sp <- s[year(s) == year-1, ]
sm <- s[year(s) == year, ]

# Function which returns size-frequency variables:
fvars <- function(x) return(names(x)[gsub("[0-9]", "", names(x)) == ""])

# Attach size-frequencies:
import(si, fill = 0) = freq(b[ix, ], by = key(s), step = 0.5)
import(sp, fill = 0) = freq(b[iy, ], by = key(s), step = 0.5)
import(sm, fill = 0) = freq(b[iz, ], by = key(s), step = 0.5)

clg()
plot(c(0, 80), c(0, 1000), type = "n", xlab = "", ylab = "", xaxs = "i", yaxs = "i")
grid()
lines(as.numeric(fvars(si)), 1000 * apply(si[fvars(si)], 2, mean, na.rm = TRUE), col = "skyblue3", lwd = 2)
lines(as.numeric(fvars(sp)), 1000 * apply(sp[fvars(sp)], 2, mean, na.rm = TRUE), col = "orange", lwd = 2)
lines(as.numeric(fvars(sm)), 1000 * apply(sm[fvars(sm)], 2, mean, na.rm = TRUE), col = "palegreen3", lwd = 2)

legend("topright", 
       legend = c(paste("Immature", year-2),
                  paste("Pubescent", year-1),
                  paste("Primiparous", year)),
       lwd = 2,
       col = c("skyblue3", "orange", "palegreen3"))

       
       
       