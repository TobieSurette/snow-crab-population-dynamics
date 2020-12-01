library(gulf.data)
library(gulf.graphics)

x <- read.nsslen(2019, survey = "regular", species = "lobster", sex = 1) 

x0 <- seq(0, 160, len = 1000)
grow(x$length, species = "snow crab") # Deterministic.
      
gbarplot(table(x$length), width = 1, border = "grey60", grid = TRUE, xlab = "Carapace width (mm)", ylab = "Frequency")
grid()
lines(x0, grow(x$length, x0, species = "lobster"), lwd = 2, col = "blue")
box()

