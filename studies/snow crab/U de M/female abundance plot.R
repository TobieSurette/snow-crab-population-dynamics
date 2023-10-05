


x <- read.csv("snow crab/U de M/female abundance 2000-2022.csv")


plot(x$year, 0.6 * x$immature, type = "l", lwd = 2, xlim = c(1997.5, 2022.5), ylim = c(0, 400), xaxs = "i", yaxs = "i")
grid()

lines(x$year-1, x$adolescent, col = "green3", lwd = 2)
lines(x$year-2, x$primiparous, col = "brown3", lwd = 2)
lines(x$year-4, 0.35 * x$multiparous, col = "skyblue3", lwd = 2)

