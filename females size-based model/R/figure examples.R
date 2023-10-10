

png(file = "C:/Users/SuretteTJ/Desktop/gulf-population-modelling/studies/females size-based model/figures/female instar example.png",
     units = "in", res = 300, height = 6, width = 7.5)

gbarplot(table(round(5 * b$carapace.width[b$stage == "immature"]) / 5),
         xlim = c(0, 65), xaxs = "i", border = "grey50", col = "grey75", grid = TRUE, lwd = 0.5)


mtext("Carapace width (mm)", 1, 2.25, cex = 1.25, font = 2)
mtext("Frequency", 2, 2.25, cex = 1.25, font = 2)

mu <- c(10.5, 14.75, 20.5, 28.5, 37.25, 48.5)
names(mu) <- as.roman(4:9)
text(mu, c(40, 185, 600, 1050, 810, 70), names(mu), font = 2, cex = 1.1)   

box(col = "grey50")

dev.off()
