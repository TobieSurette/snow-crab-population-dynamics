library(gulf.graphics)

M <- matrix(NA, nrow = 16, ncol = 16)

colnames(M) <- rownames(M) <- apply(expand.grid(c("I", "S", "N", "O"), 4:1), 1, paste, collapse = "")

M["I4", c("S4")] <- 1
M["I4", c("I3", "N3")] <- 2
M["S4", c("N3")] <- 2
M[c("N4", "O4"), c("O4")] <- 1

M["I3", c("S3")] <- 1
M["I3", c("I2", "N2")] <- 2
M["S3", c("N2")] <- 2
M[c("N3", "O3"), c("O3")] <- 1

M["I2", c("S2")] <- 1
M["I2", c("I1", "N1")] <- 2
M["S2", c("N1")] <- 2
M[c("N2", "O2"), c("O2")] <- 1

M["I1", c("S1")] <- 1
M["I1", c("I1", "N1")] <- 2
M["S1", c("N1")] <- 2
M[c("N1", "O1"), c("O1")] <- 1

png("stage-structured/male/male transition matrix figure.png",
    width = 7, height = 7, units = "in", res = 500)
image(1:16, 1:16, t(M[rev(1:nrow(M)), ]), col = fade(c("grey20", "green3")),
      xlab = "", ylab = "", xaxt = "n", yaxt = "n")
hline((1:16)-0.5, col = fade("grey50"), lwd = 0.5)
vline((1:16)-0.5, col = fade("grey50"), lwd = 0.5)
hline((1:3)*4 + 0.5, lwd = 1.5)
vline((1:3)*4 + 0.5, lwd = 1.5)

axis(1, at = (1:4)*4 - 1.5, labels = paste0("R-", 4:1), tick = FALSE, font = 2, cex.axis = 1.25)
axis(2, at = (4:1)*4 - 1.5, labels = paste0("R-", 4:1), tick = FALSE, font = 2, cex.axis = 1.25)

axis(1, at = 1:16, labels = rownames(M), tick = FALSE, cex.axis = 0.85, padj = -2)
axis(2, at = 1:16, labels = rev(apply(expand.grid(c("I", "S", "N", "O"), 4:1), 1, paste, collapse = "")), 
     tick = FALSE, cex.axis = 0.85, padj = 1.5)

box(col = "grey50")
dev.off()


