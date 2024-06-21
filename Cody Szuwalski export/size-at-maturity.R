library(gulf.data)
library(gulf.stats)
library(mgcv)

years <- 1997:2023

res <- data.frame(cw = seq(0, 150, by = 1))
for (i in 1:length(years)){
   print(years[i])
   b <- read.scsbio(years[i], valid = 1, survey = "regular")
   b$tow.id <- tow.id(b)
   b <- b[which(b$sex == 1), ]
   b <- b[which((b$carapace.width > 0) & (b$carapace.width <= 150)), ]
   b$maturity <- morphometric.maturity(b)
   b$cw <- round(b$carapace.width)
   
   tmp <- aggregate(list(k = b$maturity), by = b["cw"], sum, na.rm = TRUE)
   tmp$n <- aggregate(list(n = b$maturity), by = b["cw"], length)$n
   
   model <- gam(maturity ~ s(cw), data = b, family = "binomial")
   
   gbarplot(tmp$k / tmp$n, tmp$cw)
   
   a <- unique(data.frame(b$cw, predict(model, type = "response")))
   a <- a[order(a[,1]), ]
   lines(a[,1], a[,2], lwd = 2, col = "red2")
   #points(b$cw, predict(model, type = "response"), col = "red")
   
   rr <- data.frame(predict(model, newdata = list(cw = res$cw), type = "response"))
   names(rr) <- years[i]
   res <- cbind(res, rr)
}
