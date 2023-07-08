library(gulf.data)
library(gulf.spatial)
library(gulf.stats)

# Generate instar abundance summary:
tab <- aggregate(list(n_imm = b$maturity == "immature"), by = b[c("date", "tow.id", "group")], sum)
tab$n_pub <- aggregate(list(n = b$maturity == "pubescent"), by = b[c("date", "tow.id", "group")], sum)$n
tab$n_mat <- aggregate(list(n = b$maturity == "mature"), by = b[c("date", "tow.id", "group")], sum)$n
n_imm <- as.data.frame(r$p_imm * repvec(tab$n_imm, ncol = ncol(r$p_imm)))
names(n_imm) <- paste0("n_instar_", 1:ncol(r$p_imm), "_imm")
n_pub <- as.data.frame(r$p_pub * repvec(tab$n_pub, ncol = ncol(r$p_pub)))
names(n_pub) <- paste0("n_instar_", 1:ncol(r$p_pub), "_pub")
n_mat <- as.data.frame(r$p_mat * repvec(tab$n_mat, ncol = ncol(r$p_mat)))
names(n_mat) <- paste0("n_instar_", 1:ncol(r$p_mat), "_mat")
tab <- cbind(tab, n_imm, n_pub, n_mat)
fvars <- names(tab)[grep("instar", names(tab))]
fvars <- fvars[which(apply(tab[fvars], 2, sum) > 0)]
tab <- tab[c("date", "tow.id", "group", "n_imm", "n_pub", "n_mat", fvars)]
ix <- match(tab[key(s)], s[key(s)])
tab$grid <- s$grid[ix]
tab$longitude <- lon(s)[ix]
tab$latitude  <- lat(s)[ix]
tab$swept.area <- s$swept.area[ix]
fvars <- names(tab)[grep("instar", names(tab))]
tab <- tab[!is.na(tab$grid), ]
tab$year <- year(tab)
tab <- aggregate(tab[fvars], by = tab[c("date", "year", "tow.id", "group", "grid")], mean)

# Merge with tow data:
ix <- match(tab[key(s)], s[key(s)])
vars <- names(tab)[grep("n_inst", names(tab))]
s[vars] <- 0
s[ix, vars] <- tab[vars]

# Load kriging polygon:
p <- read.gulf.spatial("kriging polygons revised")["gulf"]

res <- NULL
for (i in 1:length(years)){
   print(years[i])
   ss <- s[year(s) == years[i], ]
   ss[vars] <- 1000000 * ss[vars] / repvec(ss$swept.area, ncol = length(vars))

   # Perform kriging:
   m <- ked(ss, variables = vars, variogram.average = 1, lag = 3, max.distance = 75)
   
   # Calculate abundance:
   tmp <- summary(m, polygon = p)
   tmp$year <- years[i]
   tmp <- tmp[c("year", setdiff(names(tmp), "year"))]
   res <- rbind(res, tmp)
}

write.csv(res, file = paste0("instar/Instar abundances ", min(years), "-", max(years), ".csv"), row.names = FALSE)


