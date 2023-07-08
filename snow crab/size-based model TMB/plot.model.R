plot.model <- function(obj, data, type = "maturity"){
   
   # Function to parse parameter vector:
   tolist <- function(x){
      names <- unique(names(x))
      r <- list(); for (i in 1:length(names)) r[[i]] <- as.numeric(x[names[i] == names(x)])
      names(r) <- names
      return(r)
   }
   
   # Parse parameter vector:
   p <- tolist(obj$par)
   
   # Generate report:
   years <- as.numeric(dimnames(data$f_imm)$year)
   r <- obj$report()
   dimnames(r$n_imm) <- list(year = (min(years)-1):max(years), cw = colnames(data$f_imm))
   dimnames(r$n_mat) <- list(year = (min(years)-1):max(years), cw = colnames(data$f_imm))
   
   if ("maturity" %in% type){
      plot(c(0, 140), c(0, 1), type = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "")
      grid()
      lines(data$x, r$p_maturity, type = "l", lwd = 2)
      mtext("Proportion", 2, 2.5, cex = 1.25)
      mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
      
      # Show parameter values:
      yy <- approx(data$x, r$p_maturity, p$xp_maturity)$y
      lines(rep(p$xp_maturity[1], 2), c(0, yy[1]), lty = "dashed", col = "red", lwd = 2)
      lines(c(0, p$xp_maturity[1]), c(yy[1], yy[1]), lty = "dashed", col = "red", lwd = 2)
      lines(rep(p$xp_maturity[2], 2), c(0, yy[2]), lty = "dashed", col = "red", lwd = 2)
      lines(c(0, p$xp_maturity[2]), c(yy[2], yy[2]), lty = "dashed", col = "red", lwd = 2)
      lines(c(0, mean(p$xp_maturity)), rep(1 / (1 + exp(p$logit_p_mix_maturity)), 2), lty = "dashed", col = "red", lwd = 2)
      
      text(p$xp_maturity[1] / 2, yy[1], paste0("p = ", round(yy[1],2)), pos = 3)
      text(p$xp_maturity[1] / 2, 2*yy[1], paste0("p = ", round(2*yy[1],2)), pos = 3)
      text(p$xp_maturity[2] / 2, yy[2], paste0("p = ", round(yy[2],2)), pos = 3)
      
      axis(1, at = round(p$xp_maturity,1), col = "red", col.axis = "red")
      axis(1, at = round(p$xp_maturity[1],1), col = "red", col.axis = "red")
      mtext("Proportion of crab that moult to maturity", 3, 1.0, cex = 1.25)
      
      box(col = "grey60")
   }
   
   if ("selectivity" %in% type){
      plot(c(0, 140), c(0, 1), type = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "")
      grid()
      lines(data$x, r$p_selectivity, type = "l", lwd = 2)
      mtext("Proportion", 2, 2.5, cex = 1.25)
      mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
      
      # Show parameter values:
      xx <- approx(r$p_selectivity, data$x, 0.5)$y
      lines(rep(xx, 2), c(0, 0.5), lty = "dashed", col = "red", lwd = 2)
      lines(c(0, xx), c(0.5, 0.5), lty = "dashed", col = "red", lwd = 2)
      
      axis(1, at = round(xx,1), col = "red", col.axis = "red")
      axis(2, at = 0.5, col = "red", col.axis = "red")
      mtext("Survey trawl selectivity", 3, 1.0, cex = 1.25)
      
      box(col = "grey60")
   }
   
   if ("skip" %in% type){
      plot(c(0, 140), c(0, 1), type = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "")
      grid()
      lines(data$x, r$p_skp, type = "l", lwd = 2)
      mtext("Proportion", 2, 2.5, cex = 1.25)
      mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
      
      # Show parameter values:
      pp <- 1 / (1 + exp(-p$logit_p_max_skp))
      xx <- approx(r$p_skp, data$x, pp / 2)$y
      lines(rep(xx, 2), c(0, pp / 2), lty = "dashed", col = "red", lwd = 2)
      lines(c(0, xx), c(pp / 2, pp / 2), lty = "dashed", col = "red", lwd = 2)
      
      axis(1, at = round(xx,1), col = "red", col.axis = "red")
      axis(2, at = round(pp / 2,2), col = "red", col.axis = "red")
      mtext("Proportion of crab that will skip a moult", 3, 1.0, cex = 1.25)
      
      box(col = "grey60")
   }
   
   if ("recruitment" %in% type){
      gbarplot(exp(p$log_scale_recruitment_mu + p$log_scale_recruitment) / 1000, years, xaxt = "n", grid = TRUE)
      
      lines(data$x, r$p_skp, type = "l", lwd = 2)
      mtext("Recruitment scale", 2, 2.5, cex = 1.25)
      mtext("Year", 1, 2.5, cex = 1.25)
      
      axis(1, at = seq(min(years), max(years), by = 4))
      axis(1, at = seq(min(years)+2, max(years), by = 4))
      box(col = "grey60")
   }
   
   if ("growth" %in% type){
      plot(c(0, 140), c(0, 25), type = "n", xaxs = "i", yaxs = "i", xlab = "", ylab = "")
      grid()
      lines(data$x, r$mu_growth, type = "l", lwd = 2)
      lines(data$x, r$mu_growth - r$sigma_growth, lty = "dashed", lwd = 1.5, col = "grey50")
      lines(data$x, r$mu_growth + r$sigma_growth, lty = "dashed", lwd = 1.5, col = "grey50")
      lines(data$x, r$mu_growth - 1.96 * r$sigma_growth, lty = "dotted", lwd = 1.0, col = "grey50")
      lines(data$x, r$mu_growth + 1.96 * r$sigma_growth, lty = "dotted", lwd = 1.0, col = "grey50")
      
      mtext("Growth increment (mm)", 2, 2.5, cex = 1.25)
      mtext("Pre-moult carapace width (mm)", 1, 2.5, cex = 1.25)
      
      yy <- approx(data$x, r$mu_growth, p$xp_growth)$y
      lines(rep(p$xp_growth, 2), c(0, yy), lty = "dashed", col = "red", lwd = 2)
      
      box(col = "grey60")
   }
   
   if ("M" %in% type){
      plot(c(min(years)-0.5, max(years)+0.5), c(0, 0.6), type = "n", xaxs = "i", xaxt = "n", yaxs = "i", xlab = "", ylab = "")
      grid()
      
      lines(years, 1/(1+ exp(-p$log_mu_M_imm - p$log_dev_M_imm)), lty = "dashed", col = "green", lwd = 2)
      lines(years, 1/(1+ exp(-p$log_mu_M_mat - p$log_dev_M_mat)), lty = "solid", col = "blue", lwd = 2)

      mtext("Proportion", 2, 2.5, cex = 1.25)
      mtext("Year", 1, 2.5, cex = 1.25)
      
      axis(1, at = seq(min(years), max(years), by = 4))
      axis(1, at = seq(min(years)+2, max(years), by = 4))
      
      legend("topright", 
             legend = c("Mature", "Immature"), 
             lwd = 2, 
             col = c("blue", "green"),
             lty = c("solid", "dashed"),
             cex = 1.25)
      
      mtext("Natural mortality", 3, 1.0, cex = 1.25)
      
      box(col = "grey60")
   }
   
   if ("F" %in% type){
      plot(c(min(years)-0.5, max(years)+0.5), c(0, 1), type = "n", xaxs = "i", xaxt = "n", yaxs = "i", xlab = "", ylab = "")
      grid()
      
      lines(years, 1 / (1+ exp(-exp(-p$log_mu_F - p$log_dev_F))), lwd = 2)

      mtext("Proportion", 2, 2.5, cex = 1.25)
      mtext("Year", 1, 2.5, cex = 1.25)
      
      axis(1, at = seq(min(years), max(years), by = 4))
      axis(1, at = seq(min(years)+2, max(years), by = 4))
      
      mtext("Maximum fishing mortality", 3, 1.0, cex = 1.25)
      
      box(col = "grey60")
   }
   
   if ("lf" %in% type){
      m <- kronecker(matrix(1:21, ncol = 3), matrix(1, ncol = 5, nrow = 5))
      m <- rbind(0, cbind(0, 0, m, 0), 0, 0)
      layout(m)
      par(mar = c(0,0,0,0))
      
      for (i in 2:length(years)){
         year <- as.character(years[i])
         
         gbarplot(cbind(data$f_imm[year, ], data$f_mat[year, ]) , data$x, 
                  legend = FALSE, grid = TRUE, xaxs = "i", xlim = c(0, 140), yaxs = "i", 
                  ylim = c(0, 500), xaxt = "n", yaxt = "n")
         lines(data$x, r$n_imm[year, ] + r$n_mat[year, ], lwd = 2, col = "blue")
         lines(data$x, r$n_imm[year, ], lwd = 2, col = "red")
         
         if (i %in% 2:8) axis(2)
         if (i %in% (c(7, 14, 21)+1)) axis(1)
         text(par("usr")[1] + 0.8 * diff(par("usr")[1:2]),
              par("usr")[3] + 0.8 * diff(par("usr")[3:4]), year)
         
         vline(c(9.4, 14, 20.3, 28.5, 38), lty = "dashed")
         box(col = "grey60")
         
         if (year == 2005) mtext("Density (n/km2)", 2, 2.5, cex = 1.25)
         if (year == 2015) mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
      }
   } 
}
