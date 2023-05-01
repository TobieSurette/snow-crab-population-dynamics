#' @title Instar Analysis
#' 
#' @description Perform instar analysis or classification.
#' 
#' @param x Target object.
#' 
#' @examples 

# Parse parameter vector and return instar moments and proportions:
instar.stats <- function(theta, species = 2526, sex = 2, maturity, log = TRUE){
     
   # Parse parameter vector:
   if (species == 2526){
      if (missing(theta)){
         theta <- c(mu4 = 2.293041, 
                    log.sigma = -2.3685149, 
                    log.sigma.pubescent = -1.322312266,
                    log.sigma.mature = -5,
                    a = 0.890443829, b = 0.645895530,            
                    a.pubescent = 0.663637366, b.pubescent = 1.515007633,          
                    a.mature = 1.040998406, b.mature = 0.004727189)
      }
     
      if (sex == 2){
         # Immature (instars I to XI):
         mu.immature <- rep(NA, 11)
         ix <- grep("mu[0-9]+", names(theta))
         iy <- as.numeric(gsub("mu", "", names(theta)[ix])) 
         mu.immature[iy] <- as.numeric(theta[ix])
         for (i in (iy+1):length(mu.immature)) mu.immature[i] <- theta["a"] * mu.immature[i-1] + theta["b"]
         for (i in (iy-1):1) mu.immature[i] <- (mu.immature[i+1] - theta["b"]) / theta["a"]
         names(mu.immature) <- as.roman(1:length(mu.immature))
         sigma.immature <- rep(exp(theta["log.sigma"]), length(mu.immature))
         names(sigma.immature) <- names(mu.immature)
         
         # Pubescent (instars VII to IX):
         mu.pubescent[1] <- theta["a.pubescent"] * ((mu.immature[1] - theta["b"]) / theta["a"]) + theta["b.pubescent"]
         mu.pubescent[2:11] <- theta["a.pubescent"] * mu.immature[as.character(as.roman(1:10))] + theta["b.pubescent"]
         names(mu.pubescent) <- as.character(as.roman(1:11))
         sigma.pubescent <- sigma.immature
        
         # Mature (instars VIII to X):
         mu.mature[1] <- theta["a.mature"] * ((mu.pubescent[1] - theta["b.pubescent"]) / theta["a.pubescent"]) + theta["b.mature"]
         mu.mature[2:11] <- theta["a.mature"] * mu.pubescent[as.character(as.roman(1:10))] + theta["b.mature"]
         names(mu.mature) <- as.character(as.roman(1:11))
         sigma.mature <- sigma.immature
         
         if (!log){
            # Immature:
            m <- exp(mu.immature + 0.5 * sigma.immature^2)
            s <- sqrt((exp(sigma.immature^2) - 1) * exp(2 * mu.immature + sigma.immature^2))
            mu.immature <- m
            sigma.immature <- s
            
            # Pubescent:
            m <- exp(mu.pubescent + 0.5 * sigma.pubescent^2)
            s <- sqrt((exp(sigma.pubescent^2) - 1) * exp(2 * mu.pubescent + sigma.pubescent^2))
            mu.pubescent <- m
            sigma.pubescent <- s
            
            # Mature:
            m <- exp(mu.mature + 0.5 * sigma.mature^2)
            s <- sqrt((exp(sigma.mature^2) - 1) * exp(2 * mu.mature + sigma.mature^2))
            mu.mature <- m
            sigma.mature <- s
         }
      }    
   }
  
   # Compile results:
   r <- list(immature  = list(mu = mu.immature, sigma = sigma.immature),
             pubescent = list(mu = mu.pubescent, sigma = sigma.pubescent),
             mature    = list(mu = mu.mature, sigma = sigma.mature))

   # Parse immature mixture proportions:
   ix <- grep("logit.p.immature", names(theta))
   if (length(ix) > 0){
      iy <- as.numeric(gsub("logit.p.immature", "", names(theta)[ix]))
      p.immature <- exp(theta[ix]) / (1 + sum(exp(theta[ix])))
      p.immature[length(p.immature)+1] <- 1 - sum(p.immature)
      names(p.immature) <- as.roman(c(iy, max(iy)+1))
      r$immature$p     <- p.immature
      r$immature$mu    <- r$immature$mu[intersect(names(r$immature$mu), names(p.immature))]
      r$immature$sigma <- r$immature$sigma[intersect(names(r$immature$sigma), names(p.immature))]
   }
   
   # Parse pubescent mixture proportions:
   ix <- grep("logit.p.pubescent", names(theta))
   if (length(ix) > 0){
      iy <- as.numeric(gsub("logit.p.pubescent", "", names(theta)[ix]))
      p.pubescent <- exp(theta[ix]) / (1 + sum(exp(theta[ix])))
      p.pubescent[length(p.pubescent)+1] <- 1 - sum(p.pubescent)
      names(p.pubescent) <- as.roman(c(iy, max(iy)+1))
      r$pubescent$p     <- p.pubescent  
      r$pubescent$mu    <- r$pubescent$mu[intersect(names(r$pubescent$mu), names(p.pubescent))]
      r$pubescent$sigma <- r$pubescent$sigma[intersect(names(r$pubescent$sigma), names(p.pubescent))]
   }
   
   # Parse mature mixture proportions:
   if (length(ix) > 0){
      ix <- grep("logit.p.mature", names(theta))
      iy <- as.numeric(gsub("logit.p.mature", "", names(theta)[ix]))
      p.mature <- exp(theta[ix]) / (1 + sum(exp(theta[ix])))
      p.mature[length(p.mature)+1] <- 1 - sum(p.mature)
      names(p.mature) <- as.roman(c(iy, max(iy)+1)) 
      r$mature$p     <- p.mature
      r$mature$mu    <- r$mature$mu[intersect(names(r$mature$mu), names(p.mature))]
      r$mature$sigma <- r$mature$sigma[intersect(names(r$mature$sigma), names(p.mature))]
   }
   
   return(r)
}

# Mixture probability density function for instars:
pdf.instar <- function(x, p, mu, sigma, theta, maturity, sum = TRUE){
   if (!missing(theta)){
      r <- instar.stats(theta, species = 2526, sex = 2)[[maturity]]
   }else{
      r <- list(p = p, mu = mu, sigma = sigma)
   } 
  
   v <- NULL
   for (i in 1:length(r$p)){
      d <- r$p[i] * dnorm(x, r$mu[i], r$sigma[i])
      d[is.na(d)] <- 0
      v <- cbind(v, d)
   } 
   colnames(v) <- names(r$p)
   rownames(v) <- rownames(x)
   
   return(v)
}

# Instar mixture log-likelihood function:
loglike.instar <- function(theta, xi, xp, xm, fixed){
   if (!missing(fixed)) theta <- c(theta, fixed)
  
   # Parse parameter vector:   
   r <- instar.stats(theta, species = 2526, sex = 2)
  
   ll <- 0 # Initialize log-likelihood function.
   
   # Immature mixture density:
   vi <- pdf.instar(as.numeric(names(xi)), theta = theta, maturity = "immature")
   ll <- ll - sum(as.numeric(xi) * log(apply(vi, 1, sum))) 
  
   # Pubescent:
   if (!missing(xp)){
     vp <- pdf.instar(as.numeric(names(xp)), theta = theta, maturity = "pubescent")
     ll <- ll - sum(as.numeric(xp) * log(apply(vp, 1, sum))) 
   }
  
   # Mature:
   if (!missing(xm)){
     vm <- pdf.instar(as.numeric(names(xm)), theta = theta, maturity = "mature")
     ll <- ll - sum(as.numeric(xm) * log(apply(vm, 1, sum))) 
   }
  
   return(ll)
}

# Instar mixture model fitting function:
fit.instar <- function(x, species, sex){
  
  
}

#' @export
instar <- function(x, ...) UseMethod("instar")

#' @export
instar.default <- function(x, group, probability = TRUE, ...){
   if (missing(group)){
      p <- pdf.instar(x, ...)
      p <- p / gulf.utils::repvec(apply(p, 1, sum), ncol = ncol(p))
   }else{
      u <- unique(group)
      r <- NULL
      iy <- NULL
      for (i in 1:length(u)){
         ix <- which(group == u[i])
         p <- pdf.instar(x[ix], ...)
         p <- apply(p, 2, sum) / sum(p)
         p <- pdf.instar(x, p = p, ...)
         if (is.null(r)){
            r <- p
         }else{
            r[setdiff(names(p), names(r))] <- 0
            p[setdiff(names(r), names(p))] <- 0
            p <- p[names(r)]
            r <- rbind(r, p)
         } 
         iy <- c(iy, ix)
      }
      p <- matrix(NA, nrow = length(x), ncol = ncol(p))
      p[iy,] <- r
   }

   return(p)
}

#' @export
instar.scsbio <- function(x, probability = TRUE, ...){
 
}

