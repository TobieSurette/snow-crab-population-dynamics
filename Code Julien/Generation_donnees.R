library(MASS)

set.seed(2024) # Résultats reproductibles

CDF_G=function(x){
  integrate(G,w=w,mu=mu,log_sigma=log_sigma,lower=1,upper=x)$value
}

inverse = function (f, lower = -100, upper = 100) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}


Genere_immature = function(n) {
  Inv_CDF_G=inverse(CDF_G)

  random = runif(n)
  immature = rep(0,length(random))

  for(i in 1:length(random)){
    immature[i]=Inv_CDF_G(random[i])[[1]]
  }
  return(immature)
}

immature = Genere_immature(10000)

hist(immature, 
     main = "Nombre de crabes pour chaque taille sur 10000 crabes",
     xlab = "Largeur de carapace en mm",
     ylab = "# de crabes",
     breaks = 200)




