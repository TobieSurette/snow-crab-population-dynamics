library(MASS)

logistique = function(x, k, c) 1/(1+exp(c*(k-x)))

VII = dnorm(intervalle, mean = mu[4], sd = sigma[4])

d_VIII_ado = function(taille,k,c) {
  logistique(taille,k,c)*dnorm(taille, mean = a*mu[4]+b, sd = sqrt(variance[4]*a^2))
}

d_VIII_imm = function(taille,k,c) {
  (1-logistique(taille,k,c))*dnorm(taille, mean = a*mu[4]+b, sd = sqrt(variance[4]*a^2))
}

Trouver_moyenne = function(mu, fonction){
  total = integrate(fonction, lower = 1, upper = 95)$value
  return((0.5-(integrate(fonction, lower = 1, upper = mu)$value/total))^2)
}

Par_logistique = function(x){
  moy_VIII_imm = optimize(function (m) Trouver_moyenne(m, fonction = function(y) d_VIII_imm(y,x[1],x[2])), intervalle)$minimum
  moy_VIII_ado = optimize(function (m) Trouver_moyenne(m, fonction = function (y) d_VIII_ado(y,x[1],x[2])), intervalle)$minimum
  return((moy_VIII_imm-37.5)^2
         +(moy_VIII_ado-41.7)^2)
}

theta=c(1,1)

par_log = optim(theta,Par_logistique)

# library(MASS)
# 
# VII = rnorm(10000, mean = mu[4], sd = sigma[4])
# 
# VIII_imm = c()
# VIII_ado = c()
# Choix = rep(0,length(VII))
# 
# logistique = function(x, moitie, k) 1/(1+exp(k*(moitie-x)))
# 
# groupes = c("Immature", "Adolescent")
# 
# Croissance_imm_ado = function(x,k){
#   for(i in 1:length(VII)){
#     Choix[i] = sample(groupes, size = 1, 
#                    prob = c(1-logistique(VII[i],x,k), logistique(VII[i],x,k)), 
#                    replace = TRUE)
#   }
#   return(Choix)
# }
# 
# F = function(x){
#   Choix=Croissance_imm_ado(x[1],x[2])
#   for(i in 1:length(Choix)){
#     if(Choix[i] == "Immature"){
#       VIII_imm = c(VIII_imm, a*VII[i]+b)
#     }else{
#       VIII_ado = c(VIII_ado, a*VII[i]+b)
#     }
#   }
#   
#   return((mean(VIII_imm)-mu[5])^2+(mean(VIII_ado)-41.7)^2)
# }
# 
# par=c(mu[4],1)
# 
# Par_logistique = optim(par,F)