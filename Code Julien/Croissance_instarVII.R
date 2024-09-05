library(MASS)

Trouver_par_VIII = function(x){
  m=c(x["mut"],x["muA"])
  
  l_s=c(x["log_sigmat"],x["log_sigmaA"])
  v=exp(2*l_s)
  
  w=x["p"]
  
  p=(exp(w)/(1+exp(w)))
  
  f=(6*(((m[1]-(1-p)*m[2])/p) -mu[5])^2
     +(((v[1]-(1-p)*v[2]-p*(1-p)*(mu[5]-m[2])^2)/p) -variance[5])^2
     +0.1*(p-0.8)^2)
  
  return(f)
}

theta = c(mut = 39, muA = 41,
          log_sigmat = 1.39, log_sigmaA = 1.39, p = 0.5)

par_VIII = optim(theta, Trouver_par_VIII)

mu_VIIIado = par_VIII$par["muA"]
mu_VIIItotal = par_VIII$par["mut"]

log_sigma_VIIIado = par_VIII$par["log_sigmaA"]
log_sigma_VIIItotal = par_VIII$par["log_sigmat"]

sigma_VIIIado = exp(log_sigma_VIIIado)
sigma_VIIItotal = exp(log_sigma_VIIItotal)

poid_VIII = exp(par_VIII$par["p"])/(1+exp(par_VIII$par["p"]))

plot(dnorm(intervalle, mean = mu_VIIItotal, sd = sigma_VIIItotal), type = "l",
     xlab = "Largeur en mm", ylab = "", ylim = c(0,0.1))
lines((1-poid_VIII)*dnorm(intervalle, mean = mu_VIIIado, sd = sigma_VIIIado), col = "blue")
lines(poid_VIII*dnorm(intervalle, mean = mu[5], sd = sigma[5]), col = "red")
abline(v=mu[5],lty=2, col="red")
abline(v=mu_VIIIado,lty=2, col="blue")
abline(v=mu_VIIItotal,lty=2, col="black")
legend(x="right", legend = c("immature", "ado", "total"), fill = c("red","blue","black"), cex=0.6)

# logistique = function(x, k=34, c=1) 1/(1+exp(c*(k-x)))
# 
# meanCW_VIIIado = data.frame(x=c(10,14.5,20.5,28.4),y=c(14.5,20.5,28.4,41.7))
# Croissance_VIIIado = lm(y~x,meanCW_VIIIado)
# 
# a_VIIIado = Croissance_VIIIado$coefficients[2]
# b_VIIIado = Croissance_VIIIado$coefficients[1]
# 
# VII = dnorm(intervalle, mean = mu[4], sd = sigma[4])
# plot(VII, type = "l")
# lines(logistique(intervalle)*VII, col="blue")
# lines((1-logistique(intervalle))*VII, col="red")
# 
# d_VIII_ado = function(taille) {
#   logistique(taille)*dnorm(taille, mean = a*mu[4]+b, sd = sqrt(variance[4]*a^2))
# }
# 
# d_VIII_imm = function(taille) {
#   (1-logistique(taille))*dnorm(taille, mean = a*mu[4]+b, sd = sqrt(variance[4]*a^2))
# }
# 
# VIII_ado = d_VIII_ado(intervalle)
# VIII_imm = d_VIII_imm(intervalle)
# 
# plot(dnorm(intervalle, mean = a*mu[4]+b, sd = sqrt(variance[4]*a^2)), type="l")
# lines(VIII_imm, col="red")
# lines(VIII_ado, col="blue")
# legend(x="right", legend=c("VIII imm", "VIII ado"), fill=c("red","blue"),cex=0.7)
# 
# Trouver_moyenne = function(mu, fonction){
#   total = integrate(fonction, lower = 1, upper = 95)$value
#   return((0.5-(integrate(fonction,lower = 1, upper = mu)$value/total))^2)
# }
# 
# moy_VIII_imm=optimize(Trouver_moyenne,intervalle,
#              fonction=d_VIII_imm)$minimum
# 
# moy_VIII_ado=optimize(Trouver_moyenne,intervalle,
#                       fonction=d_VIII_ado)$minimum


# VII = rnorm(100000, mean = mu[4], sd = sigma[4])
# 
# VIII_imm = c()
# VIII_ado = c()
# 

# 
# groupes = c("Immature", "Adolescent")
# 
# for(i in 1:length(VII)){
#   Choix = sample(groupes, size = 1, 
#          prob = c(1-logistique(VII[i]), logistique(VII[i])), 
#          replace = TRUE)
# 
#   if(Choix == "Immature"){
#       VIII_imm = c(VIII_imm, a*VII[i]+b)
#   }else{
#     VIII_ado = c(VIII_ado, a*VII[i]+b)
#   }
# }
# 
# plot(VIII_imm, col="blue",ylim=c(min(VIII_imm), max(VIII_ado)))
# points(VIII_ado, col="red")
# 
# hist(VIII_imm, breaks = 200)
# hist(VIII_ado, breaks = 200)
# 
# hist(c(a*VIII_imm+b, VIII_ado), breaks = 200)
# 
