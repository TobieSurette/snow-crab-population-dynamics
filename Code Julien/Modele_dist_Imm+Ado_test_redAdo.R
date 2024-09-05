library(MASS)
library(readxl)

rm(list=ls()) # effacer la mémoire d el'environnement de travail
ls() # vérifier la liste des objets
set.seed(2024) # Résultats reproductibles

#######################################################################
#
# Partie Régression (taille après en fonction de taille avant)
#
#######################################################################

# Tableau des moyennes de tailles avant (x) et après (y) la mue
MeanCW = data.frame(x=c(10, 14.5, 20.5), y=c(14.5, 20.5, 28.4))
MeanCWAdo = data.frame(x=c(29.55, 37.5, 49.1), y=c(41.7, 49.1, 58.1))

# Régression linéaire (y=ax+b) de la relation entre la taille avant et
# après la mue
Croissance = lm(y~x, data = MeanCW)
CroissanceAdo = lm(y~x,dat = MeanCWAdo)

# Graphique des points réelles avant et après mue
plot(c(10,14.5,20.5,28.4,37.5,49.1), c(14.5,20.5,28.4,37.5,49.1,58.1), type="o", xlim=c(0,100), ylim=c(0,100), 
     xlab = "Taille avant", ylab = "Taille après",col="orange")
points(c(10,14.5,20.5,28.4,37.5), c(14.5,20.5,28.4,37.5,49.1), type="o", col="red")
points(c(10,14.5,20.5,28.4), c(14.5,20.5,28.4,41.7), type="o", col="gold2")

# Pente de la régression
a=Croissance$coefficients[2]
a_ado=CroissanceAdo$coefficients[2]

# Ordonnée de la régression
b=Croissance$coefficients[1]
b_ado=CroissanceAdo$coefficients[1]

plot(MeanCWAdo$x,MeanCWAdo$y,"o",xlim = c(0,60), ylim = c(0,60))
# Droite de régression
axe_x = seq(0,100,length=100)
lines(axe_x, predict(CroissanceAdo, data.frame(x=axe_x)), col="blue")
grid(nx=NULL, ny=NULL)

#######################################################################
#
# Extraction et visualisation des données
#
#######################################################################


# On extrait les données du excel
donneesImm=read_excel("~/Crabe/female densities by size 1997-2023.xlsx")
donneesAdo=read_excel("~/Crabe/female densities by size 1997-2023.xlsx", 
                      sheet = "Adolescent")

# Calcul de la densité moyenne des crabes de chaque taille
# à travers les années
data_moyImm=colMeans(donneesImm[-1])
data_moyAdo=colMeans(donneesAdo[-1])
data_moy=data_moyImm+data_moyAdo

# Visualisation de data_moy
plot(data_moy, type = "l", 
     xlab = "Largeur du crabe en mm",
     ylab = "Densité moyenne")
lines(data_moyImm, col="green2")
lines(data_moyAdo, col="yellow2")
grid()

#######################################################################
#
# Modélisation et fonctions
#
#######################################################################

intervalle = 1:95

# Paramètres initiaux de l'optimisation de la fonction de vraissemblance
theta=c(div_mu_imm=28, div_var_imm=3, epsilon=0.1, p=log(0.5268567/(1-0.5268567)), 
        w1=0.1,w2=0.1,w3=1,w4=1,w5=1,w6=0.001,w7=0.1, w8=0.01)


# Fonction de distribution d'un modèle de mélange Gaussien 
# (Gaussian Mixture model)
G = function(taille,w,mu,log_sigma){
  g=0
  for(i in 1:length(mu)){
    g=g+w[i]*(1/(sqrt(2*pi)*exp(log_sigma[i])))*exp((-1/(2*exp(2*log_sigma[i])))*((taille-mu[i])^2))
  }
  return(g)
}

# Calcul des moyennes immatures et adolescents
Calcul_mu = function(p,mu_init,div_mu_imm){
  # Calcul des mu généraux pour les 5 premiers instars
  mu=c(mu_init)
  for(i in 1:5){
    mu = c(mu, a*mu[i]+b)
  }
  
  mu_imm = c(mu[1:4], a*div_mu_imm+b, a*(a*div_mu_imm+b)+b)
  
  div_mu_ado = (mu[4]-p*div_mu_imm)/(1-p)
  
  mu_ado = c(a_ado*div_mu_ado+b_ado, 
             a_ado*mu_imm[5]+b_ado, a_ado*mu_imm[6]+b_ado)
  
  names(mu_imm) = paste0("imm",1:6)
  names(mu_ado) = paste0("ado",1:3)
  names(div_mu_ado) = "div_mu_ado"
  
  return(c(mu_imm, mu_ado, div_mu_ado))
}

# Calcul des variances immatures et adolescentes
Calcul_variance = function(p,log_sigma_init,mu,mu_imm,mu_ado,epsilon,div_var_imm){
  variance = exp(2*c(log_sigma_init))
  for(i in 1:5){
    variance = c(variance, a^2*variance[i]+epsilon)
  }
  var_imm = c(variance[1:3], a^2*div_var_imm+epsilon, a^2*(a^2*div_var_imm+epsilon)+epsilon, a^2*(a^2*(a^2*div_var_imm+epsilon)+epsilon)+epsilon)
  
  div_var_ado = (variance[4]-p*div_var_imm-p*(mu-mu_imm)^2-(1-p)*(mu-mu_ado)^2)/(1-p)
  
  var_ado = c(a^2*div_var_ado+epsilon, a^2*var_imm[5]+epsilon,
              a^2*var_imm[6]+epsilon)
  
  names(var_imm) = paste0("imm",1:6)
  names(var_ado) = paste0("ado",1:3)
  
  return(c(var_imm, var_ado, div_var_ado))
}

# Fonction de vraissemblance pour une distribution de 
# Gaussian mixture model
L = function(x, y=data_moy){
  
  div_mu_imm=x["div_mu_imm"]
  
  div_var_imm=x["div_var_imm"]
  
  epsilon=x["epsilon"]
  
  p=x[paste0("p")]
  
  p=(exp(p)/(1+exp(p)))
  
  mu = Calcul_mu(p,10,div_mu_imm)
  
  mu_imm = mu[paste0("imm",1:6)]
  mu_ado = mu[paste0("ado", 1:3)]
  div_mu_ado = mu[length(mu)]
  
  variance = Calcul_variance(p,-0.210797245,mu_imm[4],div_mu_imm,div_mu_ado,epsilon,div_var_imm)
  
  var_imm = variance[paste0("imm",1:6)]
  var_ado = variance[paste0("ado", 1:3)]
  
  log_sigma_imm = log(sqrt(var_imm))
  log_sigma_ado = log(sqrt(var_ado))
  
  w=x[paste0("w",1:8)]
  w_temp=(exp(w)/(1+sum(exp(w))))
  w=c(w_temp,1-sum(w_temp))
  
  
  # la log-vraissemblance ainsi que des termes pour encourager à se 
  # comporter d'une certaine façon. Les chiffres de la variance_immature
  # viennent des variances trouvées pour le modèle immature seulement.
  l=(-1*y*log(G(intervalle,w[1:6],mu_imm,log_sigma_imm)+G(intervalle,w[7:9],mu_ado, log_sigma_ado))
     -1*data_moyImm*log(G(intervalle,w[1:6]/sum(w[1:6]),mu_imm,log_sigma_imm))
     -1*data_moyAdo*log(G(intervalle,w[7:9]/sum(w[7:9]),mu_ado, log_sigma_ado))
     +(p-0.5268567)^2
     +(mu_imm[5]-37.5)^2
     +(mu_ado[1]-41.7)^2)
  
  return(sum(l))
  
}


# On fait l'optimisation. On veut trouver le maximum de vraissemblance,
# mais optim ici trouve le minimum, donc on a multiplié par -1 dans L.
optimisation=optim(theta,L,control = list(maxit=1e4, abstol=1e-8))


#######################################################################
#
# Récupération et visualisation des résultats
#
######################################################################

epsilon = optimisation$par["epsilon"]

div_mu_imm = optimisation$par["div_mu_imm"]
div_var_imm = optimisation$par["div_var_imm"]

p=optimisation$par["p"]
p=(exp(p)/(1+exp(p)))

mu=Calcul_mu(p,10,div_mu_imm)
mu_imm = mu[paste0("imm",1:6)]
mu_ado = mu[paste0("ado", 1:3)]
div_mu_ado = mu[length(mu)]

variance = Calcul_variance(p,-0.210797245,mu_imm[4],div_mu_imm,div_mu_ado,epsilon,div_var_imm)
var_imm = variance[paste0("imm",1:6)]
var_ado = variance[paste0("ado", 1:3)]

log_sigma_imm = log(sqrt(var_imm))
log_sigma_ado = log(sqrt(var_ado))

w = exp(optimisation$par[paste0("w",1:8)])/(1+sum(exp(optimisation$par[paste0("w",1:8)])))
w = c(w, 1-sum(w))

# Graphique de la distribution représentant les données
plot(data_moy/sum(data_moy), type = "l", 
     xlab = "Largeur du crabe en mm",
     ylab = "Densité moyenne")
lines(G(intervalle,w[1:6],mu_imm,log_sigma_imm)+G(intervalle,w[7:9],mu_ado, log_sigma_ado), type = "l",
      xlab = "Largeur en mm",
      ylab = "",
      lwd = 2,
      col="cyan")
abline(v=c(mu_imm, mu_ado),lty=2)

plot(data_moyAdo/sum(data_moyAdo), type = "l", 
     xlab = "Largeur du crabe en mm",
     ylab = "Densité moyenne", col="yellow2",
     ylim = c(0,max(c(0.06,max(G(intervalle,w[7:9]/sum(w[7:9]),mu_ado, log_sigma_ado)),max(G(intervalle,w[1:6]/sum(w[1:6]),mu_imm,log_sigma_imm))))))
lines(G(intervalle,w[7:9]/sum(w[7:9]),mu_ado, log_sigma_ado), type = "l",
      xlab = "Largeur en mm",
      ylab = "",
      lwd = 2,
      col="gold2")
abline(v=mu_ado,lty=2,col="yellow2")

lines(data_moyImm/sum(data_moyImm), type = "l", 
      xlab = "Largeur du crabe en mm",
      ylab = "Densité moyenne",col = "chartreuse")
lines(G(intervalle,w[1:6]/sum(w[1:6]),mu_imm,log_sigma_imm), type = "l",
      xlab = "Largeur en mm",
      ylab = "",
      lwd = 2,
      col="green2")
abline(v=mu_imm,lty=2,col="chartreuse")



# library(MASS)
# library(readxl)
# 
# rm(list=ls()) # effacer la mémoire d el'environnement de travail
# ls() # vérifier la liste des objets
# set.seed(2024) # Résultats reproductibles
# 
# #######################################################################
# #
# # Partie Régression (taille après en fonction de taille avant)
# #
# #######################################################################
# 
# # Tableau des moyennes de tailles avant (x) et après (y) la mue
# MeanCW = data.frame(x=c(10, 14.5, 20.5), y=c(14.5, 20.5, 28.4))
# MeanCWAdo = data.frame(x=c(29, 37.5, 49.1), y=c(41.7, 49.1, 58.1))
# 
# # Régression linéaire (y=ax+b) de la relation entre la taille avant et
# # après la mue
# Croissance = lm(y~x, data = MeanCW)
# CroissanceAdo = lm(y~sqrt(x),dat = MeanCWAdo)
# 
# # Graphique des points réelles avant et après mue
# plot(c(10,14.5,20.5,28.4,37.5,49.1), c(14.5,20.5,28.4,37.5,49.1,58.1), type="o", xlim=c(0,100), ylim=c(0,100), 
#      xlab = "Taille avant", ylab = "Taille après",col="orange")
# points(c(10,14.5,20.5,28.4,37.5), c(14.5,20.5,28.4,37.5,49.1), type="o", col="red")
# points(c(10,14.5,20.5,28.4), c(14.5,20.5,28.4,41.7), type="o", col="gold2")
# 
# # Pente de la régression
# a=Croissance$coefficients[2]
# a_ado=CroissanceAdo$coefficients[2]
# 
# # Ordonnée de la régression
# b=Croissance$coefficients[1]
# b_ado=CroissanceAdo$coefficients[1]
# 
# plot(MeanCWAdo$x,MeanCWAdo$y,"o",xlim = c(0,60), ylim = c(0,60))
# # Droite de régression
# axe_x = seq(0,100,length=100)
# lines(axe_x, predict(CroissanceAdo, data.frame(x=axe_x)), col="blue")
# grid(nx=NULL, ny=NULL)
# 
# #######################################################################
# #
# # Extraction et visualisation des données
# #
# #######################################################################
# 
# 
# # On extrait les données du excel
# donneesImm=read_excel("~/Crabe/female densities by size 1997-2023.xlsx")
# donneesAdo=read_excel("~/Crabe/female densities by size 1997-2023.xlsx", 
#                       sheet = "Adolescent")
# 
# # Calcul de la densité moyenne des crabes de chaque taille
# # à travers les années
# data_moyImm=colMeans(donneesImm[-1])
# data_moyAdo=colMeans(donneesAdo[-1])
# data_moy=data_moyImm+data_moyAdo
# 
# # Visualisation de data_moy
# plot(data_moy, type = "l", 
#      xlab = "Largeur du crabe en mm",
#      ylab = "Densité moyenne")
# lines(data_moyImm, col="green2")
# lines(data_moyAdo, col="yellow2")
# grid()
# 
# #######################################################################
# #
# # Modélisation et fonctions
# #
# #######################################################################
# 
# intervalle = 1:95
# 
# # Paramètres initiaux de l'optimisation de la fonction de vraissemblance
# theta=c(div_mu_imm=28, div_var_imm=3, epsilon=0.1, p=0.9999,#p=0.6744749, 
#         w1=0.1,w2=0.1,w3=1,w4=1,w5=1,w6=0.001,w7=0.1, w8=0.01)
# 
# 
# # Fonction de distribution d'un modèle de mélange Gaussien 
# # (Gaussian Mixture model)
# G = function(taille,w,mu,log_sigma){
#   g=0
#   for(i in 1:length(mu)){
#     g=g+w[i]*(1/(sqrt(2*pi)*exp(log_sigma[i])))*exp((-1/(2*exp(2*log_sigma[i])))*((taille-mu[i])^2))
#   }
#   return(g)
# }
# 
# # Calcul des moyennes immatures et adolescents
# Calcul_mu = function(p,mu_init,div_mu_imm){
#   # Calcul des mu généraux pour les 5 premiers instars
#   mu=c(mu_init)
#   for(i in 1:5){
#     mu = c(mu, a*mu[i]+b)
#   }
#   
#   mu_imm = c(mu[1:4], a*div_mu_imm+b, a*(a*div_mu_imm+b)+b-1)
#   
#   div_mu_ado = (mu[4]-p*div_mu_imm)/(1-p)
#   
#   mu_ado = c(a_ado*sqrt(div_mu_ado)+b_ado, 
#              a_ado*sqrt(mu_imm[5])+b_ado, a_ado*sqrt(mu_imm[6])+b_ado)
#   
#   names(mu_imm) = paste0("imm",1:6)
#   names(mu_ado) = paste0("ado",1:3)
#   names(div_mu_ado) = "div_mu_ado"
#   
#   return(c(mu_imm, mu_ado, div_mu_ado))
# }
# 
# # Calcul des variances immatures et adolescentes
# Calcul_variance = function(p,log_sigma_init,mu,mu_imm,mu_ado,epsilon,div_var_imm){
#   variance = exp(2*c(log_sigma_init))
#   for(i in 1:5){
#     variance = c(variance, a^2*variance[i]+epsilon)
#   }
#   var_imm = c(variance[1:3], a^2*div_var_imm+epsilon, a^2*(a^2*div_var_imm+epsilon)+epsilon, a^2*(a^2*(a^2*div_var_imm+epsilon)+epsilon)+epsilon)
#   
#   div_var_ado = (variance[4]-p*div_var_imm-p*(mu-mu_imm)^2-(1-p)*(mu-mu_ado)^2)/(1-p)
#   
#   var_ado = c(a^2*div_var_ado+epsilon, a^2*var_imm[5]+epsilon,
#               a^2*var_imm[6]+epsilon)
#   
#   names(var_imm) = paste0("imm",1:6)
#   names(var_ado) = paste0("ado",1:3)
#   
#   return(c(var_imm, var_ado))
# }
# 
# # Fonction de vraissemblance pour une distribution de 
# # Gaussian mixture model
# L = function(x, y=data_moy){
#   
#   div_mu_imm=x["div_mu_imm"]
#   
#   div_var_imm=x["div_var_imm"]
#   
#   epsilon=x["epsilon"]
#   
#   p=x[paste0("p")]
#   
#   p=(exp(p)/(1+exp(p)))
#   
#   mu = Calcul_mu(p,10,div_mu_imm)
#   
#   mu_imm = mu[paste0("imm",1:6)]
#   mu_ado = mu[paste0("ado", 1:3)]
#   div_mu_ado = mu[length(mu)]
#   
#   variance = Calcul_variance(p,-0.210797245,mu_imm[4],div_mu_imm,div_mu_ado,epsilon,div_var_imm)
#   
#   var_imm = variance[paste0("imm",1:6)]
#   var_ado = variance[paste0("ado", 1:3)]
#   
#   log_sigma_imm = log(sqrt(var_imm))
#   log_sigma_ado = log(sqrt(var_ado))
#   
#   w=x[paste0("w",1:8)]
#   w_temp=(exp(w)/(1+sum(exp(w))))
#   w=c(w_temp,1-sum(w_temp))
#   
#   
#   # la log-vraissemblance ainsi que des termes pour encourager à se 
#   # comporter d'une certaine façon. Les chiffres de la variance_immature
#   # viennent des variances trouvées pour le modèle immature seulement.
#   l=(-1*y*log(G(intervalle,w[1:6],mu_imm,log_sigma_imm)+G(intervalle,w[7:9],mu_ado, log_sigma_ado))
#      -1*data_moyImm*log(G(intervalle,w[1:6]/sum(w[1:6]),mu_imm,log_sigma_imm))
#      -1*data_moyAdo*log(G(intervalle,w[7:9]/sum(w[7:9]),mu_ado, log_sigma_ado)))
#   # +500*(mu_ado[1]-41.7)^2
#   # +500*(mu_imm[5]-37.5)^2
#   # +500*(var_imm[2]-1.67)^2
#   # +500*(var_imm[3]-3.85)^2
#   # +1*(var_imm[4]-7.43)^2
#   # +1*(var_imm[5]-16.13)^2
#   # +1*(var_imm[6]-44.6)^2)
#   
#   return(sum(l))
#   
# }
# 
# 
# # On fait l'optimisation. On veut trouver le maximum de vraissemblance,
# # mais optim ici trouve le minimum, donc on a multiplié par -1 dans L.
# optimisation=optim(theta,L,control = list(maxit=1e4, abstol=1e-8))
# 
# 
# #######################################################################
# #
# # Récupération et visualisation des résultats
# #
# ######################################################################
# 
# epsilon = optimisation$par["epsilon"]
# 
# div_mu_imm = optimisation$par["div_mu_imm"]
# div_var_imm = optimisation$par["div_var_imm"]
# 
# p=optimisation$par["p"]
# p=(exp(p)/(1+exp(p)))
# 
# mu=Calcul_mu(p,10,div_mu_imm)
# mu_imm = mu[paste0("imm",1:6)]
# mu_ado = mu[paste0("ado", 1:3)]
# div_mu_ado = mu[length(mu)]
# 
# variance = Calcul_variance(p,-0.210797245,mu_imm[4],div_mu_imm,div_mu_ado,epsilon,div_var_imm)
# var_imm = variance[paste0("imm",1:6)]
# var_ado = variance[paste0("ado", 1:3)]
# 
# log_sigma_imm = log(sqrt(var_imm))
# log_sigma_ado = log(sqrt(var_ado))
# 
# w = exp(optimisation$par[paste0("w",1:8)])/(1+sum(exp(optimisation$par[paste0("w",1:8)])))
# w = c(w, 1-sum(w))
# 
# # Graphique de la distribution représentant les données
# plot(data_moy/sum(data_moy), type = "l", 
#      xlab = "Largeur du crabe en mm",
#      ylab = "Densité moyenne")
# lines(G(intervalle,w[1:6],mu_imm,log_sigma_imm)+G(intervalle,w[7:9],mu_ado, log_sigma_ado), type = "l",
#       xlab = "Largeur en mm",
#       ylab = "",
#       lwd = 2,
#       col="cyan")
# abline(v=c(mu_imm, mu_ado),lty=2)
# 
# plot(data_moyAdo/sum(data_moyAdo), type = "l", 
#      xlab = "Largeur du crabe en mm",
#      ylab = "Densité moyenne", col="yellow2",
#      ylim = c(0,max(c(0.06,max(G(intervalle,w[7:9]/sum(w[7:9]),mu_ado, log_sigma_ado)),max(G(intervalle,w[1:6]/sum(w[1:6]),mu_imm,log_sigma_imm))))))
# lines(G(intervalle,w[7:9]/sum(w[7:9]),mu_ado, log_sigma_ado), type = "l",
#       xlab = "Largeur en mm",
#       ylab = "",
#       lwd = 2,
#       col="gold2")
# abline(v=mu_ado,lty=2,col="yellow2")
# 
# lines(data_moyImm/sum(data_moyImm), type = "l", 
#       xlab = "Largeur du crabe en mm",
#       ylab = "Densité moyenne",col = "chartreuse")
# lines(G(intervalle,w[1:6]/sum(w[1:6]),mu_imm,log_sigma_imm), type = "l",
#       xlab = "Largeur en mm",
#       ylab = "",
#       lwd = 2,
#       col="green2")
# abline(v=mu_imm,lty=2,col="chartreuse")
# 
# 
