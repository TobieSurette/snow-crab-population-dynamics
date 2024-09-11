library(MASS)
library(readxl)

# rm(list=ls()) # effacer la mémoire d el'environnement de travail
# ls() # vérifier la liste des objets
set.seed(2024) # Résultats reproductibles

#######################################################################
#
# Partie Régression (taille après en fonction de taille avant)
#
#######################################################################

# Tableau des moyennes de tailles avant (x) et après (y) la mue
MeanCW = data.frame(x=c(10, 14.5, 20.5), y=c(14.5, 20.5, 28.4))
MeanCWImm = data.frame(x=c(10, 14.5, 20.5, 28.4), y=c(14.5, 20.5, 28.4, 37.5))

# Régression linéaire (y=ax+b) de la relation entre la taille avant et
# après la mue
Croissance = lm(y~x, data = MeanCW)
CroissanceImm = lm(y~x, data = MeanCWImm)

# Graphique des points réelles avant et après mue
plot(c(10,14.5,20.5,28.4,37.5,49.1), c(14.5,20.5,28.4,37.5,49.1,58.1), type="o", xlim=c(0,100), ylim=c(0,100), 
     xlab = "Taille avant", ylab = "Taille après",col="orange")
points(c(10,14.5,20.5,28.4,37.5), c(14.5,20.5,28.4,37.5,49.1), type="o", col="red")
points(c(10,14.5,20.5,28.4), c(14.5,20.5,28.4,41.7), type="o", col="gold2")

# Pente de la régression
a=Croissance$coefficients[2]
a_imm=CroissanceImm$coefficients[2]

# Ordonnée de la régression
b=Croissance$coefficients[1]
b_imm=CroissanceImm$coefficients[1]

# Droite de régression
axe_x = seq(0,100,length=100)
lines(axe_x, predict(Croissance, data.frame(x=axe_x)), col="blue")
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

plot(data_moyImm, col="green3", lwd = 2, type = "l",main = "Distribution des immatures et adolescents", xlab = "Largeur du crabe en mm", ylab = "Densité moyenne (#/km2)")
polygon(c(1:length(data_moyImm), rev(1:length(data_moyImm))), c(data_moyImm, rep(0, length(data_moyImm))), col = rgb(0,1,0,0.4), border = NA)
lines(data_moyAdo, col="yellow3", lwd = 2)
polygon(c(1:length(data_moyAdo), rev(1:length(data_moyAdo))), c(data_moyAdo, rep(0, length(data_moyAdo))), col = rgb(249/256,214/256,148/256,0.4), border = NA)
abline(v=c(10, 14.5, 20.5, 28.4, 37.5, 49.1), lty = 2, col = "green", lwd = 3)
abline(v = c(41.7, 49.1, 58.1), lty = 2, col = "yellow", lwd = 3)
mtext(c("IV", "V", "VI", "VII", "VIII", "VIII", "IX"), side = 3, line = 0.5, at = c(10, 14.5, 20.5, 28.4, 37.5, 41.7, 49.1, 58.1))

#######################################################################
#
# Modélisation et fonctions
#
#######################################################################

intervalle = 1:95

# Paramètres initiaux de l'optimisation de la fonction de vraissemblance
theta=c(epsilon=0.1, p=0.5, 
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
Calcul_mu = function(p,mu_init){
  # Calcul des mu généraux pour les 5 premiers instars
  mu=c(mu_init)
  for(i in 1:5){
    mu = c(mu, a*mu[i]+b)
  }
  # Les 4 premiers instars ont les mêmes moyennes que général
  # On utilise les coefficients de la régression immature pour les autres
  mu_imm = c(mu[1:4], a_imm*mu[4]+b_imm, a_imm*(a_imm*mu[4]+b_imm)+b_imm)
  
  #Pas oublié que j'ai ajouté un -5.5
  # On utilise le fait que mu=p*mu_imm+(1-p)*mu_ado pour trouver le mu
  # de l'instar 8. On utilise la croissance immature pour les deux autres
  # instars
  mu_ado = c((mu[5]-p*mu_imm[5])/(1-p), mu_imm[6], a_imm*mu_imm[6]+b_imm)
  
  names(mu_imm) = paste0("imm",1:6)
  names(mu_ado) = paste0("ado",1:3)
  
  return(c(mu_imm, mu_ado))
}

# Calcul des variances immatures et adolescentes
Calcul_variance = function(p,log_sigma_init,mu_imm,mu_ado,epsilon){
  variance = exp(2*c(log_sigma_init))
  for(i in 1:5){
    variance = c(variance, a^2*variance[i]+epsilon)
  }
  var_imm = c(variance[1:4], a_imm^2*variance[4]+epsilon, a_imm^2*(a_imm^2*variance[4]+epsilon)+epsilon)

  var_ado = c((variance[5]-p*var_imm[5]-p*(1-p)*(mu_imm[5]-mu_ado[1])^2)/(1-p),
              var_imm[6], a_imm^2*var_imm[6]+epsilon)
  
  names(var_imm) = paste0("imm",1:6)
  names(var_ado) = paste0("ado",1:3)
  
  return(c(var_imm, var_ado))
}

# Fonction de vraissemblance pour une distribution de 
# Gaussian mixture model
L = function(x, y=data_moy){
  
  epsilon=x["epsilon"]
  
  p=x[paste0("p")]
  
  p=(exp(p)/(1+exp(p)))
  
  mu = Calcul_mu(p,10)
  
  mu_imm = mu[paste0("imm",1:6)]
  mu_ado = mu[paste0("ado", 1:3)]
  
  variance = Calcul_variance(p,-0.210797245,mu_imm,mu_ado,epsilon)
  
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
     -50*data_moyImm*log(G(intervalle,w[1:6]/sum(w[1:6]),mu_imm,log_sigma_imm))
     -50*data_moyAdo*log(G(intervalle,w[7:9]/sum(w[7:9]),mu_ado, log_sigma_ado))
     +10*(mu_ado[1]-41.7)^2
     +500*(var_imm[2]-1.67)^2
     +500*(var_imm[3]-3.85)^2
     +1*(var_imm[4]-7.43)^2
     +1*(var_imm[5]-16.13)^2
     +1*(var_imm[6]-44.6)^2)
  
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

p=optimisation$par["p"]
p=(exp(p)/(1+exp(p)))

mu=Calcul_mu(p,10)
mu_imm = mu[paste0("imm",1:6)]
mu_ado = mu[paste0("ado", 1:3)]

variance = Calcul_variance(p,-0.210797245,mu_imm,mu_ado,epsilon)
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


