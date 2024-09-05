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
MeanCW = data.frame(x=c(10, 14.5, 20.5, 28.4), y=c(14.5, 20.5, 28.4, 37.5))

# Régression linéaire (y=ax+b) de la relation entre la taille avant et
# après la mue
Croissance = lm(y~x, data = MeanCW)

# Graphique des points réelles avant et après mue
plot(MeanCW$x, MeanCW$y, type="o", xlim=c(0,40), ylim=c(0,40), 
     xlab = "Taille avant", ylab = "Taille après")

# Pente de la régression
a=Croissance$coefficients[2]
# Ordonnée de la régression
b=Croissance$coefficients[1]

# Droite de régression
axe_x = seq(0,40,length=100)
lines(axe_x, predict(Croissance, data.frame(x=axe_x)), col="blue")
grid(nx=NULL, ny=NULL)

#######################################################################
#
# Extraction et visualisation des données
#
#######################################################################


# On extrait les données du excel
donnees=read_excel("~/Crabe/female densities by size 1997-2023.xlsx")

# Calcul de la densité moyenne des crabes de chaque taille
# à travers les années
data_moy=colMeans(donnees[-1])

# Visualisation de data_moy
plot(1:length(data_moy), data_moy, type = "l", 
     xlab = "Largeur du crabe en mm",
     ylab = "Densité moyenne")
grid()

#######################################################################
#
# Modélisation et fonctions
#
#######################################################################

intervalle = 1:95

# Fonction de vraissemblance d'une distribution normale
L_simple = function(x, n=length(data_moy), f=data_moy, taille=1:length(data_moy)){
  mu=x["mu"]
  log_sigma=x["log_sigma"]
  
  -((-n/2)*log(2*pi)-n*log_sigma-(1/(2*exp(2*log_sigma)))*sum(f*(taille-mu)^2))
}

# Fonction de distribution d'un modèle de mélange Gaussien 
# (Gaussian Mixture model)
G = function(taille,w,mu,log_sigma){
  g=0
  for(i in 1:length(mu)){
    g=g+w[i]*(1/(sqrt(2*pi)*exp(log_sigma[i])))*exp((-1/(2*exp(2*log_sigma[i])))*((taille-mu[i])^2))
  }
  return(g)
}

# Fonction de vraissemblance pour une distribution de 
# Gaussian mixture model
L = function(x, cst_reg_moy, cst_reg_var, y=data_moy){
  mu=x[paste0("mu",1:6)]
  log_sigma=x[paste0("log_sigma", 1:6)]
  w=x[paste0("w",1:5)]
  epsilon=x["epsilon"]
  variance=(exp(log_sigma))^2
  
  k=cst_reg_moy
  c=cst_reg_var
  
  #top = y[c(10,15,21,28,37)]
  #top_proportion = top/sum(top)
  
  # Pour s'assurer que les poids sont positifs et que leur somme donne 1
  p=(exp(w)/(1+sum(exp(w))))
  p=c(p, 1-sum(p))
  
  # la log-vraissemblance ainsi que des termes de régularisation qui suit
  # ce qu'on sait (la croissance suit une droite ax+b donc mu2=a*mu1+b
  # et var2=a^2*var1)
  l=(-y*log(G(intervalle,p,mu,log_sigma))
     +k*(a*mu[1]+b-mu[2])^2 + c*(a^2*variance[1]+epsilon-variance[2])^2
     +k*(a*mu[2]+b-mu[3])^2 + c*(a^2*variance[2]+epsilon-variance[3])^2
     +k*(a*mu[3]+b-mu[4])^2 + c*(a^2*variance[3]+epsilon-variance[4])^2
     +k*(a*mu[4]+b-mu[5])^2 + c*(a^2*variance[4]+epsilon-variance[5])^2
     +k*(a*mu[5]+b-mu[6])^2 + c*(a^2*variance[5]+epsilon-variance[6])^2
     #+(p-top_proportion)^2
  )
  
  return(sum(l))
  
}

# theta utilisé pour la distribution normale simple
theta_simple=c(mu=10, log_sigma=0)

# Paramètres initiaux de l'optimisation de la fonction de vraissemblance
theta=c(mu1=10, log_sigma1=1, w1=1,
        mu2=14, log_sigma2=1, w2=1,
        mu3=20, log_sigma3=1, w3=1,
        mu4=28, log_sigma4=1, w4=1, 
        mu5=37, log_sigma5=1, w5=1,
        mu6=45, log_sigma6=1, epsilon=0)

# On fait l'optimisation. On veut trouver le maximum de vraissemblance,
# mais optim ici trouve le minimum, donc on a multiplié par -1 dans L.
optimisation=optim(theta,L,cst_reg_moy=3,cst_reg_var=0.001,control = list(maxit=1e4, abstol=1e-8))


#######################################################################
#
# Récupération et visualisation des résultats
#
######################################################################

mu=rep(0,6)
log_sigma=rep(0,6)
w=rep(0,6)

mu=optimisation$par[paste0("mu",1:6)]

log_sigma=optimisation$par[paste0("log_sigma", 1:6)]
sigma = exp(log_sigma)
variance = sigma^2

w[1:5]=exp(optimisation$par[paste0("w", 1:5)])/(1+sum(exp(optimisation$par[paste0("w", 1:5)])))
w[6]=1-sum(w[1:5])

# Graphique de la distribution représentant les données
plot(intervalle,G(intervalle,w,mu,log_sigma), type = "l",
     xlab = "Largeur en mm",
     ylab = "",
     ylim = c(0, 0.06))
abline(v=mu,lty=2, col="green")
# Calcul de l'intégrale pour s'assurer qu'elle est bien 1 
integrale=integrate(G,w=w,mu=mu,log_sigma=log_sigma,lower=1,upper=95)$value
text(x=80, y=0.047, labels = paste("Intégrale : ", integrale),cex=0.7)


lines(intervalle, data_moy[intervalle]/sum(data_moy[intervalle]), type="l", col="blue")



