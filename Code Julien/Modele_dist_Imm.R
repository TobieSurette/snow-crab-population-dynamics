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
MeanCW = data.frame(x=c(10, 14.5, 20.5, 28.4, 37.5), y=c(14.5, 20.5, 28.4, 37.5, 49.1))

# Régression linéaire (y=ax+b) de la relation entre la taille avant et
# après la mue
Croissance = lm(y~x, data = MeanCW)

# Graphique des points réelles avant et après mue
plot(MeanCW$x, MeanCW$y, type="o", xlim=c(0,50), ylim=c(0,50), 
     xlab = "Taille avant", ylab = "Taille après", lwd = 2)

# Pente de la régression
a=Croissance$coefficients[2]
# Ordonnée de la régression
b=Croissance$coefficients[1]

# Droite de régression
axe_x = seq(0,40,length=100)
lines(axe_x, predict(Croissance, data.frame(x=axe_x)), col="blue", lwd = 2)
grid(nx=NULL, ny=NULL)
legend(x = "bottomright", legend = c("Données réelles", "Régression"), fill = c("black", "blue"))

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
     ylab = "Densité moyenne",
     xlim = c(0,70), ylim = c(0,max(data_moy)),
     main = "Distribution des immatures",
     col = "green3", lwd = 2)
polygon(c(1:length(data_moy), rev(1:length(data_moy))), c(data_moy, rep(0, length(data_moy))), col = rgb(0,1,0,0.4), border = NA)
grid()
abline(v=c(10, 14.5, 20.5, 28.4, 37.5, 49.1), lty = 2, col = "green", lwd = 3)

#######################################################################
#
# Modélisation et fonctions
#
#######################################################################

intervalle = 1:95

# Meilleur configuration est nb_instars = 6, cst_reg_moy = 3,
# cst_reg_var = 0.0001, reg_proportion=FALSE
# Mettre cst_reg_var = 0.001 pour un fit un peu meilleur pour 
# les plus gros crabes

nb_instars = 6

cst_reg_moy=3

cst_reg_var=0.0001

reg_proportion=FALSE

# theta utilisé pour la distribution normale simple
theta_simple=c(mu=10, log_sigma=0)

# Paramètres initiaux de l'optimisation de la fonction de vraissemblance
theta=c(mu1=10, log_sigma1=1, w1=1,
        mu2=14, log_sigma2=1, w2=1,
        mu3=20, log_sigma3=1, w3=1,
        mu4=28, log_sigma4=1, w4=1, 
        mu5=37, log_sigma5=1, w5=1,
        mu6=45, log_sigma6=1, epsilon=0)

# Fonction de vraissemblance d'une distribution normale
L_simple = function(x, n=length(data_moy), y=data_moy){
  mu=x["mu"]
  log_sigma=x["log_sigma"]

  -((-n/2)*log(2*pi) - n*log_sigma - (1/(2*exp(2*log_sigma)))*sum((y-mu)^2))
  
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
L = function(x, k=cst_reg_moy, c=cst_reg_var, y=data_moy, n=nb_instars, tailles = 1:95 ){
  mu=x[paste0("mu",1:n)]
  log_sigma=x[paste0("log_sigma", 1:n)]
  w=x[paste0("w",1:(n-1))]
  epsilon=x["epsilon"]
  variance=(exp(log_sigma))^2
  
  top = y[floor(mu)]
  top_proportion = top/sum(top)
  
  # Pour s'assurer que les poids sont positifs et que leur somme donne 1
  p=(exp(w)/(1+sum(exp(w))))
  p=c(p, 1-sum(p))
  
  # la log-vraissemblance ainsi que des termes de régularisation qui suit
  # ce qu'on sait (la croissance suit une droite ax+b donc mu2=a*mu1+b
  # et var2=a^2*var1)
  l=-y*log(G(tailles,p,mu,log_sigma))
  for(i in 1:(n-1)){
    l=l+k*(a*mu[i]+b-mu[i+1])^2 + c*(a^2*variance[i]+epsilon-variance[i+1])^2
  }
  if(reg_proportion) l=l+(p-top_proportion)^2
  
  return(sum(l))
  
}

# On fait l'optimisation. On veut trouver le maximum de vraissemblance,
# mais optim ici trouve le minimum, donc on a multiplié par -1 dans L.
optimisation=optim(theta, L, control = list(maxit=1e4, abstol=1e-8))


#######################################################################
#
# Récupération et visualisation des résultats
#
######################################################################

mu=rep(0,nb_instars)
log_sigma=rep(0,nb_instars)
w=rep(0,nb_instars)

mu=optimisation$par[paste0("mu",1:nb_instars)]

log_sigma=optimisation$par[paste0("log_sigma", 1:nb_instars)]
sigma = exp(log_sigma)
variance = sigma^2

w[1:(nb_instars-1)]=exp(optimisation$par[paste0("w", 1:(nb_instars-1))])/(1+sum(exp(optimisation$par[paste0("w", 1:(nb_instars-1))])))
w[nb_instars]=1-sum(w[1:(nb_instars-1)])

# Graphique de la distribution représentant les données
plot(G(intervalle,w,mu,log_sigma), type = "l",
     xlab = "Largeur en mm",
     ylab = "",
     main = "Fonction de distribution des crabes des neiges immatures",
     ylim = c(0, 0.06), xlim=c(0,70),
     lwd = 2)
lines(w[1]*dnorm(intervalle, mean = mu[1], sd = sigma[1]), lty = 2, col="blue")
lines(w[2]*dnorm(intervalle, mean = mu[2], sd = sigma[2]), lty = 2, col="red")
lines(w[3]*dnorm(intervalle, mean = mu[3], sd = sigma[3]), lty = 2, col="green")
lines(w[4]*dnorm(intervalle, mean = mu[4], sd = sigma[4]), lty = 2, col="orange")
lines(w[5]*dnorm(intervalle, mean = mu[5], sd = sigma[5]), lty = 2, col="purple")
lines(w[6]*dnorm(intervalle, mean = mu[6], sd = sigma[6]), lty = 2, col="yellow3")
abline(v=mu,lty=2, col="grey")
# Calcul de l'intégrale pour s'assurer qu'elle est bien 1 
integrale=integrate(G,w=w,mu=mu,log_sigma=log_sigma,lower=1,upper=95)$value
#text(x=80, y=0.047, labels = paste("Intégrale : ", integrale),cex=0.7)


lines(intervalle, data_moy[intervalle]/sum(data_moy[intervalle]), type="l", col="cyan4",lwd = 2)
legend(x="topright", legend = c("Données réelles", "Données modélisées"), fill = c("cyan4", "black"))

mu_croissance=rep(0,length(mu))
mu_croissance[1]=theta["mu1"]
for(i in 1:(length(mu)-1)){
  mu_croissance[i+1] = a*mu_croissance[i]+b
}

plot(mu,type = "o", col = "blue")
lines(mu_croissance, type = "o", col = "red")
legend(x="bottomright", legend = c("mu fitté", "mu croissance (ax+b)"), 
       fill = c("blue","red"),cex = 0.7)
grid()


