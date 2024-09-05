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
#MeanCW = data.frame(x=c(10, 14.5, 20.5), y=c(14.5, 20.5, 28.4))
MeanCW = data.frame(x=c(10, 14.5, 20.5, 28.4), y=c(14.5, 20.5, 28.4, 37.5))

# Régression linéaire (y=ax+b) de la relation entre la taille avant et
# après la mue
Croissance = lm(y~x, data = MeanCW)

# Graphique des points réelles avant et après mue
plot(c(10,14.5,20.5,28.4,37.5,49.1), c(14.5,20.5,28.4,37.5,49.1,58.1), type="o", xlim=c(0,100), ylim=c(0,100), 
     xlab = "Taille avant", ylab = "Taille après",col="orange")
points(c(10,14.5,20.5,28.4,37.5), c(14.5,20.5,28.4,37.5,49.1), type="o", col="red")
points(c(10,14.5,20.5,28.4), c(14.5,20.5,28.4,41.7), type="o", col="gold2")

# Pente de la régression
a=Croissance$coefficients[2]
#a_imm=CroissanceImm$coefficients[2]

# Ordonnée de la régression
b=Croissance$coefficients[1]
#b_imm=CroissanceImm$coefficients[1]

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

#######################################################################
#
# Modélisation et fonctions
#
#######################################################################

intervalle = 1:95

# Paramètres initiaux de l'optimisation de la fonction de vraissemblance
theta=c(t=c(0.3,0.01), epsilon=0.1, p=c(log(0.63/(1-0.63)),0.001), 
        w1=0.1,w2=0.1,w3=1,w4=1,w5=1,w6=0.001,w7=0.1, w8=0.01)

#*******************************************
selectivite = function(taille){
  log(exp(1)/(1+exp(1/19*(10-taille))))
}
#*******************************************

# Fonction de distribution d'un modèle de mélange Gaussien 
# (Gaussian Mixture model)
G = function(taille,w,mu,log_sigma){
  g=0
  for(i in 1:length(mu)){
    g=g+w[i]*(1/(sqrt(2*pi)*exp(log_sigma[i])))*exp((-1/(2*exp(2*log_sigma[i])))*((taille-mu[i])^2))
  }
  return(g)
}

Instar_split = function(mu,variance,p,t){
  sigma=sqrt(variance)
  
  #mu = mu_a + t * qrt((1-p)/p) * sigma
  mu_a = mu - t * sqrt((1-p)/p) * sigma   # Mean size of first instar component (immature).
  mu_b = (mu - p * mu_a)/(1-p)            # Mean size of second instar component (adolescent).
  
  # Calculate instar variances, assuming that they are equal:
  V = sigma^2 - p * (1-p) * (mu_b-mu_a)^2
  sigma_a <- sigma_b <- sqrt(V)
  
  names(mu_a) = "mu_imm"
  names(mu_b) = "mu_ado"
  names(sigma_a) = "var_imm"
  names(sigma_b) = "var_ado"
  
  return(c(mu_a, mu_b, sigma_a^2, sigma_b^2))
}

# Calcul des moyennes et variances
Calcul_mu_var = function(mu_init,var_init,p,t,epsilon){
  
  #epsilon = exp(epsilon)/(1+exp(epsilon))
  
  mu=c(mu_init)
  variance=c(var_init)
  for(i in 1:5){
    mu = c(mu, a*mu[i]+b)
    variance = c(variance, a^2*variance[i]+epsilon)
  }
  
  split_VII = Instar_split(mu[4], variance[4], p[1], t[1])
  split_VIII = Instar_split(mu[5], variance[5], p[2], t[2])
  
  # On fait grandir le mu général et prends les valeurs du split sauf
  # pour l'instar 9 où on ne fait pas de split donc on fait grandir 
  # mu_imm précédent.
  # mu_imm = c(mu[1:3], split_VII["mu_imm"],
  #            split_VIII["mu_imm"], a*split_VIII["mu_imm"]+b)
  
  mu_imm = c(mu[1:4], a*split_VII["mu_imm"]+b, a*split_VIII["mu_imm"]+b)
  
  
  mu_ado = c(a*split_VII["mu_ado"]+b, a*split_VIII["mu_ado"]+b, 
             a*mu_imm[6]+b)
  
  
  var_imm = c(variance[1:3], split_VII["var_imm"], 
              split_VIII["var_imm"], a^2*split_VIII["var_imm"]+epsilon)

  
  var_ado = c(a^2*split_VII["var_ado"]+epsilon, 
              a^2*split_VIII["var_ado"]+epsilon,
              a^2*var_imm[6]+epsilon)
  
  names(var_imm) = paste0("var_imm",1:6)
  names(var_ado) = paste0("var_ado",1:3)
  
  names(mu_imm) = paste0("mu_imm",1:6)
  names(mu_ado) = paste0("mu_ado",1:3)
  
  return(c(mu_imm, mu_ado, var_imm, var_ado))
}


# Fonction de vraissemblance pour une distribution de 
# Gaussian mixture model
L = function(x, y=data_moy){
  
  epsilon=x["epsilon"]
  
  p=x[paste0("p",1:2)]
  
  p=(exp(p)/(1+sum(exp(p))))
  
  t=x[paste0("t",1:2)]
  
  t=(exp(t)/(1+sum(exp(t))))
  
  parametre = Calcul_mu_var(10,0.656,p,t,epsilon)
  
  mu_imm = parametre[paste0("mu_imm",1:6)]
  mu_ado = parametre[paste0("mu_ado", 1:3)]
  
  var_imm = parametre[paste0("var_imm",1:6)]
  var_ado = parametre[paste0("var_ado", 1:3)]
  
  log_sigma_imm = log(sqrt(var_imm))
  log_sigma_ado = log(sqrt(var_ado))
  
  w=x[paste0("w",1:8)]
  w=(exp(w)/(1+sum(exp(w))))
  
  w[5] = p[1]*w[4]
  w[6] = p[2]*w[5]

  w[7] = (1-p[1])*w[4]
  w[8] = (1-p[2])*w[5]
  
  w=c(w,1-sum(w))
  
  
  # la log-vraissemblance ainsi que des termes pour encourager à se 
  # comporter d'une certaine façon. Les chiffres de la variance_immature
  # viennent des variances trouvées pour le modèle immature seulement.
  l=(-1*y*log(G(intervalle,w[1:6],mu_imm,log_sigma_imm)+G(intervalle,w[7:9],mu_ado, log_sigma_ado))
     -1*data_moyImm*log(G(intervalle,w[1:6]/sum(w[1:6]),mu_imm,log_sigma_imm))
     -1*data_moyAdo*log(G(intervalle,w[7:9]/sum(w[7:9]),mu_ado, log_sigma_ado)))

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

p=optimisation$par[paste0("p",1:2)]
p=(exp(p)/(1+sum(exp(p))))

t=optimisation$par[paste0("t",1:2)]
t=(exp(t)/(1+sum(exp(t))))

parametre = Calcul_mu_var(10,0.656,p,t,epsilon)

mu_imm = parametre[paste0("mu_imm",1:6)]
mu_ado = parametre[paste0("mu_ado", 1:3)]

var_imm = parametre[paste0("var_imm",1:6)]
var_ado = parametre[paste0("var_ado", 1:3)]

log_sigma_imm = log(sqrt(var_imm))
log_sigma_ado = log(sqrt(var_ado))

w = exp(optimisation$par[paste0("w",1:8)])/(1+sum(exp(optimisation$par[paste0("w",1:8)])))
w[5] = p[1]*w[4]
w[6] = p[2]*w[5]
w[7] = (1-p[1])*w[4]
w[8] = (1-p[2])*w[5]
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
     ylab = "Densité moyenne", col="yellow4", lwd=2,
     ylim = c(0,max(c(0.06,max(G(intervalle,w[7:9]/sum(w[7:9]),mu_ado, log_sigma_ado)),max(G(intervalle,w[1:6]/sum(w[1:6]),mu_imm,log_sigma_imm))))))
lines(G(intervalle,w[7:9]/sum(w[7:9]),mu_ado, log_sigma_ado), type = "l",
      xlab = "Largeur en mm",
      ylab = "",
      lwd = 2,
      col="gold")
abline(v=mu_ado,lty=2,col="yellow2",lwd=3)

lines(data_moyImm/sum(data_moyImm), type = "l", 
      xlab = "Largeur du crabe en mm",
      ylab = "Densité moyenne",col = "chartreuse4", lwd=2)
lines(G(intervalle,w[1:6]/sum(w[1:6]),mu_imm,log_sigma_imm), type = "l",
      xlab = "Largeur en mm",
      ylab = "",
      lwd = 2,
      col="green")
abline(v=mu_imm,lty=2,col="chartreuse",lwd=3)

legend(x="topright", legend = c("Vrai immature", "Vrai adolescent", "Modèle immature", "Modèle adolescent"), 
       fill = c("chartreuse4", "yellow4", "green", "yellow"))

mtext(round(c(mu_imm[-6], mu_ado), digits = 2), side = 3, at = c(mu_imm[-6], mu_ado), cex = 0.7)


plot(data_moyAdo/sum(c(data_moyImm,data_moyAdo)), type = "l", 
     xlab = "Largeur du crabe en mm",
     ylab = "Densité moyenne", col="yellow4", lwd=2,
     ylim = c(0,max(c(0.06,max(G(intervalle,w[7:9],mu_ado, log_sigma_ado)),max(G(intervalle,w[1:6],mu_imm,log_sigma_imm))))))
lines(G(intervalle,w[7:9],mu_ado, log_sigma_ado), type = "l",
      xlab = "Largeur en mm",
      ylab = "",
      lwd = 2,
      col="gold")
abline(v=mu_ado,lty=2,col="yellow2",lwd=3)

lines(data_moyImm/sum(c(data_moyImm,data_moyAdo)), type = "l", 
      xlab = "Largeur du crabe en mm",
      ylab = "Densité moyenne",col = "chartreuse4", lwd=2)
lines(G(intervalle,w[1:6],mu_imm,log_sigma_imm), type = "l",
      xlab = "Largeur en mm",
      ylab = "",
      lwd = 2,
      col="green")
abline(v=mu_imm,lty=2,col="chartreuse",lwd=3)

legend(x="topright", legend = c("Vrai immature", "Vrai adolescent", "Modèle immature", "Modèle adolescent"), 
       fill = c("chartreuse4", "yellow4", "green", "yellow"))

mtext(round(c(mu_imm[-6], mu_ado), digits = 2), side = 3, at = c(mu_imm[-6], mu_ado), cex = 0.7)

