library(MASS)

rm(list=ls()) # effacer la mémoire d el'environnement de travail
ls() # vérifier la liste des objets
set.seed(2024) # Résultats reproductibles

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

# Nombre total de crabes
n=1000

CW <- seq(from=0,to=70,length.out=1000)

# Fonction de la variance et de la moyenne qui dépend de la variance
# et moyenne de l'instar précédent
var_func = function(x){return(a*a*x)} #V(ax+b) = a^2*V(x)
moy_func = function(x){return(a*x+b)} #E(ax+b) = a*E(x)+b

# Variance et moyenne du premier instar
var_init=2/3
moy_init=10
nb_init=15

# On initialise et ensuite calcul les variances et moyennes 
# des autres instances
moyenne=rep(0,5)
moyenne[1]=moy_init
variance=rep(0,5)
variance[1]=var_init
for(i in 1:4){
variance[i+1]=var_func(variance[i])
moyenne[i+1]=moy_func(moyenne[i])
}

total = 264480.095833587 + 1960326.20709249 + 11343632.3036489 
+ 74844309.8641673 + 97584033.168346
#Dist_instars = dnorm(CW, mean=moyenne, sd=variance)

# Fonction logistique de sélectivité
select_log = function(x){
  1/(1+exp(-(1/6)*(x-25))) #milieu à 25mm
}

curve(select_log(x), from = 0, to=70, xlab = "x", ylab = "y")
title(main = "Courbe logistique de la sélectivité du filet")

# On crée les distributions des instars

#Imm=matrix(,nrow=n,ncol = 5)
#Imm[,1]=nb_init*dnorm(CW[,1], mean=moyenne[1], sd=variance[1])


IV = select_log(CW)*dnorm(CW, mean=moyenne[1], sd=sqrt(variance[1]))
V = select_log(CW)*dnorm(CW, mean = moyenne[2], sd=sqrt(variance[2]))
VI = select_log(CW)*dnorm(CW, mean = moyenne[3], sd=sqrt(variance[3]))
VII = select_log(CW)*dnorm(CW, mean = moyenne[4], sd=sqrt(variance[4]))
VIII = select_log(CW)*dnorm(CW, mean = moyenne[5], sd=sqrt(variance[5]))

# Graphique de la distribution des crabes immatures.
plot(CW, IV+V+VI+VII+VIII, type = "l", col = "green", 
   xlab = "Carapace width (mm)", ylab="", ylim=c(0,0.5))

abline(v=moyenne[1],lty=2, col="green")
abline(v=moyenne[2],lty=2, col="green")
abline(v=moyenne[3],lty=2, col="green")
abline(v=moyenne[4],lty=2, col="green")
abline(v=moyenne[5],lty=2, col="green")

################################################
#                                              #
#   Simulation des données                     #
#                                              #
################################################

Imm = matrix(,nrow=n, ncol=5)
Imm1=Imm

# Première méthode en utilisant 5 distributions normales 
for(i in 1:5){
  Imm[,i]=rnorm(n, mean=moyenne[i], sd=sqrt(variance[i]))
}

# Deuxième méthode en utilisant la fonction de croissance
Imm1[,1]=rnorm(n, mean=moyenne[1], sd=sqrt(variance[1]))
for(i in 2:5){
  Imm1[,i]=a*Imm1[,i-1]+b*rep(1,n)
}

hist(Imm, 
     main = "Nombre de crabes pour chaque taille sur 5000 crabes",
     xlab = "Largeur de carapace en mm",
     ylab = "# de crabes",
     breaks = 200)

hist(Imm1, 
     main = "Nombre de crabes pour chaque taille sur 5000 crabes",
     xlab = "Largeur de carapace en mm",
     ylab = "# de crabes",
     breaks = 200)

P=c(0.6, 0.6, 0.6, 0.6, 0.6)
G=rep(0.3, 5)

M = matrix(c(P[1], G[1], 0, 0, 0, #(colonnes x rangées)
             0, P[2], G[2], 0, 0,
             0, 0, P[3], G[3], 0,
             0, 0, 0, P[4], G[4],
             0, 0, 0, 0, P[5]),
           nrow = 5, 
           ncol = 5)

N_init = c(1000,1000,1000,1000,1000)

N = matrix(,nrow = 5, ncol = 100)

N[,1] = N_init
for(t in 1:99){
  N[,t+1] = M%*%N[,t]+c(1000,0,0,0,0)
}

plot(1:100,N[1,], type = "l", col = "blue", 
     xlab = "temps", ylab="# de crabes", xlim = c(0,100), ylim = c(0,10000))
lines(1:100,N[2,],col = "cyan")
lines(1:100,N[3,],col = "orange")
lines(1:100,N[4,], col = "green")
lines(1:100,N[5,], col = "red")

legend(x="topright", legend=c("IV","V","VI","VII","VIII"), 
       fill=c("blue","cyan","orange","green","red"), cex=0.7)


