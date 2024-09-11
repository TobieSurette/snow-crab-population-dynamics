library(MASS)

set.seed(2024)

Instars = matrix(,nrow=nb_instars,ncol=95)

for (i in 1:nb_instars){
  Instars[i,] = w[i]*dnorm(1:95,mean = mu[i], sd = sigma[i])/G(1:95, w, mu, log_sigma)
}

plot(1:95, Instars[1,], type = "l", col="blue",
     xlab = "Largeur en mm",
     ylab = "Probabilité d'appartenir à l'instar")
lines(1:95, Instars[2,], type = "l", col="red")
lines(1:95, Instars[3,], type = "l", col="yellow")
lines(1:95, Instars[4,], type = "l", col="green")
lines(1:95, Instars[5,], type = "l", col="purple")
lines(1:95, Instars[6,], type = "l", col="black")
grid()
legend(x="right",legend=c("IV","V","VI","VII","VIII","IX"), 
       fill=c("blue","red","yellow","green","purple","black"), cex=0.7)


# Fonction qui affiche le % de chance d'appartenir dans les instars
# pour une taille données
Instar_predictor = function(taille){
  Instars = matrix(,nrow=5,ncol=95)
  
  for (i in 1:5){
    Instars[i,] = w[i]*dnorm(1:95,mean = mu[i], sd = sigma[i])/G(1:95, w, mu, log_sigma)
  }
  
  print(paste0(paste("Instar IV:", 100*Instars[1,taille]),"%"))
  print(paste0(paste("Instar V:", 100*Instars[2,taille]),"%"))
  print(paste0(paste("Instar VI:", 100*Instars[3,taille]),"%"))
  print(paste0(paste("Instar VII:", 100*Instars[4,taille]),"%"))
  print(paste0(paste("Instar VIII:", 100*Instars[5,taille]),"%"))
  }


