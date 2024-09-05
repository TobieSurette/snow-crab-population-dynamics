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
L = function(x, y=image_grise){
  
  intervalle = seq(0,1,len=length(y))
  
  w=x["w"]
  
  w=c(w,(exp(w)/(1+sum(exp(w)))))
  
  mu=x[paste0("mu",1:2)]
  
  log_sigma = x[paste0("ls", 1:2)]
  
  l = -c(y)*log(G(intervalle,w,mu,log_sigma))
  
  return(sum(l))
  
}

TrouverContour2 = function(image){
  im = as.vector(image)
  
  div = matrix(kmeans(im,3)$cluster, nrow(image), ncol(image))
  
  centre = matrix(c(max(A1[,1])-min(A1[,1]), max(A1[,2])-min(A1[,2]),
                    max(A2[,1])-min(A2[,1]), max(A2[,2])-min(A2[,2]),
                    max(A3[,1])-min(A3[,1]), max(A3[,2])-min(A3[,2])),
                  3,2)
  
  
  x=c()
  y=c()
  angles = seq(0,2*pi, length.out = 1000)
  for(theta in angles){
    r=0
    k=100
    while(TRUE){
      if(div[round(centre[1]+(r+k)*cos(theta)),round(centre[2]+(r+k)*sin(theta))] == div[round(centre[1]),round(centre[2])]){
        r=r+k
      } else{
        k=k*0.1
      }
      if(k<1e-8) break
    }
    x=c(x,round(centre[1]+r*cos(theta)))
    y=c(y,round(centre[2]+r*sin(theta)))
  }
  
  Contour=cbind(x,y)
  
  return(Contour)
}

TrouverContour = function(image){
  centre = c(nrow(image)/2, ncol(image)/2)
  im = as.vector(image)
  
  div = matrix(kmeans(im,2)$cluster, nrow(image), ncol(image))
  
  A1 = which(div == 1, arr.ind = TRUE)
  A2 = which(div == 2, arr.ind = TRUE)
  
  if(mean(image[A1]) > mean(image[A2])){
    O = A1
  } else{
    O = A2
  }
  
  
  x=c()
  y=c()
  angles = seq(0,2*pi, length.out = 1000)
  for(theta in angles){
    r=0
    k=100
    while(TRUE){
      if(div[round(centre[1]+(r+k)*cos(theta)),round(centre[2]+(r+k)*sin(theta))] == div[round(centre[1]),round(centre[2])]){
        r=r+k
      } else{
        k=k*0.1
      }
      if(k<1e-8) break
    }
    x=c(x,round(centre[1]+r*cos(theta)))
    y=c(y,round(centre[2]+r*sin(theta)))
  }
  
  Contour=cbind(x,y)
  
  # for(j in min(O[,2]):max(O[,2])){
  #   for(i in min(O[,1]):max(O[,1]))
  #   if(div[i-1,j] != div[i,j] || div[i+1,j] != div[i,j]){
  #     if(is.na(Contour[1,1])){
  #       Contour[1,] = c(i, j)
  #     } else{
  #       Contour = rbind(Contour, c(i,j))
  #     }
  #   }
  # }
  
  # init = c(mu1 = 0.15, mu2 = 0.6, ls1 = -10, ls2 = -10, w=0.5)
  # 
  # param = optim(init,L)$par
  
  return(Contour)
}