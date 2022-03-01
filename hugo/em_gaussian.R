##EM algorithme 

##needed package 
install.packages("DirichletReg")

#load package 
library(DirichletReg)
library(mvtnorm)



lse <- function(x){
  c <- max(x)
  return(c + (log(sum(exp(x - c)))))
}


#Algorithme EM 
etape_E <- function(data, K, parameters){
  #step E of EM algorithm
  prop <- parameters$prop
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  
  n<- nrow(data); p <- ncol(data)
  
  tk <- matrix(0, n, K)
  temp <- matrix(0, n, K)
  for (k in 1:K){
    temp[, k] <- prop[k]*dmvnorm(data, moyenne[,k], variance[,,k])
  }
  tk <- temp/rowSums(temp)
  
  return (tk)
}


etape_M <- function(data, K, tk, parameters){
  #step M of EM algorithm
  prop <- parameters$prop
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  n<- nrow(data); p <- ncol(data)
  
  for (k in 1:K){
    prop[k] <- sum(tk[, k])/ n
    moyenne[,k] <- (1/sum(tk[, k]))*colSums( apply(data,2, function(x) x*tk[, k])  )
    variance[,,k] <- cov.wt((data),  tk[,k])$cov
    #a <- (1/sum(tk[, k]))
    #b <- t((tk[, k])*(data - moyenne[,k]))%*% (data - moyenne[,k])
    #variance[,,k] <- a*b
    
    
  }
  return (list(prop=prop, moyenne=moyenne, variance=variance))
}





EM <- function(data, K, epislon, init = "random"){
  
  n <- nrow(data); p <- ncol(data)
  #initialize variables
  if (init == "kmeans") {
    res_kmeans <- kmeans(data, K)
    prop <- numeric(K)
    for (k in 1:K){
      prop[k] <- length(which(res_kmeans$cluster==k))/n
    }
    moyenne <- t(as.matrix(res_kmeans$centers))
    
  } else if (init == 'random'){
    prop <- rdirichlet(1, rep(1:K)) #K
    moyenne <- t(data[sample(nrow(data), size = K),]) #p*K
    variance <- array(0, dim=c(p,p,K)) 
    var_temp <- var(data)
    variance <- rWishart(K, 4, var_temp)
    # for (k in seq(K)){
    #   var_temp <- apply(var_temp, 1, function(x) x * sample(c(1.1,0.9), 1))
    #   variance[,,k] <- var_temp
    # }
  }
 
  
  #Loop 
  i <- 2
  log_likelihood <- numeric(0)
  log_likelihood[1] <- 0
  
  parameters <- list(prop=prop, moyenne=moyenne, variance=variance)
  while (i > 0){
    tk <- etape_E(data, K, parameters)
    parameters <- etape_M(data, K, tk, parameters)
    
    
    prop <- parameters$prop
    moyenne <- parameters$moyenne
    variance <- parameters$variance
    
    #print
    
    log_likelihood[i] <- 0
    log_k <- matrix(0,n,K)
    for (k in 1:K){
      log_k[,k] <- prop[, k]*dmvnorm(data, moyenne[,k], variance[,,k], log = TRUE)
    }
    test <- colSums(log_k)
    log_n <- apply(log_k, 1, lse)
    
    log_likelihood[i] <- sum(log_n)
    
    
    #if (abs(log_likelihood[i] - log_likelihood[i-1]) < epislon){
      #break
    #}
    if (i == 30){
     break
    }
    i <- i + 1
  }
  plot(log_likelihood[-1])
  return(log_likelihood)
}




K <- 3

clus_1 <- matrix(rnorm(20*4, 0, 1), 20, 4)
clus_2 <- matrix(rnorm(20*4, 5, 2), 20, 4)
clus_3 <- matrix(rnorm(20*4, 10, 2), 20, 4)

data <- rbind(clus_1, clus_2, clus_3)
epislon <- 10^-6

undebug(etape_E)
undebug(etape_M)
undebug(EM)
log_likelihood <- EM(data, K, epislon)

