#Lia Furtado
#Model- based Learning

# EM algorithm for Gaussian mixture model estimation

etape_CE <- function(data, K, parameters){ 
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  variance <- parameters$variance

  tk <- matrix(0,  nrow(data), K)
  for (k in 1:K){
    tk[, k] <- pk[k]*dmvnorm(as.matrix(data), moyenne[k,], variance[,,k])
  }
  tk <- tk/rowSums(tk)
  
  #rouding the TK
  max_values <- apply(tk, 1, which.max)  
  j <- 1
  for (i in max_values){
    
    tk[j,] <- 0
    tk[j, i] <- 1
    j <- j + 1
  }
  print(tk)
  return (tk)
}
etape_M <- function(data, K, tk, parameters){
  
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  N <- nrow(data)
  for (k in 1:K){
    nk <- sum(tk[, k])
    moyenne[k,] <- (1/nk)*colSums(tk[, k] * data)
    
  }

  return (list(pk=pk, moyenne=moyenne, variance=variance))
  }


kmeans <- function(data, K, number_iterations){
 
  #initialize variables
  #same proportions
  pk <- rep(1/K, K)
  moyenne <- matrix(0, ncol=ncol(data), nrow = K)
  
  for (i in 1:ncol(data)){
    for (k in 1:K) {
      
      moyenne[k,i] <- sample(apply(data, 2, min)[i]:apply(data, 2, max)[i],1)
    }
  }
  #same covariance matrix
  variance <- array(0, dim=c(ncol(data),ncol(data),K))
  matrix_temp <- matrix(runif((ncol(data)^2)), (ncol(data)), ncol(data))
  matrix_temp[lower.tri(matrix_temp)] = t(matrix_temp)[lower.tri(matrix_temp)]
  matrix_temp <- t(matrix_temp) %*% matrix_temp

  variance[,,] <- matrix_temp
  #Loop 
  i <- 2
  log_likelihood <- numeric(0)
  log_likelihood[1] <- 0
  parameters <- list(pk=pk, moyenne=moyenne, variance=variance)
  while (i > 0){
    
    tk <- etape_CE(data, K, parameters)
    #print(tk)
    parameters <- etape_M(data, K, tk, parameters)
    #print(parameters)
  
    pk <- parameters$pk
    moyenne <- parameters$moyenne
    variance <- parameters$variance
    
    
    log_k <- numeric(K)
    for (k in 1:K){
      
      log_k[k] <- sum(tk[, k]*dmvnorm(as.matrix(data), moyenne[k,], variance[,,k], log = TRUE))

    
    }
    log_likelihood[i] <- sum(log_k)

    if (i == number_iterations){
      break
    }
    i <- i + 1
  }
  y_pred <- apply(tk, 1, which.max)
  return(list(log_likelihood = log_likelihood[-1], y_pred = y_pred))
  
}
K <- 3
data <- iris[-sample(1:150, 50), -5]
number_iterations <- 50
log_likelihood <- kmeans(data, K, number_iterations)
plot(log_likelihood$log_likelihood)
log_likelihood$y_pred
tkres<- replicate(5, clustering(data, K, number_iterations)) 

BIC <- numeric(0)
i <- 1
number_clusters <- c(2,3,4,5)
for (k in number_clusters){
  log_likelihood <- clustering(data, k, epislon)
  BIC[i] <- -2*max(log_likelihood$log_likelihood[-1]) + k*log(N)
  i <- i + 1
}
plot(number_clusters, BIC)

  
