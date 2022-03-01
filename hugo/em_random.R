#Lia Furtado
#Model- based Learning

# EM algorithm for Gaussian mixture model estimation
library("stats")
library("emdbook")
etape_E <- function(data, K, parameters){ 
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  variance <- parameters$variance

  tk <- matrix(0,  nrow(data), K)
  for (k in 1:K){
    tk[, k] <- pk[k]*dmvnorm(as.matrix(data), moyenne[k,], variance[,,k])
  }
  tk <- tk/rowSums(tk)
  return (tk)
}
etape_M <- function(data, K, tk, parameters){
  
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  N <- nrow(data)
  for (k in 1:K){
    nk <- sum(tk[, k])
    pk[k] <- nk / nrow(data)
    moyenne[k,] <- (1/nk)*colSums(tk[, k] * data)
    
    #for (i in nrow(data)){
      #variance[,,k] <- variance[,,k] + tk[i, k] * t(as.matrix(data[i,] - moyenne[k,]))%*%as.matrix((data[i,] - moyenne[k,]))}
    #variance[,,k] <- 1/nk*variance[,,k] 
    
    variance[,,k] = cov.wt(data, wt = tk[, k], method = "ML")$cov
      }
  return (list(pk=pk, moyenne=moyenne, variance=variance))
  }


clustering <- function(data, K, epislon){
 
  #initialize variables
  pk <- rep(1/K, K)
  moyenne <- matrix(0, ncol=ncol(data), nrow = K)
  
  for (i in 1:ncol(data)){
    for (k in 1:K) {
      
      moyenne[k,i] <- sample(apply(data, 2, min)[i]:apply(data, 2, max)[i],1)
    }
  }
  variance <- array(0, dim=c(ncol(data),ncol(data),K))
  for (i in seq(K)){
    matrix_temp <- matrix(runif((ncol(data)^2)), (ncol(data)), ncol(data))
    matrix_temp[lower.tri(matrix_temp)] = t(matrix_temp)[lower.tri(matrix_temp)]
    matrix_temp <- t(matrix_temp) %*% matrix_temp
    variance[,,i] <- matrix_temp
  }
  #Loop 
  i <- 2
  log_likelihood <- numeric(0)
  log_likelihood[1] <- 0
  parameters <- list(pk=pk, moyenne=moyenne, variance=variance)
  while (i > 0){
    
    tk <- etape_E(data, K, parameters)
    parameters <- etape_M(data, K, tk, parameters)

    pk <- parameters$pk
    moyenne <- parameters$moyenne
    variance <- parameters$variance
    
    

    log_k <- numeric(nrow(data))
    for (j in 1:nrow(data)){
      log_k[j] <- logSumExp(lapply(1:K, function(k) log(pk[k]) + dmvnorm(as.matrix(data[j,]), moyenne[k,], variance[,,k], log = TRUE)))
    
    }
    print(log_k)
    log_likelihood[i] <- sum(log_k)
    if (abs(log_likelihood[i] - log_likelihood[i-1]) < epislon){
        break
    }

    i <- i + 1
  }
  y_pred <- apply(tk, 1, which.max)
  return(list(log_likelihood = log_likelihood[-1], y_pred = y_pred))
  
}
K <- 3
data <- iris[-sample(1:150, 50), -5]
epislon <- 10^-6
res <- clustering(data, K, epislon)
plot(res$log_likelihood)
plot(res$y_pred)
tkres<- replicate(5, clustering(data, K, epislon)) 

  
