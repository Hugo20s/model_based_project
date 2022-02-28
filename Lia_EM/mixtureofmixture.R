

# EM algorithm for Gaussian mixture model estimation
library("stats")
library("emdbook")
library("matrixStats")

etape_E <- function(data, K, parameters){ 
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  
  #Sigma -> a square variance-covariance matrix
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
  nk <- apply(tk,2,sum)
  pk <- nk/N
  #columns are features and the lines are the classes
  moyenne <- t(t(data)%*%tk/nk)
  data <- as.matrix(data)
  
  for (k in 1:K){
    temp <- replace(tk[, k], tk[, k]==0, 10^-8) 
    variance[,,k] = cov.wt(data, wt = temp , method = "ML")$cov
  }
  # for(k in 1:K) {
  #   temp <- 0
  #   for(n in 1:N) {
  #     temp <- temp + tk[n, k] * t(data[n,] - moyenne[k,]) * (data[n,] - moyenne[k,])}
  #   print(temp)
  #   #variance[,,k] <- temp/nk[k]
  # } 
  print("moyenne:")
  print(moyenne)
  print(" ")
  print("variance: ")
  print(variance)
  print(" ")
  print("tk: ")
  print(tk)
  print(" ")
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
      log_k[j] <- logSumExp(lapply(1:K, function(k) {
        print(log(pk[k]))
        print(dmvnorm(data[j,], moyenne[k,], variance[,,k], log = TRUE))
        log(pk[k]) + dmvnorm(data[j,], moyenne[k,], variance[,,k], log = TRUE)
      }))
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