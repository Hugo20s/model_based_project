# EM algorithm for Gaussian mixture model estimation
library("stats")
library("emdbook")
library("matrixStats")
mysum <- function(x) {
  sum(x[is.finite(x)])
}
init_random <- function(data, K){
  pk <- rep(1/K, K)
  
  moyenne <- matrix(0, ncol=ncol(data), nrow = K)
  
  for (k in 1:K) {
    rand_1 <- runif(1, min(data[,1]), max(data[,1]))
    rand_2 <- runif(1, min(data[,2]), max(data[,3]))
    rand_3 <-   runif(1, min(data[,3]), max(data[,3]))
    rand_4 <- runif(1, min(data[,4]), max(data[,4]))
    moyenne[k,] <- c(rand_1, rand_2, rand_3, rand_4)
  }
  
  covariance <- array(0, dim=c(ncol(data),ncol(data),K))
  for (i in seq(K)){
    matrix_temp <- matrix(runif((ncol(data)^2)), (ncol(data)), ncol(data))
    matrix_temp[lower.tri(matrix_temp)] = t(matrix_temp)[lower.tri(matrix_temp)]
    matrix_temp <- t(matrix_temp) %*% matrix_temp
    covariance[,,i] <- matrix_temp
  }
  return (list(pk=pk, moyenne=moyenne, variance=covariance))
}

etape_E <- function(data, K, parameters){ 
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  
  
  tk <- matrix(0,  nrow(data), K)
  for (k in 1:K){
    # print("density")
    # print(k)
    # print(pk[k])
    # print(dmvnorm(data, moyenne[k,], variance[,,k]))
    tk[, k] <- pk[k]*dmvnorm(data, moyenne[k,], variance[,,k])
  }
  tk <- tk/apply(tk, 1, mysum)
  
  return (tk)
}

etape_M <- function(data, K, tk, parameters){
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  N <- nrow(data)
  f <- ncol(data)
  nk <- colSums(tk)
  pk <- nk/N
  moyenne <- t(t(data) %*% tk) / nk #dimension(k,f)
  covariance_aux <- lapply(1:K, function(j) matrix(
    apply(
      sapply(1:N, function(i) {
        temp <- replace(tk[i, j], tk[i, j]==0, 10^-8)
        temp * (data[i, ] - moyenne[[j]]) %*%
          t(data[i, ] - moyenne[[j]])}
      ), 1, mysum
    ), f, f)/mysum(tk[,j]))
  variance <- array(unlist(covariance_aux), dim=c(f,f,K))
  
  return (list(pk=pk, moyenne=moyenne, variance=variance))
}


clustering <- function(data, K, epislon){
  set_inf <- FALSE
  #--------initialize variables
  n <- nrow(data); p <- ncol(data)
  if (FALSE){
    pk <- c(0.2, 0.5, 0.3)
    moyenne <- matrix(0, ncol=ncol(data), nrow = K)
    variance <- array(0, dim=c(ncol(data),ncol(data),K))
    for (k in 1:K) {
      rand_1 <- runif(1, min(data[,1]), max(data[,1]))
      rand_2 <- runif(1, min(data[,2]), max(data[,3]))
      rand_3 <-   runif(1, min(data[,3]), max(data[,3]))
      rand_4 <- runif(1, min(data[,4]), max(data[,4]))
      moyenne[k,] <- c(rand_1, rand_2, rand_3, rand_4)
      
      n_sample <- 0.7*n
      print(data[sample(1:n, n_sample),])
      variance[,,k] <- cov(data[sample(1:n, n_sample),])
    }
    parameters <- list(pk=pk, moyenne=moyenne, variance=variance)
  }
  
  parameters <- init_random(data, K)
  # #initialize variables
  # km.res <- kmeans(data, K)
  # pk <- rep(1/K, K)
  # moyenne <- as.matrix(km.res$centers)
  # 
  # variance <- array(0, dim=c(ncol(data),ncol(data),K))
  # for (i in seq(K)){
  #   variance[,,i] <- cov(data[which(as.matrix(km.res$cluster) == i),])
  # }

  #Loop 
  i <- 2
  log_likelihood <- numeric(0)
  log_likelihood[1] <- -Inf
  
  
  while (i > 0) {
    tk <- etape_E(data, K, parameters)
    
    for (k in 1:K ){
      if (sum(tk[,k]) == 0){
        return("a")
      }
    }
    
    parameters <- etape_M(data, K, tk, parameters)
    
    pk <- parameters$pk
    moyenne <- parameters$moyenne
    variance <- parameters$variance
    
    for (k in 1:K){
      if (det(variance[,,k]) < 10^-5){
        return(0)
      }
    }
    
  
    LL <- sum(log(apply(sapply(lapply(1:K, 
                                      function(k) pk[k] * dmvnorm(data, moyenne[k,], 
                                                                  variance[,,k])), cbind), 1, sum)))
    
    log_likelihood[i] <- LL
    
    if (abs(log_likelihood[i] - log_likelihood[i-1]) < epislon){
      y_pred <- apply(tk, 1, which.max)
      if (set_inf){
        print("KO")
      }
      return(list(log_likelihood = log_likelihood[-1], y_pred = y_pred))
    }
    
    if (log_likelihood[i] < log_likelihood[i-1]){
      set_inf <- TRUE
    }
    i <- i + 1
  }
  
}



K <- 3
data <- iris[-sample(1:150, 50), -5]
epislon <- 10^-6
data <- as.matrix(data)

c <- 0
for (i in 1:20){
  res <- clustering(data, K, epislon)
  print(typeof(res))
  if (typeof(res) == "list"){ 
    print(res$log_likelihood[length(res$log_likelihood)])
    }
}
