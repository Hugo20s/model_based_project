# EM algorithm for Gaussian mixture model estimation
library("stats")
library("emdbook")
library("matrixStats")
mysum <- function(x) {
  sum(x[is.finite(x)])
}
etape_E <- function(data, K, parameters){ 
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  
  
  tk <- matrix(0,  nrow(data), K)
  for (k in 1:K){
    #print("density")
    #print(k)
    #print(pk[k])
    #print(dmvnorm(data, moyenne[k,], variance[,,k]))
    #print(dmvnorm(data, moyenne[k,], variance[,,k]))
    tk[, k] <- pk[k]*dmvnorm(data, moyenne[k,], variance[,,k])
  }
  tk <- tk/apply(tk, 1, mysum)
  # print('-----')
  # print(tk[, k])
  # print('-----')
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
  #columns are features and the lines are the classes
  #moyenne <- t(t(data)%*%tk/nk)
  #Sigma -> a square variance-covariance matrix
  
  moyenne2 <- t(t(data) %*% tk) / nk #dimension(k,f)
  print("moyenne2")
  print(moyenne2)
  for(k in 1:K) {
    moyenne[k,] <- (1/mysum(tk[, k]))*colSums(apply(data,2, function(x) x*tk[, k])  )
  }
  
  covariance_aux <- lapply(1:K, function(j) matrix(
    apply(
      sapply(1:N, function(i) {
        temp <- replace(tk[i, j], tk[i, j]==0, 10^-8)
        temp * (data[i, ] - moyenne[[j]]) %*%
          t(data[i, ] - moyenne[[j]])}
      ), 1, mysum
    ), f, f)/mysum(tk[,j]))

  variance <- array(unlist(covariance_aux), dim=c(f,f,K))

  # for(k in 1:K) {
  #   temp <- 0
  #   for(n in 1:N) {
  #     temp <- temp + tk[n, k] * t(data[n,] - moyenne[k,]) * (data[n,] - moyenne[k,])}
  #   print(temp)
  #   #variance[,,k] <- temp/nk[k]
  # }
  
  return (list(pk=pk, moyenne=moyenne, variance=variance))
}
  

clustering <- function(data, K, epislon){
  
  #--------initialize variables
  pk <- c(0.2, 0.5, 0.3)

  moyenne <- matrix(0, ncol=ncol(data), nrow = K)

  for (k in 1:K) {
    rand_1 <- runif(1, min(data[,1]), max(data[,1]))
    rand_2 <- runif(1, min(data[,2]), max(data[,3]))
    rand_3 <-   runif(1, min(data[,3]), max(data[,3]))
    rand_4 <- runif(1, min(data[,4]), max(data[,4]))
    moyenne[k,] <- c(rand_1, rand_2, rand_3, rand_4)
  }

  variance <- array(0, dim=c(ncol(data),ncol(data),K))
  for (i in seq(K)){
    matrix_temp <- matrix(runif((ncol(data)^2)), (ncol(data)), ncol(data))
    matrix_temp[lower.tri(matrix_temp)] = t(matrix_temp)[lower.tri(matrix_temp)]
    matrix_temp <- t(matrix_temp) %*% matrix_temp
    variance[,,i] <- matrix_temp
  }
  
  #initialize variables
  # km.res <- kmeans(data, K)
  # pk <- rep(1/K, K)
  # moyenne <- as.matrix(km.res$centers)
  # 
  # variance <- array(0, dim=c(ncol(data),ncol(data),K))
  # 
  # for (i in seq(K)){
  #   variance[,,i] <- cov(data[which(as.matrix(km.res$cluster) == i),])
  # }
  
  print("PK")
  print(pk)
  print("sum pk")
  print(sum(pk))
  print("moyenne:")
  print(moyenne)
  print("variance: ")
  print(variance)
  
  print("-------LOOP------")
  
  
  #Loop 
  i <- 2
  log_likelihood <- numeric(0)
  log_likelihood[1] <- 0
  parameters <- list(pk=pk, moyenne=moyenne, variance=variance)
  while (i > 0){
    
    print(paste("ETAPE E", i))
    tk <- etape_E(data, K, parameters)
    print("tk")
    print(tk)
    print(paste("ETAPE M", i))
    parameters <- etape_M(data, K, tk, parameters)
    
    pk <- parameters$pk
    moyenne <- parameters$moyenne
    variance <- parameters$variance
    
    log_k <- numeric(nrow(data))
    print("PK")
    print(pk)
    print("moyenne:")
    print(moyenne)
    print("variance: ")
    print(variance)
    
    print(paste("ETAPE LOGLIKELIHOOD", i))
    LL <- sum(log(apply(sapply(lapply(1:K, 
                                function(k) pk[k] * dmvnorm(data, moyenne[k,], 
                                                            variance[,,k])), cbind), 1, sum)))
    # for (j in 1:nrow(data)){
    #   # for (k in 1:K){
    #   #   log_k[j] <- lof(logSumExp(
    #   # }
    #   
    #   log_k[j] <- logSumExp(lapply(1:K, function(k) {
    #     pk_modif <- replace(pk[k], pk[k]==0, 10^-8)
    #     log(pk_modif) + dmvnorm(data[j,], moyenne[k,], variance[,,k], log = TRUE)
    #   }))
      #   # print("log")
      #   pk_modif <- replace(pk[k], pk[k]==0, 10^-8)
      #   # print(log(pk_modif))
      #   # print("densité")
      #   #print(dmvnorm(data[j,], moyenne[k,], variance[,,k], log = TRUE))
      #   log(pk_modif) + dmvnorm(data[j,], moyenne[k,], variance[,,k], log = TRUE)
      # }))
    #}
    #log_likelihood[i] <- sum(log_k)
    log_likelihood[i] <- LL
    print("log_likelihood")
    print(log_likelihood[i])
    
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
data <- as.matrix(data)
res <- clustering(data, K, epislon)

plot(res$log_likelihood)
plot(res$y_pred)
tkres<- replicate(5, clustering(data, K, epislon)) 

mean((res - data$pred)^2)
