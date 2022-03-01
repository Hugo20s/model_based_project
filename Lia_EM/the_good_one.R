
mysum <- function(x) {
  sum(x[is.finite(x)])
}
etape_E <- function(data, K, parameters){ 
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  covariance <- parameters$covariance
  
  #Sigma -> a square variance-covariance matrix
  tk <- matrix(0,  nrow(data), K)
  for (k in 1:K){
    tk[, k] <- pk[k]*dmvnorm(as.matrix(data), as.numeric(moyenne[k,]), covariance[,,k])
  }
  tk <- tk/apply(tk, 1, mysum)
  return (tk)
}

etape_M <- function(data, K, tk, parameters){
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  covariance <- parameters$covariance
  
  P <- ncol(data)
  N <- nrow(data)
  
  nk <- apply(tk,2,mysum)
  pk <- nk/N
  #columns are features and the lines are the classes
  
  for(k in 1:K) {
    moyenne[k,] <- (1/mysum(tk[, k]))*colSums(apply(data,2, function(x) x*tk[, k])  )
  }
  covariance_aux <- lapply(1:K, function(j) matrix(
    apply(
      sapply(1:N, function(i) 
        tk[i, j] * (data[i, ] - moyenne[[j]]) %*% 
          t(data[i, ] - moyenne[[j]])
      ), 1, mysum
    ), P, P)/mysum(tk[,j]))
  
  covariance <- array(unlist(covariance_aux), dim=c(P,P,K))
  
  return (list(pk=pk, moyenne=moyenne, covariance=covariance))
}


clustering <- function(data, K, epislon){
  
  data <- as.matrix(data)
  
  km.res <- kmeans(data, K)
  pk <- rep(1/K, K)
  moyenne <- as.matrix(km.res$centers)
  
  covariance <- array(0, dim=c(ncol(data),ncol(data),K))
  
  for (i in seq(K)){
    covariance[,,i] <- cov(data[which(as.matrix(km.res$cluster) == i),])
  }
  #Loop 
  i <- 2
  log_likelihood <- numeric(0)
  log_likelihood[1] <- 0
  parameters <- list(pk=pk, moyenne=moyenne, covariance=covariance)
  while (i > 0){
    
    tk <- etape_E(data, K, parameters)
    parameters <- etape_M(data, K, tk, parameters)
    print(parameters)
    pk <- parameters$pk
    moyenne <- parameters$moyenne
    covariance <- parameters$covariance
    
    
    log_k <- numeric(nrow(data))
    for (j in 1:nrow(data)){
      log_k[j] <- logSumExp(lapply(1:K, function(k) {
        log(pk[k]) + logdmvnorm(data[j,], moyenne[k,], covariance[,,k])
      }))
    }
    log_likelihood[i] <- mysum(log_k)
    if (abs(log_likelihood[i] - log_likelihood[i-1]) < epislon){
      break
    }
    print(log_likelihood)
    
    i <- i + 1
  }
  y_pred <- apply(tk, 1, which.max)
  return(list(log_likelihood = log_likelihood[-1], y_pred = y_pred))
  
}
set.seed(3)
K <- 3
data <- iris[-sample(1:150, 50), -5]
epislon <- 10^-6
res <- clustering(data, K, epislon)
plot(res$log_likelihood)
plot(res$y_pred)
tkres<- replicate(5, clustering(data, K, epislon)) 