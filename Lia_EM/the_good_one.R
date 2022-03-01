
library("MLmetrics")
library("stats")
library("emdbook")
library("matrixStats")
library("mixtools")

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
      sapply(1:N, function(i) {
        temp <- replace(tk[i, j], tk[i, j]==0, 10^-8)
        temp * (data[i, ] - moyenne[[j]]) %*% 
          t(data[i, ] - moyenne[[j]])}
      ), 1, mysum
    ), P, P)/mysum(tk[,j]))
  
  covariance <- array(unlist(covariance_aux), dim=c(P,P,K))
  
  return (list(pk=pk, moyenne=moyenne, covariance=covariance))
}

LL <- function(data, K, parameters){
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  covariance <- parameters$covariance
  
  # Calculate log-likelihood for mixture model
  LL <- sum(log(apply(sapply(lapply(1:K, 
                                    function(i) pk[i] * dmvnorm(data, moyenne[i,], 
                                                                covariance[,,i])), cbind), 1, sum)))
  return(LL)
  
}
init <- function(data, K){
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
  return (list(pk=pk, moyenne=moyenne, covariance=covariance))
  
}
main <- function(data, K, epislon){
  
  data <- as.matrix(data)
  
  #initialize variables
  parameters <- init(data, K)
  #Loop 
  i <- 2
  log_likelihood <- numeric(0)
  log_likelihood[1] <- 0
  while (i > 0){
    
    print("Etape E")
    tk <- etape_E(data, K, parameters)
    parameters <- etape_M(data, K, tk, parameters)

    LL <- LL(data, K, parameters)
    log_likelihood <- c(log_likelihood, LL)

        if (abs(log_likelihood[i] - log_likelihood[i-1]) < epislon){
      break
    }
    i <- i + 1
  }
  y_pred <- apply(tk, 1, which.max)
  return(list(log_likelihood = log_likelihood[-1], y_pred = y_pred))
  
}


#set.seed(3)
K <- 3
data <- iris[-sample(1:150, 50), -5]
y <- iris[-sample(1:150, 50), 5]
epislon <- 10^-6
res <- clustering(data, K, epislon)
plot(res$log_likelihood)
plot(res$y_pred)
tkres<- replicate(20, clustering(data, K, epislon)) 

res <- do.call(rbind, tkres[1,])
best_iteration <- which.max(res[,ncol(res)])

y_pred <- tkres[2, best_iteration]$y_pred

plot(1:length(as.list(y_pred)), y_pred)
#Accuracy(y_pred, y)


l <- res$log_likelihood
for (i in seq(length(l)-1)){
  if (l[i] < l[i]+1) {
    print('ok')
  } else {
    print('not ok')
  }
}

