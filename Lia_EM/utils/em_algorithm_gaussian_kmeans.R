#Lia Furtado
#Model- based Learning

# EM algorithm for Gaussian mixture model estimation
library("stats")
library("emdbook")
library("matrixStats")
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
    variance[,,k] = cov.wt(data, wt = tk[, k], method = "ML")$cov
  }
  return (list(pk=pk, moyenne=moyenne, variance=variance))
  }


clustering <- function(data, K, epislon){
 
  #initialize variables
  km.res <- kmeans(data, K)
  pk <- rep(1/K, K)
  moyenne <- as.matrix(km.res$centers)
  
  variance <- array(0, dim=c(ncol(data),ncol(data),K))
  
  for (i in seq(K)){
    variance[,,i] <- cov(data[which(as.matrix(km.res$cluster) == i),])
  }
  #Loop 
  i <- 2
  log_likelihood <- numeric(0)
  log_likelihood[1] <- 0
  parameters <- list(pk=pk, moyenne=moyenne, variance=variance)
  while (i > 0){
    
    tk <- etape_E(data, K, parameters)
    #print(tk)
    parameters <- etape_M(data, K, tk, parameters)
    #print(parameters)
  
    pk <- parameters$pk
    moyenne <- parameters$moyenne
    variance <- parameters$variance
    
    
    
    log_k <- numeric(nrow(data))
    for (j in 1:nrow(data)){
      log_k[j] <- logSumExp(lapply(1:K, function(k) log(pk[k]) + dmvnorm(as.matrix(data[j,]), moyenne[k,], variance[,,k], log = TRUE)))
      
    }
    log_likelihood[i] <- sum(log_k)
    #log_likelihood[i] <- LL
    if (abs(log_likelihood[i] - log_likelihood[i-1]) < epislon){
        break
    }
    #if (i == 100){
      #break
    #}
    i <- i + 1
  }
  y_pred <- apply(tk, 1, which.max)
  return(list(log_likelihood = log_likelihood[-1], y_pred = y_pred, maxll = max(log_likelihood[-1])))
  
}
K <- 3
data <- iris[-sample(1:150, 50), -5]
epislon <- 10^-6
log_likelihood <- clustering(data, K, epislon)
plot(log_likelihood$log_likelihood)
log_likelihood
plot(log_likelihood$y_pred)

tkres<- replicate(5, clustering(data, K, epislon)) 
n <- nrow(data)
p <- ncol(data)
BIC <- numeric(0)
i <- 1
number_clusters <- c(1, 2,3,4,5)
for (k in number_clusters){
  res <- clustering(data, k, epislon)
  BIC[i] <- -2*res$maxll + (k*(p*(p + 1)/2 + p + 1) - 1)*log(n)
  i <- i + 1
}
plot(number_clusters, BIC)

  
