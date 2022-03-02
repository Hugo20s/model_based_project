# EM algorithm for Gaussian mixture model estimation
library("stats")
library("emdbook")
library("matrixStats")
mysum <- function(x) {
  sum(x[is.finite(x)])
}
init_random <- function(data, K){
  f <- ncol(data)
  pk <- runif(K,3/10,7/10)
  pk <- pk / sum(pk)
  
  moyenne <- matrix(0, ncol=ncol(data), nrow = K)
  for (k in 1:K) {
    for (cf in 1:f){
      rand_1 <- runif(1, min(data[,cf]), max(data[,cf]))
      moyenne[k,cf] <- rand_1
    }
  }
  # variance <- array(0, c(f, f, K))
  # for (k in 1:K){
  #   split <- sample(0.40*nrow(data))
  #   print(split)
  #   variance[,,k] <- cov(data[split,])  
  #   print("variance")
  #   print(variance[,,k])
  #   print("___")
  # }
  variance <- fun_var(data, f, K, moyenne)
  
  # print("variance")
  # print(variance)
  
  return (list(pk=pk, moyenne=moyenne, variance=variance))
}


init_kmeans <- function (data, K) {
  km.res <- kmeans(data, K)
  pk <- rep(1/K, K)
  moyenne <- as.matrix(km.res$centers)

  variance <- array(0, dim=c(ncol(data),ncol(data),K))
  for (i in seq(K)){
    variance[,,i] <- cov(data[which(as.matrix(km.res$cluster) == i),])
  }
  
  return (list(pk=pk, moyenne=moyenne, variance=variance))
}


etape_E <- function(data, K, parameters){ 
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  
  tk <- matrix(0,  nrow(data), K)
  for (k in 1:K){
    tk[, k] <- pk[k]* emdbook::dmvnorm(data, moyenne[k,], variance[,,k])
  }
  tk <- tk/apply(tk, 1, sum)
  
  return (tk)
}

etape_M <- function(data, K, tk){
  N <- nrow(data)
  f <- ncol(data)
  nk <- colSums(tk)
  
  pk <- nk/N
  moyenne  <- fun_moyenne_vec(data, tk,nk,K, N,  f)
  variance  <- fun_variance_pkg(data, tk, nk, K, N, f, moyenne)

  return (list(pk=pk, moyenne=moyenne, variance=variance))
}


clustering <- function(data, K, epsilon, type_init = "kmeans", parameters = 0){
  set_inf <- FALSE
  #--------initialize variables
  n <- nrow(data); p <- ncol(data)

  if (type_init == "kmeans"){
    parameters <- init_kmeans(data, K)
  }
  if ((type_init == "random") || (type_init == "small"))  {
    parameters <- init_random(data, K)  
  }

  #Loop 
  i <- 2
  log_likelihood <- numeric(0)
  log_likelihood[1] <- -Inf
  
  
  while (i > 0) {
    tk <- etape_E(data, K, parameters)
    # 
    # for (k in 1:K ){
    #   if (sum(tk[,k]) == 0){
    #     return("a")
    #   }
    # }
    
    parameters <- etape_M(data, K, tk)
    
    pk <- parameters$pk
    moyenne <- parameters$moyenne
    variance <- parameters$variance
    
    # for (k in 1:K){
    #   if (det(variance[,,k]) < 10^-100){
    #     print(det(variance[,,k]))
    #     return(0)
    #   }
    # }

    LL <- sum(log(apply(sapply(lapply(1:K, 
                                      function(k) pk[k] * emdbook::dmvnorm(data, moyenne[k,], 
                                                                  variance[,,k])), cbind), 1, sum)))
    #LL <- fun_ll_all(data, moyennne, variance, pk, n)
    
    log_likelihood[i] <- LL
    
    if (abs(log_likelihood[i] - log_likelihood[i-1]) < epsilon){
      y_pred <- apply(tk, 1, which.max)
      if (set_inf){
        # print(variance)
        # print(moyenne)
        # print(pk)
        print("KO")
        print(log_likelihood[-2])
        print(abs(log_likelihood[i] - log_likelihood[i-1]))
      }
      print(i)
      print(abs(log_likelihood[i] - log_likelihood[i-1]))
      return(list(log_likelihood = log_likelihood[-1], y_pred = y_pred))
    }
    
    if (LL < log_likelihood[i-1]){
      set_inf <- TRUE
      #print(paste("pb", i))
    }
    
    
    if ((type_init == "small") && (i == 5)) {
      print("init fini")
      return( clustering(data, K, epsilon, "", parameters ) )
    }
    
    
    i <- i + 1
  }
  
}



K <- 3
split <- sample(1:150, 50)
y <- iris[-split, 5]
data <- iris[-split, -5]
epsilon <- 10^-6
data <- as.matrix(data)

res <- clustering(data, K, epsilon, "small")
table(res$y_pred, y)


plot(res$log_likelihood)  
plot(res$y_pred)

c <- 0
for (i in 1:100){
  res <- clustering(data, K, epsilon, "small")
  print(typeof(res))
  if (typeof(res) == "list"){ 
    print(res$log_likelihood[length(res$log_likelihood)])
  }
  print("-----")
}

table(res$y_pred, y)
