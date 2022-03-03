# EM algorithm for Gaussian mixture model estimation

library("mvtnorm")

init_random <- function(data, K){
  f <- ncol(data)
  N <- nrow(data)
  pk <- runif(K,3/10,7/10)
  pk <- pk / sum(pk)
  
  moyenne <- matrix(0, ncol=ncol(data), nrow = K)
  for (k in 1:K) {
    for (cf in 1:f){
      rand_1 <- runif(1, min(data[,cf]), max(data[,cf]))
      moyenne[k,cf] <- rand_1
    }
  }
  variance <- array(0, c(f,f,K))
  
  for(k in 1:K) {
    x_moyenne <- sweep(data, 2, moyenne[k, ], FUN = '-', check.margin=FALSE)
    trans_x_moyenne <- crossprod(x_moyenne)
    variance[,,k] <- trans_x_moyenne/ N
  }
  
  return (list(pk=pk, moyenne=moyenne, variance=variance))
}

init_kmeans <- function (data, K) {
  km.res <- kmeans(data, K)
  pk <- rep(1/K, K)
  moyenne <- as.matrix(km.res$centers)
  
  
  variance <- array(0, dim=c(ncol(data),ncol(data),K))
  for (i in seq(K)){
    
    data_cls <- data[which(as.matrix(km.res$cluster) == i),]
    if (length(data_cls) == ncol(data)){
      data_cls <- rbind(data_cls, data_cls+runif(1,-0.1, 0.1))
    } 
    variance[,,i] <- cov(data_cls)
  }
  
  return (list(pk=pk, moyenne=moyenne, variance=variance))
}

etape_E <- function(data, K, parameters){ 
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  
  tk <- matrix(0,  nrow(data), K)
  for (k in 1:K){
    tk[, k] <- pk[k]* mvtnorm::dmvnorm(data, moyenne[k,], variance[,,k])
  }
  tk <- tk/apply(tk, 1, sum)
  
  return (tk)
}

etape_M <- function(data, K, tk){
  N <- nrow(data)
  f <- ncol(data)
  nk <- colSums(tk)
  
  pk <- nk/N
  moyenne <- matrix(0, K, f)
  for(k in 1:K) {
    moyenne[k,] <- (1/sum(tk[,k])) * colSums(apply(data,2, function(x) x*tk[, k])  )
  }
  
  variance <- array(0, c(f, f, K))
  for(k in 1:K) {
    temp <- replace(tk[, k], tk[, k]==0, 10^-8) 
    # variance[,,k] = cov.wt(data, wt = temp , method = "ML")$cov
     
    x_moyenne <- sweep(data, 2, moyenne[k, ], FUN = '-', check.margin=FALSE)
    ex <- temp  * x_moyenne
    trans_x_moyenne <- t(ex) %*% x_moyenne
    variance[,,k] <- trans_x_moyenne/ nk[k]
  
  }

  return (list(pk=pk, moyenne=moyenne, variance=variance))
}


main <- function(data, K, epislon, type_init = "kmeans", parameters = 0){
  set_inf <- FALSE
  #--------initialize variables
  n <- nrow(data); p <- ncol(data)
  data <- as.matrix(data)
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
  
  
  while( (i > 0) || (i > 100))  {
    tk <- etape_E(data, K, parameters)
    
    parameters <- etape_M(data, K, tk)
    
    pk <- parameters$pk
    moyenne <- parameters$moyenne
    variance <- parameters$variance
    LL <- sum(log(apply(sapply(lapply(1:K, 
                                      function(k) pk[k] * mvtnorm::dmvnorm(data, moyenne[k,], 
                                                                           variance[,,k])), cbind), 1, sum)))
    
    log_likelihood[i] <- LL
    
    #print(LL)
    if (is.nan(LL)){
      set_inf <- TRUE
      print("------ LL KO")
      return(list(invalid = set_inf))}
    if (LL < log_likelihood[i-1]){
      set_inf <- TRUE
      print("------ debug KO")
    }
    
    
    if (abs(log_likelihood[i] - log_likelihood[i-1]) < epislon){
      y_pred <- apply(tk, 1, which.max)
      if (set_inf){
        print("___________________________________________________________KO")
        plot(log_likelihood)
      }
      return(list(log_likelihood = log_likelihood[-1], y_pred = y_pred, maxll = max(log_likelihood[-1]), parameters = parameters, invalid = set_inf))
    }
    
    if ((type_init == "small") && (i == 5)) {
      return( main(data, K, epislon, "", parameters ) )
    }
    
    
    i <- i + 1
  }
  
}



# K <- 3
# split <- sample(1:150, 50)
# y <- iris[-split, 5]
# data <- iris[-split, -5]
# epsilon <- 10^-6
# data <- as.matrix(data)
# #
# # res <- main(data, K, epsilon, "small")
# # table(res$y_pred, y)
# #
# #
# # plot(res$log_likelihood)
# # plot(res$y_pred)
# #
# for (i in 1:100){
#   res <- main(data, 4, epsilon, "random")
#   if (typeof(res) == "list"){
#     print(res$log_likelihood[length(res$log_likelihood)])
#   }
#   print("-----")
# }
# #
# table(res$y_pred, y)
