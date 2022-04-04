
#Function to initialize the EM algo with random values
init_random <- function(data, K){
  f <- ncol(data)
  N <- nrow(data)
  #set.seed(1)
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
#Function to initialize the EM algo with Kmeans algo
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

#Costum function to compute the dmvnorm function
dmultinorm_all <- function(data, moyenne, variance) {
  N <- nrow(data)
  f <- ncol(data)
  coeff <- (1/ ((2*pi)^(f/2) * (det(variance))^0.5))
  variance_inv <- solve.default(variance, tol = 1e-500)
  res <- numeric(N)
  diff <- sweep(data, 2, moyenne, FUN = '-', check.margin=FALSE)
  for (i in 1:N){
    res[i] <- (exp( (-1/2) * (t(diff[i,]) %*% variance_inv %*% diff[i,]))) * coeff
  }

  return(res)
}
loss_function <- function(y_true, y_pred){

  # calculate the metric
  loss <- Accuracy(y_true ,y_pred)

  # convert to tensor
  return(loss)
}
