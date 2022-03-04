# EM algorithm for Gaussian mixture model estimation
library("mvtnorm")
library("emdbook")
#install.packages('matrixStats')
library("matrixStats")

####### INIT
init_0 <- function(data, K){
  f <- ncol(data)
  N <- nrow(data)
  set.seed(1)
  pk <- rep(1/K, K)
  #pk <- pk / sum(pk)
  
  moyenne <- matrix(0, ncol=ncol(data), nrow = K)
  
  variance <- array(0, c(f,f,K))
  for (k in 1:K){
    moyenne[k,] <- rnorm(f, mean = 0)
    variance[,,k] <- diag(f)
  }
  return (list(pk=pk, moyenne=moyenne, variance=variance))
  
}

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


#######FONCTION MOYENNE VARIANCE DMVNORM
fun_variance_all <- function(data, tk, nk, K, N, f, moyenne){
  variance <- array(0, c(f, f, K))
  for(k in 1:K) {
    x_moyenne <- sweep(data, 2, moyenne[k, ], FUN = '-', check.margin=FALSE)
    ex <- tk[, k]  * x_moyenne
    trans_x_moyenne <- t(ex) %*% x_moyenne
    variance[,,k] <- trans_x_moyenne/ nk[k]
  }
  
  return(variance)
}

fun_moyenne_vec <- function(data, tk, nk, K, N,  f){
  moyenne <- matrix(0, K, f)
  for(k in 1:K) {
    moyenne[k,] <- (1/sum(tk[,k])) * colSums(apply(data,2, function(x) x*tk[, k])  )
  }
  return(moyenne)
}

dmultinorm_all <- function(data, moyenne, variance) {
  N <- nrow(data)
  f <- ncol(data)
  coeff <- (1/ ((2*pi)^(f/2) * (det(variance))^0.5))
  variance_inv <- solve.default(variance)
  res <- numeric(N)
  diff <- sweep(data, 2, moyenne, FUN = '-', check.margin=FALSE)
  for (i in 1:N){
    res[i] <- (exp( (-1/2) * (t(diff[i,]) %*% variance_inv %*% diff[i,]))) * coeff
  }
  
  return(res)
}


######STEP E 
etape_E <- function(data, K, parameters){ 
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  
  tk <- matrix(0,  nrow(data), K)
  for (k in 1:K){
    tk[, k] <- pk[k]* dmultinorm_all(data, moyenne[k,], variance[,,k])
  }
  tk <- tk/apply(tk, 1, sum)

  return (tk)
}
######STEP M 
etape_M <- function(data, K, tk){
  N <- nrow(data)
  f <- ncol(data)
  nk <- colSums(tk)
  
  pk <- nk/N
  moyenne  <- fun_moyenne_vec(data, tk,nk,K, N,  f)
  variance<- fun_variance_all(data, tk, nk, K, N, f, moyenne)
  
  return (list(pk=pk, moyenne=moyenne, variance=variance))
}


main <- function(data, K, epsilon, type_init = "kmeans", parameters = 0){
  set_inf <- FALSE
  an.error.occured <- FALSE
  #--------initialize variables
  n <- nrow(data); p <- ncol(data)
  data <- as.matrix(data)
  if (type_init == "kmeans"){
    parameters <- init_kmeans(data, K)
  }
  if ((type_init == "random") || (type_init == "small"))  {
    parameters <- init_random(data, K)  
  }
  if (type_init == "init0"){
    parameters <- init_0(data, K)
  }
  
  #Loop 
  i <- 2
  log_likelihood <- numeric(0)
  log_likelihood[1] <- -Inf
  
  while( (i > 0) || (i > 100))  {
    #print(i)
    tk <- etape_E(data, K, parameters)
    
    parameters <- etape_M(data, K, tk)
    pk <- parameters$pk
    moyenne <- parameters$moyenne
    variance <- parameters$variance
    LL1 <- sum(log(apply(sapply(lapply(1:K,
                                      function(k) pk[k] * mvtnorm::dmvnorm(data, moyenne[k,],
                                                                  variance[,,k])), cbind), 1, sum)))
    log_lik <- numeric(n)
    for (j in 1:nrow(data)){
      log_lik[j] <- logSumExp(lapply(1:K, function(k) 
        log(pk[k]) + dmvnorm(data[j,], moyenne[k,], variance[,,k], log = TRUE)
      ))
    }

    LL2 <- sum(log_lik)
    
    an.error.occured <- FALSE
    tryCatch( { LL <- lapply(1:K, function(k) pk[k] * dmultinorm_all(data, moyenne[k,], variance[,,k])) }
              , error = function(e) {an.error.occured <<- TRUE; set_inf=FALSE ; print("ko")})
    if (an.error.occured) {
      return(list(maxll = LLsauv, invalid = set_inf))
    }
    #LL <- lapply(1:K, function(k) pk[k] * dmultinorm_all(data, moyenne[k,], variance[,,k]))
    LL <- Reduce('+', LL)
    LL <- sum(log(LL))
    # print("LL")
    # print(LL)
    # print("LL1")
    # print(LL1)
    # print("LL2")
    # print(LL)
    if (is.na(LL)){
      print(paste('itration',i))
      print("tk")
      print(tk)
      print("sum tk---------")
      print("new")
      print(colSums(tk))
      # print("old")
      # print(colSums(tk_sauv))
      print("pk------------")
      print("new")
      print(pk)
      # print("old")
      # print(param_sauv$pk)
      print("moyenne--------")
      print("new")
      print(moyenne)
      # print("old")
      # print(param_sauv$moyenne)
      print("variance-----------")
      print("new")
      print(variance)
      # print("old")
      # print(param_sauv$variance)
      for (k in 1:K){
        print(paste("dmvnorm",k))
        print("new")
        print(dmultinorm_all(data, moyenne[k,], variance[,,k]))
        # print("old")
        # print(list_densite[[k]])
      }
    }
    
    # print(paste("LL", LL))
    # print(paste("LLgood", LL1))
    log_likelihood[i] <- LL
    
    if (LL < log_likelihood[i-1]){
      set_inf <- TRUE
      print("------ debug KO")
      return(list(invalid = TRUE))
      # print("tk")
      # print(tk)
      # print("sum tk---------")
      # print("new")
      # print(colSums(tk))
      # print("old")
      # print(colSums(tk_sauv))
      # print("pk------------")
      # print("new")
      # print(pk)
      # print("old")
      # print(param_sauv$pk)
      # print("moyenne--------")
      # print("new")
      # print(moyenne)
      # print("old")
      # print(param_sauv$moyenne)
      # print("variance-----------")
      # print("new")
      # print(variance)
      # print("old")
      # print(param_sauv$variance)
      # for (k in 1:K){
      #   print(paste("dmvnorm",k))
      #   print("new")
      #   print(dmultinorm_all(data, moyenne[k,], variance[,,k]))
      #   print("old")
      #   print(list_densite[[k]])
      # }
    }
    
    if (abs(log_likelihood[i] - log_likelihood[i-1]) < epsilon){
      y_pred <- apply(tk, 1, which.max)
      if (set_inf){
        print("___________________________________________________________KO")
        plot(log_likelihood)
      }
      return(list(log_likelihood = log_likelihood[-1], y_pred = y_pred, maxll = max(log_likelihood[-1]), parameters = parameters, invalid = set_inf))
    }
    
    if ((type_init == "small") && (i == 5)) {
      return( main(data, K, epsilon, "", parameters ) )
    }
    
    LLsauv <- LL
    parameters_sauv <- parameters
    # tk_sauv <- tk 
    # list_densite <- list()
    # list_densite_all <- list()
    # for (k in 1:K){
    #   list_densite[[k]] <- mvtnorm::dmvnorm(data, moyenne[k,], variance[,,k])
    #   #list_densite_all[[k]] <- dmultinorm_all(data, moyenne[k,], variance[,,k], p)
    # }
    # for (k in 1:K){
    #   print("emd")
    #   print(emdbook::dmvnorm(data, moyenne[k,], variance[,,k]))
    #   print("all")
    #   #print(dmultinorm_all(data, moyenne[k,], variance[,,k], p))
    #   print("mv")
    #   print(dmultinorm_all(data, moyenne[k,], variance[,,k]))
    # }
    i <- i + 1
    #return()

  }
  
}



K <- 3
split <- sample(1:150, 100)
y <- iris[-split, 5]
data <- iris[-split, -5]
epsilon <- 10^-6
data <- as.matrix(data)

res <- main(data, K, epsilon, "init0")
table(res$y_pred, y)
plot(res$log_likelihood)  
plot(res$y_pred)


data_test <- data[which(y=="setosa"),]


for (i in 1:30){
  res <- main(data, K, epsilon, "init0")
  if (typeof(res) == "list"){ 
    print(res$maxll)
    print(length(res$log_likelihood))
  }
  print("-----")
}






