#Lia Furtado - Hugo Vinson
#Model- based Learning

#install.packages("emdbook")
library(emdbook)
#install.packages("stats")
library(stats)

# EM algorithm for Gaussian mixture model estimation

etape_E <- function(data, K, parameters){
  #step E of EM algorithm
  prop <- parameters$prop
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  
  n<- nrow(data); p <- ncol(data)
  
  tk <- matrix(0, n, K)
  temp <- matrix(0, n, K)
  for (k in 1:K){
    temp[, k] <- prop[k]*dmvnorm(data, moyenne[,k], variance[,,k])
  }
  tk <- temp/rowSums(temp)
  return (tk)
}

etape_M <- function(data, K, tk, parameters){
  #step M of EM algorithm
  n<- nrow(data); p <- ncol(data)
  
  prop <- parameters$prop
  moyenne <- parameters$moyenne
  variance <- parameters$variance

  for (k in 1:K){
    prop[k] <- sum(tk[, k])/ nrow(data)
    moyenne[,k] <- (1/sum(tk[, k]))*colSums( apply(data,2, function(x) x*tk[, k])  )
    #moyenne[,k] <- weighted.mean(data, t[, k])
    #a <- (1/sum(tk[, k]))
    #b <- t((tk[, k])*(data - moyenne[,k]))%*% (data - moyenne[,k])
    #variance[,,k] <- a*b
    variance[,,k] <- cov.wt((data),  tk[,k])$cov
  
  }
  return (list(prop=prop, moyenne=moyenne, variance=variance))
  }



clustering <- function(data, K, epislon){
  
  n <- nrow(data); p <- ncol(data)
  
  #initialize variables
  init <- kmeans(data, K)
  prop <- numeric(K)
  for (k in 1:K){
    prop[k] <- length(which(init$cluster==k))/n
  }
  moyenne <- t(as.matrix(init$centers))
  variance <- array(0, dim=c(p,p,K))
  for (i in seq(K)){
    matrix_temp <- matrix(runif(p^2, 0, 1), p, p)
    matrix_temp <- t(matrix_temp) %*% matrix_temp
    variance[,,i] <- matrix_temp
  }
  
  #Loop 
  i <- 2
  log_likelihood <- numeric(0)
  log_likelihood[1] <- 0
  parameters <- list(prop=prop, moyenne=moyenne, variance=variance)
  while (i > 0){
    
    tk <- etape_E(data, K, parameters)
    parameters <- etape_M(data, K, tk, parameters)
    
  
    prop <- parameters$prop
    moyenne <- parameters$moyenne
    variance <- parameters$variance
  
    #print
    print('----')
    print(prop)
    print(moyenne)
    print(variance)
    print(tk)
    print(colSums(tk))
    
    log_likelihood[i] <- 0
    log_k <- numeric(K)
    #print(variance)
    for (j in 1:nrow(data)){
      log_k[j] <- logSumExp(lapply(1:K, function(k) 
      {
        print(as.matrix(data[j,]))
        print(variance[,,k])
        log(pk[k]) + dmvnorm(as.matrix(data[j,]), moyenne[k,], variance[,,k], log = TRUE)
      }))
    }
    log_likelihood[i] <- sum(log_k)
    
    
    #if (abs(log_likelihood[i] - log_likelihood[i-1]) < epislon){
        #break
    #}
    if (i == 30){
      break
    }
    i <- i + 1
  }
  plot(log_likelihood[-1])
  return(list(ll = log_likelihood[-1], pred = apply(tk, 1, which.max )))
}

K <- 3
data <- iris[-sample(1:150, 50), -5]
data <- as.matrix(data)

epislon <- 10^-6

log_likelihood <- clustering(as.matrix(data), K, epislon)

plot(log_likelihood$pred)


clus_1 <- matrix(rnorm(20*4, 0, 1), 20, 4)
clus_2 <- matrix(rnorm(20*4, 5, 2), 20, 4)
clus_3 <- matrix(rnorm(20*4, 10, 2), 20, 4)

data <- rbind(clus_1, clus_2, clus_3)



is.unsorted()

undebug(clustering)
undebug(etape_M)
undebug(etape_E)

dmvnorm(data, moyenne[, k], variance[, , k])


test_1 <- matrix(runif(30), 6, 5)
test_2 <- c(1,2,2,2,2) 
2*test_2

test_2*test_1

rep(3, c(1,2))
help(rep)
