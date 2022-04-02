

tk <- rbind(c(0.2, 0.5, 0.3), c(0.1, 0.8, 0.1), c(0.4, 0.4, 0.2), c(0.2, 0.2, 0.6))
ma <- cbind(c(3,4,3, 6), c(5,5,8, 1), c(5,6,7, 5), c(14,9,10,2), c(11,23,13,9))
nk <- colSums(tk)  

K <- ncol(tk)
f <- ncol(ma)
moyenne <- matrix(0, K, f)
for(k in 1:K) {
  moyenne[k,] <- (1/sum(tk[,k])) * colSums(apply(ma,2, function(x) x*tk[, k])  )
}
moyenne
variance1 <- array(0, c(f, f, K))

for(k in 1:K) {
  x_moyenne <- sweep(ma, 2, moyenne[k, ], FUN = '-', check.margin=FALSE)
  ex <- tk[, k]  * x_moyenne
  trans_x_moyenne <- t(ex) %*% x_moyenne
  variance1[,,k] <- trans_x_moyenne/ nk[k]
}


variance1 <- array(0, c(f, f, K))

for(k in 1:K) {
  x_moyenne <- apply(ma, 1, function(x) x - moyenne[k,])
  ex <- tk[, k]  * x_moyenne
  trans_x_moyenne <- t(ex) %*% x_moyenne
  variance1[,,k] <- trans_x_moyenne/ nk[k]
}
for(k in 1:K) {
  new <- tk[, k]*ne/ nk[k]
  print(new)
}

km.res <- kmeans(ma, K)
pk <- rep(1/K, K)
moyenne <- as.matrix(km.res$centers)


variance <- array(0, dim=c(ncol(data),ncol(data),K))
for (i in seq(K)){
  
  data_cls <- data[which(as.matrix(km.res$cluster) == i),]
  print(length(data_cls))
  print(data_cls)
  print(ncol(data))
  if (length(data_cls) == ncol(data)){
    data_cls <- rbind(data_cls, data_cls+runif(1,-0.1, 0.1))
  } 
  print(data_cls)
  variance[,,i] <- cov(data_cls)}


dmultinorm_all <- function(data, moyenne, variance) {
  N <- nrow(data)
  f <- ncol(data)
  K <- nrow(moyenne)
  
  coeff <- (1/ ((2*pi)^(f/2) * (det(variance))^1/2))
  variance_inv <- solve(variance)
  
  diff <- sweep(data, 2, moyenne, FUN = '-', check.margin=FALSE)
  print(diff)
  res <- numeric(N)
  for (i in 1:N){
    #diff <- data[i,] -  moyenne
    res[i] <- (exp( (-1/2) * (t(diff) %*% variance_inv %*% diff))) * coeff
  }
  
  return(res)
}
library("pracma")
variance = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3)
mu <- c(1,2,3)
x <- matrix(1:6,nrow=2)
res3 <- dmultinorm_all(x, mu, variance)

res1 <- emdbook::dmvnorm(x,mu,M)
res2 <- mvtnorm::dmvnorm(x ,mu,M)
res3 <- dmultinorm_all(x, mu, variance)

LL <-NaN

an.error.occured <- FALSE
tryCatch( {
  if(is.nan(LL) == FALSE){
    
  }
  }, error = function(e) {
    print(e)
    
    an.error.occured <<- TRUE; set_inf=FALSE ; print("Log-likelihood is null ")})
#print(an.error.occured)
if (an.error.occured) {
  print("Yes")
}

BIC <- c(51.61223, 60.93807, 63.84906, 74.24222, 80.07096)
AUX <- c(60.93807, 63.84906, 74.24222 ,80.07096)
which.min(BIC== min(AUX))
which.min(BIC[BIC== ])

min(AUX)
which(BIC== min(AUX))

