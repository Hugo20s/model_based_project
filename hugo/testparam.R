K <- 3
P <- 2
N <- 3
#n = 3
#f = 2
#k = 3


d <- cbind(c(1,2,3), c(4,5,6))    #n;f
t <- cbind(c(0.1,0.2,0.3), c(0.7,0.7,0.4), c(0.2, 0.1, 0.3))         #n;k
st <- colSums(t) #k
moyenne <- t(t(d) %*% t) / st #k,f


#TEST MOYENNE 
moyenne <- matrix(0, K, P)
for(k in 1:K) {
  moyenne[k,] <- (1/mysum(t[, k]))*colSums(apply(d,2, function(x) x*t[, k])  )
}

moyenne <- t(t(d) %*% t) / st #k,f



#TEST VARIANCE 
variance <- array( 0, c(2, 2, K))   # f;f;k
for (k in 1:K) {
  temp1 <- (d - moyenne[k,])  #n;f    
  temp2 <- t(t[,k] * temp1) %*% temp1     # t(n;f * n)  %*% n;f
  variance[,,k] <- temp2 / st[k]
}

covariance_aux <- lapply(1:K, function(j) matrix(
  apply(
    sapply(1:N, function(i) 
      t[i, j] * (d[i, ] - moyenne[j]) %*% 
        t(d[i, ] - moyenne[j])
    ), 1, mysum
  ), P, P)/mysum(t[,j]))

covariance <- array(unlist(covariance_aux), dim=c(P,P,K))

variance <- array( 0, c(2, 2, K)) 
for (k in 1:K){
  res <- matrix(0, 2, 2)
  for (i in 1:N){
    res <- res + t[i, k] * (d[i, ] - moyenne[k]) %*% t(d[i, ] - moyenne[k])
  }
  variance[,,k] <- res / mysum(t[,k])
}

for (k in 1:K){
  
}

#TEST LIKELIHOOD

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

pk <- c(0.4, 0.3, 0.3)
res <- matrix(0, 3, 3)
for (i in 1:K){
  res[i,] <- pk[i] * dmvnorm(d, moyenne[i,], variance[,,i])
}
res <- rowSums(res)
res <- log(res)
ll <- sum(res)



i <- 1
k <- 1

apply(temp1,2, function(x) x*t[, k])



###ANCIENNE METHODE
for (k in 1:K){
  #moyenne[k,] <- (1/sum(tk[, k]))*colSums( apply(data,2, function(x) x*tk[, k])  )
  moyenne[k,] <- (1/nk)*colSums(tk[, k] * data)
  wt <- replace(tk[, k], tk[, k]==0, 10^-8) 
  variance[,,k] = cov.wt(data, wt = wt , method = "ML")$cov
}
