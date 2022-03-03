Mstep <- function(data,K,tk,parameters ) {
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  covariance <- parameters$covariance
  
  P <- ncol(data)
  N <- nrow(data)
  
  
  nk <- apply(tk,2,sum)
  pk <- nk/N;
  moyenne <- t(t(t(data) %*% tk) / nk);
  covariance <- array(0, dim=c(P, P, K));
  for(k in 1:K) {
    sig <- matrix(0,P,P);
    for(n in 1:N) {
      sig <- sig + tk[n, k] * (data[n,] %*% t(data[n,]));
    }
    covariance[,,k] <- sig / nk[k] - moyenne[,k] %*% t(moyenne[,k])
  }  
  return (list(pk=pk, moyenne=t(moyenne), covariance=covariance))
}
# Maximization step: update parameters mixtures, means, covariances, precisions
# Calculating mixtures update.
etape_M2 <- function(data, K, tk, parameters){
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  covariance <- parameters$covariance
  
  P <- ncol(data)
  N <- nrow(data)
  
  pk <- colMeans(tk)
  
  # Calculating means update.
  moyenne <- lapply(1:K, function(j) 
    sapply(1:ncol(data), function(i) 
      apply(tk * data[, i], 2, sum)) 
    [j,]/sum(tk[, j]
    )
  )
  
  # Calculating sample covariances update.
  covariance <- lapply(1:K, function(j) matrix(
    apply(
      sapply(1:N, function(i) 
        tk[i, j] * (data[i, ] - moyenne[[j]]) %*% 
          t(data[i, ] - moyenne[[j]])
      ), 1, sum
    ), P, P)/sum(tk[,j]))
  return (list(pk=pk, moyenne=moyenne, covariance=covariance))
}
etape_M <- function(data, K, tk, parameters){
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  covariance <- parameters$covariance
  
  N <- nrow(data)
  nk <- apply(tk,2,mysum)
  pk <- nk/N
  #columns are features and the lines are the classes
  
  for(k in 1:K) {
    moyenne[k,] <- (1/sum(tk[, k]))*colSums(apply(data,2, function(x) x*tk[, k])  )
    covariance[,,k] = cov.wt(data, wt = tk[, k] , method = "ML")$cov
  }
  return (list(pk=pk, moyenne=moyenne, covariance=covariance))
}

soft.assign <- function(data,  K, parameters){
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  covariance <- parameters$covariance
  
  tauTable <- matrix(NA, nrow = nrow(data), ncol = K)
  
  for (c in 1:ncol(tauTable)){
    tauTable[,c] <- dmvnorm(data, mean=as.numeric(moyenne[c,]), 
                            sigma=covariance[,,c])
    tauTable[,c] <- pk[c] * tauTable[,c]
  }
  tauTable <- tauTable / rowSums(tauTable)  
  return(tauTable)
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