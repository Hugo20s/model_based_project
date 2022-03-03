fun_variance_all <- function(data, tk, nk, K, N, f, moyenne){
  variance <- array(0, c(f,f,K))
  for (k in 1:K){
    sig <- matrix(0, f, f)
    for (i in 1:N) {
      
      sig <- sig + tk[i,k] * (data[i,] - moyenne[k,]) %*% t(data[i,] - moyenne[k,])
    }
    variance[,,k] <- sig / nk[k]
  }
  return(variance)
}


fun_var <- function(data, f, K, moyenne){
  N <- nrow(data)
  variance <- array(0, c(f,f,K))
  for (k in 1:K){
    sig <- matrix(0, f, f)
    for (i in 1:N) {
      
      sig <- sig + (data[i,] - moyenne[k,]) %*% t(data[i,] - moyenne[k,])
    }
    variance[,,k] <- sig / N
  }
  return(variance)
}

fun_variance_guy <- function(data, tk, nk, K, N, f, moyenne){
  covariance_aux <- lapply(1:K, function(j) matrix(
    apply(
      sapply(1:N, function(i) {
        temp <- replace(tk[i, j], tk[i, j]==0, 10^-8)
        temp * (data[i, ] - moyenne[[j]]) %*%
          t(data[i, ] - moyenne[[j]])}
      ), 1, mysum
    ), f, f)/mysum(tk[,j]))
  variance <- array(unlist(covariance_aux), dim=c(f,f,K))
  return(variance)
}


fun_variance_pkg <- function(data, tk, nk, K, N, f, moyenne){
  variance <- array(0, c(f, f, K))
  for(k in 1:K) {
    variance[,,k] = cov.wt(data, wt = tk[,k] , cor = TRUE, method = "ML")$cov
  }
  return(variance)
}


fun_moyenne_all <- function(data, tk, nk, K, N, f){
  moyenne <- matrix(0, K, f)
  for (k in 1:K){
    moy <- numeric(f)
    for (i in 1:N){
      moy <- moy + tk[i,k] * data[i]
    }  
    moyenne[k,] <- moy / nk[k]
  }
  return(moyenne)
}


fun_moyenne_mat <- function(data, tk, nk, K, N, f){
  moyenne <- t(t(data) %*% tk) / nk 
  return(moyenne)
}

fun_moyenne_vec <- function(data, tk, nk, K, N,  f){
  moyenne <- matrix(0, K, f)
  for(k in 1:K) {
    moyenne[k,] <- (1/sum(tk[,k])) * colSums(apply(data,2, function(x) x*tk[, k])  )
  }
  return(moyenne)
}

dmultinorm_emd <- function(data, moyenne, variance) {
  res <- emdbook::dmvnorm(data, moyenne, variance)
  return(res)
}

dmultinorm_mvt <- function(data, moyenne, variance) {
  res <- mvtnorm::dmvnorm(data, moyenne, variance)
  return(res)
}

dmultinorm_all <- function(data, moyenne, variance) {
  data2 <- apply(data, 1, function(x) x - moyenne)
  res <- (1/ ((2*pi)^(f/2) * det(sigma)^1/2)) * exp(-1/2 * t(data2) %*% solve(variance) %*% data2)
  return(res)
}



fun_ll_all <- function(data, moyennne, variance, pk, N){
  ll <- numeric(0)
  for (i in 1:N){
    for (k in 1:K){
      ll <- ll + pk[k] * dmvnorm(data[i,], moyenne[k,],variance[,,k]) 
    }
  }
  return(ll)
}


fun_ll_guy <- function(data, moyenne, variance, pk){
  LL <- sum(log(apply(sapply(lapply(1:K, 
                                    function(k) pk[k] * dmvnorm(data, moyenne[k,], 
                                                                variance[,,k])), cbind), 1, sum)))
  return(LL)
}
  