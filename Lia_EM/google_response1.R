library("glasso")
library("mvtnorm")

# Initialize parameters in order to start the EM algorithm
init <- function(data, nclust){  
  mixtures <- c()
  means <- list()
  covariances <- list()
  
  init.assign <- sample(1:nclust, nrow(data), replace=T)
  
  for(k in 1:nclust){
    mixtures[k] <- table[names(table)==k] / nrow(data)
    means[[k]] <- as.vector(colMeans(data[which(init.assign %in% k),]))
    covariances[[k]] <- cov(data[which(init.assign %in% k),])
  }
  return(list(mixtures = mixtures, means = means, covariances = covariances))
}

# Calculate posterior probabilities tau_ik
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

# Calculate (penalized) log-likelihood given model parameters
LL <- function(data, lambda, mixtures, means, covariances, penalize = F){
  # Calculate log-likelihood for mixture model
  LL <- sum(log(apply(sapply(lapply(1:length(mixtures), 
                                    function(i) mixtures[i] * dmvnorm(data, means[[i]], 
                                                                      covariances[[i]])), cbind), 1, sum)))
  
  # Calculate and apply penalty through parameter lambda
  if(penalize){
    O <- lapply(covariances, solve)
    sum <- sum(sapply(1:length(mixtures), function(i) mixtures[i] * 
                        sum(abs(O[[i]]))))
    
    penalty <- (nrow(data)/2) * lambda * sum
    pLL <- LL - penalty
    return(pLL)
  }
  else{ return(LL) }
}

# Main EM function. Returns parameters at convergence.
EM <- function(data, lambda, mixtures, means, covariances, 
               iterations = 100, threshold = 10^-4){
  data <- as.matrix(data)
  
  K <- length(mixtures)
  P <- ncol(data)
  N <- nrow(data)
  
  deltaLL <- 1
  
  prevLL <- LL(data, lambda, mixtures, means, covariances, penalize=T)
  LLpath <- c(prevLL)
  
  mixtures <- mixtures
  means <- means
  covariances <- covariances
  precisions <- list()
  
  increment <- 1
  while(increment < iterations && deltaLL > threshold){
    
    # Expectation step: calculate posterior probabilities tau_{ik}
    tauTable <- etape_E(data, mixtures, means, covariances)
    
    # Maximization step: update parameters mixtures, means, covariances, precisions
    # Calculating mixtures update.
    newMixtures <- colMeans(tauTable)
    
    # Calculating means update.
    newMeans <- lapply(1:K, function(j) 
      sapply(1:18, function(i) 
        apply(tauTable * data[, i], 2, sum)) 
      [j,]/sum(tauTable[, j]
      )
    )
    
    # Calculating sample covariances update.
    sampleCov <- lapply(1:K, function(j) matrix(
      apply(
        sapply(1:N, function(i) 
          tauTable[i, j] * (data[i, ] - newMeans[[j]]) %*% 
            t(data[i, ] - newMeans[[j]])
        ), 1, sum
      )
      , P, P)/sum(tauTable[,j]))
    
    
    # Calculating updates for covariance estimates & precision matrix estimates.
    newEstimates <- list()
    newPrecisions <- list()
    for (k in 1:K){
      gl <- glasso(s = sampleCov[[k]], rho = lambda, penalize.diagonal = T)
      newEstimates[[k]] <- gl$w
      newPrecisions[[k]] <- gl$wi
    }
    
    mixtures <- newMixtures
    means <- newMeans
    covariances <- newEstimates
    precisions <- newPrecisions
    
    # Calculate model penalized log-likelihood and compare with previous
    
    pLL <- LL(data, lambda, mixtures, means, covariances, penalize=T)
    print(pLL)
    deltaLL <- abs((pLL / prevLL) - 1)
    prevLL <- pLL
    LLpath <- c(LLpath, pLL)
    
    increment <- increment + 1
    
  }
  return(list(
    LL = prevLL, 
    LLpath = LLpath, 
    soft.assign = tauTable, 
    mixtures = mixtures, 
    means = means, 
    covariances = covariances, 
    precisions = precisions))
}

K <- 3
data <- iris[-sample(1:150, 50), -5]
epislon <- 10^-6
result <- init(data, K)
mixtures
EM(data, lambda, result$mixtures, result$means, result$covariances, iterations = 100, threshold = 10^-4)
