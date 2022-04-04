
#' @importFrom stats cov
#' @importFrom stats kmeans

## File: EM.r
## Students : Lia Furtado and Hugo Vinson
## Projet Model Based 2021 - Code improvements

## Description : EM algorithm

## Date : 23 February 2022

#Function to calculate/update the variance
fun_variance_all <- function(data, tk, nk, K, N, f, moyenne){
  variance <- array(0, c(f, f, K))
  for(k in 1:K) {
    x_moyenne <- sweep(data, 2, moyenne[k, ], FUN = '-', check.margin=FALSE)
    ex <- tk[,k]  * x_moyenne
    trans_x_moyenne <- t(ex) %*% x_moyenne
    variance[,,k] <- trans_x_moyenne/ nk[k]
  }

  return(variance)
}
#Function to calculate/update the mean
fun_moyenne_vec <- function(data, tk, nk, K, N,  f){
  moyenne <- matrix(0, K, f)
  for(k in 1:K) {
    moyenne[k,] <- (1/sum(tk[,k])) * colSums(apply(data,2, function(x) x*tk[, k])  )
  }
  return(moyenne)
}


######STEP E
etape_E <- function(data, K, parameters){
  pk <- parameters$pk
  moyenne <- parameters$moyenne
  variance <- parameters$variance
  #calculate the probabilities tk
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
  #update the parameters of the model
  pk <- nk/N
  moyenne  <- fun_moyenne_vec(data, tk,nk,K, N,  f)
  variance<- fun_variance_all(data, tk, nk, K, N, f, moyenne)

  return (list(pk=pk, moyenne=moyenne, variance=variance))
}




#' @title EM_algo: Compute the EM algorithm in a dataset
#' @description Implementation of the EM algorithm for the use for a Mixture Model
#' problem
#
#' @param data : main dataset that will be used for the algorithm
#' @param K : number of clusters
#' @param epsilon : stop criteria for a loss decrease lower then this value
#' @param type_init : choose parameter initialization.
#' There are three options: ["random", "kmeans"].
#' "random" is random initialization, "kmeans" is with the kmeans algo
#'
#' @param parameters : model parameters
#'
#' @return : log_likelihood : log likelihood values of the iterations
#'           y_pred : cluster label of each data point
#'           maxll: maximum log likelihood
#'          parameters: the model parameters (mean, variance, tk etc)
#' @export
#'
#' @examples  K <- 5
#' split <- sample(1:150, 100)
#' y <- iris[-split, 5]
#' data <- iris[-split, -5]
#' epsilon <- 10^-6
#' data <- as.matrix(data)
#'
#' res <- EM_algo(data, K, epsilon, "random")
#
#
EM_algo <- function(data, K, epsilon, type_init = "kmeans", parameters = 0){
  set_inf <- FALSE
  an.error.occured <- FALSE
  #--------initialize variables
  n <- nrow(data); p <- ncol(data)
  data <- as.matrix(data)
  if (type_init == "kmeans"){
    parameters <- init_kmeans(data, K)
  }
  if ((type_init == "random"))  {
    parameters <- init_random(data, K)
  }
  i <- 2
  log_likelihood <- numeric(0)
  log_likelihood[1] <- -Inf
  #iterate 100 or when loss decrease is smaller then episilon
  while( (i > 0) || (i > 100))  {
    tk <- etape_E(data, K, parameters)

    parameters <- etape_M(data, K, tk)
    pk <- parameters$pk
    moyenne <- parameters$moyenne
    variance <- parameters$variance
    an.error.occured <- FALSE

    #catch when determinant is equal to null
    tryCatch( {
      LL <- lapply(1:K, function(k) pk[k] * dmultinorm_all(data, moyenne[k,], variance[,,k])) }
      , error = function(e) {an.error.occured <<- TRUE; })
    if (an.error.occured) {
      return(list(maxll = LLsauv, invalid = TRUE))
    }
    LL <- Reduce('+', LL)
    LL <- sum(log(LL))
    #when log-likelihood is equals to null
    if(is.nan(LL)){
      set_inf=TRUE
      return(list(maxll = LLsauv, invalid = set_inf))

    }


    log_likelihood[i] <- LL

    if(LL < log_likelihood[i-1]){
      set_inf=TRUE
      #warning(print("Log-likelihood decreases"))
      return(list(maxll = LLsauv, invalid = set_inf))

    }
    #if loss is smaller then epsilon stops iteration
    if (abs(log_likelihood[i] - log_likelihood[i-1]) < epsilon){
      y_pred <- apply(tk, 1, which.max)
      return(list(log_likelihood = log_likelihood[-1], y_pred = y_pred, maxll = max(log_likelihood[-1]), parameters = parameters, invalid = set_inf))
    }

    LLsauv <- LL

    i <- i + 1

  }

}



