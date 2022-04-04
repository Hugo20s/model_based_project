
#' @importFrom MLmetrics  Accuracy
#' @importFrom emdbook dmvnorm
#' @importFrom stats kmeans
#' @importFrom stats runif
#'
## File: classification.R
## Students : Lia Furtado and Hugo Vinson
## Projet Model Based 2021 - Code improvements

## Description : train and predict a classification model
## using a Mixture of Mixture model

## Date : 23 February 2022

#-------------- Functions ----------------------

#Choose the best possible models configurations for the smallest BIC
choose_best_mixture <- function(data_part, number_clusters, n_interaction, epislon ){
  #----------------------------------- Clustering Train  -------------------
  n <- nrow(data_part)
  p <- ncol(data_part)
  BIC <- numeric(0)
  results <- list()
  i <- 1

  for (k in number_clusters){
    list_ll <- numeric()
    list_res <- list()
    try <- 1
    for (interaction in  1:n_interaction){
      res <- EM_algo(data_part, k, epislon, "random")
      # ALGO EM
      if (res$invalid == FALSE){
        list_ll[try] <- res$maxll
        list_res[[try]] <- res
        try <- try + 1
      }
    }
    if (!(length(list_res) == 0) ){
      res <- list_res[[which.max(list_ll)]]

      BIC[i] <- -2*res$maxll + (k*(p*(p + 1)/2 + p + 1) - 1)*log(n)

      results[[i]] <- res
      i <- i + 1
    }

  }
  return (list(BIC = BIC, results = results))

}

#' @title MixtureMixture_train
#' @description  Train a classification model by picking the model parameters
#' that give as output a highest accuracy and a number of cluster with
#' the smallest BIC value
#' @param X_train : matrix of data of dimension individus x feature
#' @param y_train : vector of class of dimension individus x feature
#' @param number_clusters : max number of cluster to find model inside class
#' @param n_interaction : number of gaussian model create to find max likelihood
#' @param epislon : difference between two likelihood to stop the model
#' @return list_model: returns the gaussian model to each class
#' @export
#'
#' @examples rand<- sample(1:150, 150)
#' X_train <- iris[rand[1:100], -5]
#' y_train <- iris[rand[1:100], 5]
#' y_train <- as.numeric(factor(y_train))
#' X_test <- iris[rand[101:150], -5]
#' y_test <- iris[rand[101:150], 5]
#' y_test <- as.numeric(factor(y_test))
#' epislon <- 10^-6
#' number_clusters <- 1:3
#' n_interaction <- 50
#' results <- MixtureMixture.train (X_train, y_train, number_clusters, n_interaction, epislon)
#'
MixtureMixture.train <- function(X_train, y_train, number_clusters, n_interaction, epislon){

  classes <- unique(y_train)
  best_models <- list()
  map_cluster_class <- array(NA)
  old_k <- 0
  for (n_class in classes){
    best_k <- numeric(0)
    data_part <- X_train[which(y_train == classes[n_class]),]

    #---------------- Choosing best mixture model ----------------- #

    best_mixture <- choose_best_mixture(data_part, number_clusters, n_interaction, epislon)
    results <- best_mixture$results
    BIC <- best_mixture$BIC

    #choose 2 models with the smallest BIC
    smallest_value <- which.min(BIC)

    AUX_BIC <- BIC[BIC!=min(BIC)]
    second_smallest <- which(BIC== min(AUX_BIC))
    #get  the model parameters of the two smallest BIC
    best_k[1] <- number_clusters[smallest_value]
    best_k[2] <- number_clusters[second_smallest]

    #----------------------------------- Classification Train  -------------------

    j <- 1
    loss <- list()

    models <- rbind(results[[smallest_value]]$parameters, results[[second_smallest]]$parameters)
    #Compute the accuracy of predicting for the two possible models
    for (n_model in 1:nrow(models)){
      k <- 1
      model <- models[n_model,]
      n_clusters <- length(model$pk)
      parcial_results <- matrix(0,  nrow(X_train), n_clusters)
      for (i in 1:n_clusters){

        parcial_results[,k] <- model$pk[i]*dmultinorm_all(as.matrix(X_train), model$moyenne[i,] , model$variance[,,i])
        k <- k + 1
      }
      parcial_results <- parcial_results/rowSums(parcial_results)
      y_pred <- apply(parcial_results, 1, which.max)
      loss[[j]] <- loss_function(y_train, y_pred)
      j <- j + 1
    }
    loss <- rapply( loss, f=function(x) ifelse(is.na(x),0,x), how="replace" )
    #Get the models with the highest accuracy
    best_models[[n_class]] <- models[which.max(loss), ]
    choosen_k <- best_k[which.max(loss)]
    map_cluster_class[(old_k+1): (old_k + choosen_k)] <- n_class

    old_k <- old_k + choosen_k

  }

  return(list(map_cluster_class = map_cluster_class, best_models = best_models))
}
#' @title MixtureMixture_predict
#' @description  predict using the classification model
#' @param X_test : matrix of data of dimension individus x feature
#' @param n_classes : number of class to predict
#' @param map_cluster_class : maping of the cluster and its corresponding class
#' @param best_models : result of MixtureMixture_train
#' @return y_pred: prediction to each individuals
#' @export
#'
#' @examples rand<- sample(1:150, 150)
#' X_train <- iris[rand[1:100], -5]
#' y_train <- iris[rand[1:100], 5]
#' y_train <- as.numeric(factor(y_train))
#' X_test <- iris[rand[101:150], -5]
#' y_test <- iris[rand[101:150], 5]
#' y_test <- as.numeric(factor(y_test))
#' epislon <- 10^-6
#' number_clusters <- 1:3
#' n_interaction <- 50
#' results <- MixtureMixture.train (X_train, y_train, number_clusters, n_interaction, epislon)
#' n_classes <- length(unique(y_train))
#' y_pred <- MixtureMixture.predict(X_test, n_classes, results$map_cluster_class, results$best_models)
MixtureMixture.predict <- function(X_test, n_classes, map_cluster_class, best_models){

  n_clu <- length(map_cluster_class)

  parcial_results <- matrix(0,  nrow(X_test), n_clu)
  k <- 1
  for (n_class in 1:n_classes){
    model<- best_models[[n_class]]
    for (i in 1:length(model$pk)){
      parcial_results[,k] <- model$pk[i]*dmvnorm(as.matrix(X_test), model$moyenne[i,] , model$variance[,,i])
    k <- k + 1
    }
  }
  parcial_results <- parcial_results/rowSums(parcial_results)
  y_pred <- map_cluster_class[apply(parcial_results, 1, which.max)]
  return(y_pred)

}
