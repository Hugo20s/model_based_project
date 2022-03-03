

library("MLmetrics")

MixtureMixture.train <- function(X_train, y_train, number_clusters, n_interaction, epislon){

  classes <- unique(y_train)
  best_models <- list()
  n_clusters <- array(NA)
  old_k <- 0
  for (n_class in seq(classes)){

    data_part <- X_train[which(y_train == classes[n_class]),]
    
    n <- nrow(data_part)
    p <- ncol(data_part)
    BIC <- numeric(0)
    results <- list()
    i <- 1
    
    for (k in number_clusters){
      list_ll <- numeric()
      list_res <- list()
      try <- 1
      #for (interaction in n_interaction){
      print("Number of Clusters")
      print(k)
      for (interaction in  1:n_interaction){
        print(interaction)
        res <- main(data_part, k, epislon, "random") 
        # ALGO EM   
        if (res$invalid == FALSE){
          print(res$maxll)
          list_ll[try] <- res$maxll
          list_res[[try]] <- res  
          try <- try + 1
        }
      }
      if (!(length(list_res) == 0) ){
        res <- list_res[[which.max(list_ll)]] 
        print("-----Max Log-likelihood-------")
        print(res$maxll)
        print("-----Cluster Number------- "); print(k)
        print("Class "); print(n_class)
        
        BIC[i] <- -2*res$maxll + (k*(p*(p + 1)/2 + p + 1) - 1)*log(n)
        
        results[[i]] <- res 
        i <- i + 1
        }
        
    }
    print("----------BIC-----------")
    print(BIC)
    plot(number_clusters, BIC)
    #choose the smallest bic
    best_k <- number_clusters[which.min(BIC)]
    
    best_models[[n_class]] <- results[[which.min(BIC)]]$parameters
    
    n_clusters[(old_k+1): (old_k + best_k)] <- n_class 
    
    old_k <- old_k + best_k 
  }
  
  return(list(n_clusters = n_clusters, best_models = best_models))
}
MixtureMixture.predict <- function(X_test, n_classes, n_clusters, best_models){

  n_clu <- length(n_clusters)
  
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

  y_pred <- n_clusters[apply(parcial_results, 1, which.max)] 
  return(y_pred)

}

# rand<- sample(1:150, 150)
# 
# X_train <- iris[rand[1:100], -5]
# y_train <- iris[rand[1:100], 5]
# y_train <- as.numeric(factor(y_train))
# X_test <- iris[rand[101:150], -5]
# y_test <- iris[rand[101:150], 5]
# y_test <- as.numeric(factor(y_test))
# epislon <- 10^-6
# 
# number_clusters <- 1:5
# n_interaction <- 50
# results <- MixtureMixture.train(X_train, y_train, number_clusters, n_interaction, epislon)
# print(results)
# 
# n_classes <- length(unique(y_train))
# y_pred <- MixtureMixture.predict(X_test, n_classes, results$n_clusters, results$best_models)
# 
# mat <- confusionMatrix(data=as.factor(y_pred),reference=as.factor(y_test))

