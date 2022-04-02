

library("MLmetrics")
library("caret")
loss_function <- function(y_true, y_pred){

    # calculate the metric
    loss <- Accuracy(y_true ,y_pred)
    
    # convert to tensor
    return(loss)
}

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
    #for (interaction in n_interaction){
    print("Number of Clusters")
    print(k)
    for (interaction in  1:n_interaction){
      res <- main(data_part, k, epislon, "random") 
      # ALGO EM   
      if (res$invalid == FALSE){
        #print(res$maxll)
        list_ll[try] <- res$maxll
        list_res[[try]] <- res  
        try <- try + 1
      }
    }
    print(list_ll)
    if (!(length(list_res) == 0) ){
      res <- list_res[[which.max(list_ll)]] 
      print("-----Max Log-likelihood-------")
      #print(res$maxll)
      print("-----Cluster Number------- "); print(k)
      
      BIC[i] <- -2*res$maxll + (k*(p*(p + 1)/2 + p + 1) - 1)*log(n)
      
      results[[i]] <- res 
      i <- i + 1
    }
    
  }
  return (list(BIC = BIC, results = results))
  
}
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
    
    print("BIC")
    print(BIC)

    smallest_value <- which.min(BIC)
    print(smallest_value)
    
    AUX_BIC <- BIC[BIC!=min(BIC)]
    print(AUX_BIC)
    second_smallest <- which(BIC== min(AUX_BIC))
    print(second_smallest)
    best_k[1] <- number_clusters[smallest_value]
    best_k[2] <- number_clusters[second_smallest]
    
    print("Best K")
    print(best_k)
    #----------------------------------- Classification Train  -------------------
    
    j <- 1
    loss <- list()
    
    models <- rbind(results[[smallest_value]]$parameters, results[[second_smallest]]$parameters)
    
    for (n_model in 1:nrow(models)){
      k <- 1
      model <- models[n_model,]
      n_clusters <- length(model$pk)
      parcial_results <- matrix(0,  nrow(X_train), n_clusters)
      for (i in 1:n_clusters){
        
        print(model$pk[i])
        print(model$moyenne[i,])
        print(model$variance[,,i])
        print(dmultinorm_all(as.matrix(X_train), model$moyenne[i,] , model$variance[,,i]))
        parcial_results[,k] <- model$pk[i]*dmultinorm_all(as.matrix(X_train), model$moyenne[i,] , model$variance[,,i])
        k <- k + 1
      }
      print(parcial_results)
      parcial_results <- parcial_results/rowSums(parcial_results)
      print(parcial_results)
      y_pred <- apply(parcial_results, 1, which.max)
      
      loss[[j]] <- loss_function(y_train, y_pred)
      j <- j + 1
    }
    print("Accuracy")
    print(loss)
    best_models[[n_class]] <- models[which.max(loss), ]
    choosen_k <- best_k[which.max(loss)]
    print("Number of class")
    print(n_class)
    map_cluster_class[(old_k+1): (old_k + choosen_k)] <- n_class 
    
    old_k <- old_k + choosen_k 
    
  }
  
  return(list(map_cluster_class = map_cluster_class, best_models = best_models))
}
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
  print(parcial_results)
  print(apply(parcial_results, 1, which.max))
  y_pred <- map_cluster_class[apply(parcial_results, 1, which.max)] 
  print(y_pred)
  return(y_pred)

}

rand<- sample(1:150, 150)

X_train <- iris[rand[1:100], -5]
y_train <- iris[rand[1:100], 5]
y_train <- as.numeric(factor(y_train))
X_test <- iris[rand[101:150], -5]
y_test <- iris[rand[101:150], 5]
y_test <- as.numeric(factor(y_test))
epislon <- 10^-6

number_clusters <- 1:5
n_interaction <- 50
results <- MixtureMixture.train(X_train, y_train, number_clusters, n_interaction, epislon)
print(results)

n_classes <- length(unique(y_train))
y_pred <- MixtureMixture.predict(X_test, n_classes, results$map_cluster_class, results$best_models)

mat <- confusionMatrix(data=as.factor(y_pred),reference=as.factor(y_test))

