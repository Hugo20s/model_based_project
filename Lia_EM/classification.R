

library("MLmetrics")

MixtureMixture.train <- function(X_train, y_train, epislon){

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
    
    #LOG-LIKE KO
    #sometime it doesn't converge (log-likelihood smaller then epislon )

    #Error in if (LL < log_likelihood[i - 1]) { : 
    #missing value where TRUE/FALSE needed missing value where TRUE/FALSE needed
    
    number_clusters <- c(1, 2,3,4,5)
    for (k in number_clusters){
      res <- main(data_part, k, epislon, "small")
      print("-----Max Log-likelihood-------")
      print(res$maxll)
      print("-----Cluster Number------- "); print(k)
      print("Class "); print(n_class)
      
      BIC[i] <- -2*res$maxll + (k*(p*(p + 1)/2 + p + 1) - 1)*log(n)
      
      results[[i]] <- res 
      i <- i + 1
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
    print(length(model$pk))
    for (i in 1:length(model$pk)){
      parcial_results[,k] <- model$pk[i]*dmvnorm(as.matrix(X_test), model$moyenne[i,] , model$variance[,,i])
    k <- k + 1
    }
  }
  parcial_results <- parcial_results/rowSums(parcial_results)

  print(parcial_results)
  y_pred <- n_clusters[apply(parcial_results, 1, which.max)] 
  return(y_pred)

}


X_train <- iris[sample(1:150, 100), -5]
y_train <- iris[sample(1:150, 100), 5]
y_train <- as.numeric(factor(y_train))
X_test <- iris[sample(1:150, 50), -5]
y_test <- iris[sample(1:150, 50), 5]
y_test <- as.numeric(factor(y_test))
epislon <- 10^-6

n_classes <- length(unique(y_train))

results <- MixtureMixture.train(X_train, y_train, epislon)
print(results)

y_pred <- MixtureMixture.predict(X_test, n_classes, results$n_clusters, results$best_models)
table(y_pred, y_test)
Accuracy(y_pred, y_test)

# # 
# # 
# # plot(res$log_likelihood)  
# # plot(res$y_pred)
# # 
# # c <- 0
# # for (i in 1:100){
# #   res <- clustering(data, K, epsilon, "small")
# #   if (typeof(res) == "list"){ 
# #     # print(res$log_likelihood[length(res$log_likelihood)])
# #   }
# #   # print("-----")
# # }
# # 
# # table(res$y_pred, y)
# 
# # 
# res <- main(data, K, epislon)
# plot(res$log_likelihood)
# plot(res$y_pred)
# tkres<- replicate(20, main(data, K, epislon))
# 
# result <- do.call(rbind, tkres[1,])
# best_iteration <- which.max(result[,ncol(result)])
# 
# y_pred <- tkres[2, best_iteration]$y_pred
# 
# plot(1:length(as.list(y_pred)), y_pred)
# plot(tkres[1,1]$log_likelihood)
