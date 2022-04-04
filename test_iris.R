
library(MixtureClass)
library(caret)
rand<- sample(1:150, 150)

X_train <- iris[rand[1:100], -5]
y_train <- iris[rand[1:100], 5]
y_train <- as.numeric(factor(y_train))
X_test <- iris[rand[101:150], -5]
y_test <- iris[rand[101:150], 5]
y_test <- as.numeric(factor(y_test))
epislon <- 10^-6

number_clusters <- 1:3
n_interaction <- 50
results <- MixtureClass::MixtureMixture.train(X_train, y_train, number_clusters, n_interaction, epislon)
print(results$map_cluster_class)

print(results$best_models)# 
n_classes <- length(unique(y_train))
y_pred <- MixtureClass::MixtureMixture.predict(X_test, n_classes, results$map_cluster_class, results$best_models)
print("Y_pred")
print(y_pred)
print("Y_test")
print(y_test)
mat <- confusionMatrix(data=as.factor(y_pred),reference=as.factor(y_test))
print("Confusion Matrix")
print(mat)
