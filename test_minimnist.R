## 4. Exercise d'application ####
library(MixtureClass)
library("caret")
library("FactoMineR")

#load data
#minist with 1.000 observations
#28x28 pixels
images    <- as.matrix(read.table("minimnist/data.txt"))
labels <- as.matrix(read.table("minimnist/labels.txt", colClasses = 'integer'))

#perform  PCA
pca <- PCA(images, scale.unit = TRUE, ncp = 40, graph = FALSE)$ind$coord

df <- data.frame(pca, labels)

trainIndex  <- createDataPartition(df$V1, p = .7, 
                                   list = FALSE, 
                                   times = 1)
#separate train and test
df_train <- df[trainIndex, ]
df_test <- df[-trainIndex, ]


# Normalizing the RGB codes by dividing it to the max RGB value.
# images_normalized <- images/255

showDigit <- function(line) {
  p <- sqrt(length(line))
  mat <- matrix(as.numeric(line), p)[, p:1] # inverse column order
  image(mat, col = grey(256:1/256))
}

i <- 8; showDigit(images[i, ]); labels[i]; rm(i)


layout(matrix(1:16, ncol = 4))
par(mai = c(0, 0, 0, 0))
for (j in sample(nrow(images), 16)) showDigit(images[j, ])

# Fitting models

result.models <- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("Model", "Accuracy", "Precision", "Recall", "F1-score")
colnames(result.models) <- x

fitControl <- trainControl(method = "none")
result.models$Model
y_test <- as.factor(df_test$V1)
models <- c("rf", "svmLinear", "knn")

for (model in models){
  print(model)
  result <- train(as.factor(V1) ~ ., df_train, method = model , trControl = fitControl)
  y_pred <- predict(result ,newdata=df_test)

  mat <- confusionMatrix(data=y_pred,reference=y_test)


  accuracy <- mat$overall["Accuracy"]
  precision <- mean(mat[["byClass"]][ , "Precision"])
  recall <- mean(mat[["byClass"]][ , "Recall"])
  f1 <- mean(mat[["byClass"]][ , "F1"])


  result.models[nrow(result.models) + 1,] = c(model, accuracy, precision, recall, f1)

}
print(result.models)

epislon <- 10^-6

X_train <- df_train[,-ncol(df_train)]
y_train <- df_train$labels
X_test <-  df_test[,-ncol(df_test)]
y_test <- df_test$labels

number_clusters <- 1:30
n_interaction <- 40

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


model <- "Mixture of Mixture"
accuracy <- mat$overall["Accuracy"]
precision <- mean(mat[["byClass"]][ , "Precision"])
recall <- mean(mat[["byClass"]][ , "Recall"])
f1 <- mean(mat[["byClass"]][ , "F1"])


result.models[nrow(result.models) + 1,] = c(model, accuracy, precision, recall, f1) 

print(result.models)
