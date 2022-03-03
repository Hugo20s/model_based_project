## 4. Exercise d'application ####
library(caret)

#rm(list = ls())
#set.seed(3)

#SVM, Random Forest

#minist with 1.000 observations
#28x28 pixels

#separar melhor esses dados
images    <- as.matrix(read.table("minimnist/data.txt"))[1:1000, ]
labels <- as.matrix(read.table("minimnist/labels.txt", colClasses = 'integer'))[1:1000]

# Normalizing the RGB codes by dividing it to the max RGB value.
# images_normalized <- images/256

images_2 <- images[which(labels == 2), ]
cov(images_2)[500:700,500:700]

showDigit <- function(line) {
  p <- sqrt(length(line))
  mat <- matrix(as.numeric(line), p)[, p:1] # inverse column order
  image(mat, col = grey(256:1/256))
}

i <- 8; showDigit(images[i, ]); labels[i]; rm(i)


layout(matrix(1:16, ncol = 4))
par(mai = c(0, 0, 0, 0))
for (j in sample(nrow(images), 16)) showDigit(images[j, ])

df <- data.frame(images, labels)

trainIndex  <- createDataPartition(df$labels, p = .7, 
                                   list = FALSE, 
                                   times = 1)

df_train <- df[trainIndex, ]
df_test <- df[-trainIndex, ]


#c = lapply(b,image_array_resize, height = 100, width = 100) #resize
#Split into train and test
# A stratified random split of the data
# Fitting SVM to the Training set

result.models <- data.frame(matrix(ncol = 5, nrow = 0))
x <- c("Model", "Accuracy", "Precision", "Recall", "F1-score")
colnames(result.models) <- x

fitControl <- trainControl(method = "none")
result.models$Model
y_test <- as.factor(df_test$labels)
models <- c("rf", "svmLinear", "knn")

for (model in models){
  print(model)
  result <- train(as.factor(labels) ~ ., df_train, method = model , trControl = fitControl)
  y_pred <- predict(result ,newdata=df_test)

  mat <- confusionMatrix(data=y_pred,reference=y_test)
  
  
  accuracy <- mat$overall["Accuracy"]
  precision <- mean(mat[["byClass"]][ , "Precision"])
  recall <- mean(mat[["byClass"]][ , "Recall"])
  f1 <- mean(mat[["byClass"]][ , "F1"])
  
  
  result.models[nrow(result.models) + 1,] = c(model, accuracy, precision, recall, f1) 
  
}
epislon <- 10^-6

X_train <- df_train[,-ncol(df_train)]
y_train <- df_train$labels
X_test <-  df_test[,-ncol(df_test)]
y_test <- df_test$labels

number_clusters <- 1:30
n_interaction <- 40
results <- MixtureMixture.train(X_train, y_train, number_clusters, n_interaction, epislon)
print(results)

n_classes <- length(unique(y_train))
y_pred <- MixtureMixture.predict(X_test, n_classes, results$n_clusters, results$best_models)

mat <- confusionMatrix(data=as.factor(y_pred),reference=as.factor(y_test))

model <- "Mixture of Mixture"
accuracy <- mat$overall["Accuracy"]
precision <- mean(mat[["byClass"]][ , "Precision"])
recall <- mean(mat[["byClass"]][ , "Recall"])
f1 <- mean(mat[["byClass"]][ , "F1"])


result.models[nrow(result.models) + 1,] = c(model, accuracy, precision, recall, f1) 

print(result.models)
