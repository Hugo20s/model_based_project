
library("FactoMineR")

images    <- as.matrix(read.table("minimnist/data.txt"))
labels <- as.matrix(read.table("minimnist/labels.txt", colClasses = 'integer'))


images_3 <- images[which(labels == 3), ]
labels_3 <- labels[which(labels == 3)]
images_3 <- images_3/255.0

pca_3 <- PCA(images_3, scale.unit = TRUE, ncp = 40, graph = FALSE)$ind$coord

df <- data.frame(pca_3, labels_3)

trainIndex  <- createDataPartition(df$labels, p = .7, 
                                   list = FALSE, 
                                   times = 1)

df_train <- df[trainIndex, ]
df_test <- df[-trainIndex, ]

epislon <- 10^-6

X_train <- df_train[,-ncol(df_train)]
y_train <- df_train$labels
X_test <-  df_test[,-ncol(df_test)]
y_test <- df_test$labels

number_clusters <- 1:10
n_interaction <- 1000

list_ll <- numeric()
list_res <- list()
for (k in number_clusters){
  try <- 1
  #for (interaction in n_interaction){
  print("Number of Clusters")
  print(k)
  for (interaction in  1:n_interaction){
    print(interaction)
    sample(X_train)
    res <- main(X_train[sample(1:350, 250), ], 20, epislon, "random") 
    # ALGO EM   
    if (res$invalid == FALSE){
      print(res$maxll)
      list_ll[try] <- res$maxll
      list_res[[try]] <- res  
      try <- try + 1
    }
  }
}
print(list_ll)
print(list_res)
