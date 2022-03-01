#Test generate data
clus_1 <- matrix(rnorm(20*4, 0, 1), 20, 4)
clus_2 <- matrix(rnorm(20*4, 5, 2), 20, 4)
clus_3 <- matrix(rnorm(20*4, 10, 2), 20, 4)

X <- rbind(clus_1, clus_2, clus_3)
y <- c(rep(1,20), rep(2,20), rep(3,20))



#Test iris 
K <- 3
data <- iris[sample(1:150, 100),]
X <- data[,-5]
y <- data[,5]


#Execution
debug(clustering)
epislon <- 10^-6
log_likelihood <- clustering(X, K, epislon)


#test error 
y2 <- as.integer(y)
error <- sum((y2 - res$y_pred)^2) / nrow(X)



#test strcitement croissant
l <- res$log_likelihood
for (i in seq(length(l)-1)){
  if (l[i] < l[i]+1) {
    print('ok')
  } else {
    print('not ok')
  }
}