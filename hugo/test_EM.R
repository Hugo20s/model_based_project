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



a <- rbind(c(0.8424933,  0.3240539 , 2.4627042,  1.3610497),
           c(  0.3240539,0.1246449, 0.9472499, 0.5235202),
          c(2.4627042,  0.9472499,   7.1987719, 3.9785222),
          c( 1.3610497, 0.5235202, 3.9785222, 2.1988365))

qr.solve(a)

[,1]      [,2]      [,3]      [,4]
[1,] 0.8424933 0.3240539 2.4627042 1.3610497
[2,] 0.3240539 0.1246449 0.9472499 0.5235202
[3,] 2.4627042 0.9472499 7.1987719 3.9785222
[4,] 1.3610497 0.5235202 3.9785222 2.1988365


#log_k <- numeric(nrow(data))
# for (j in 1:nrow(data)){
#   # for (k in 1:K){
#   #   log_k[j] <- lof(logSumExp(
#   # }
#   
#   log_k[j] <- logSumExp(lapply(1:K, function(k) {
#     pk_modif <- replace(pk[k], pk[k]==0, 10^-8)
#     log(pk_modif) + dmvnorm(data[j,], moyenne[k,], variance[,,k], log = TRUE)
#   }))
#   # print("log")
#   pk_modif <- replace(pk[k], pk[k]==0, 10^-8)
#   # print(log(pk_modif))
#   # print("densité")
#   #print(dmvnorm(data[j,], moyenne[k,], variance[,,k], log = TRUE))
#   log(pk_modif) + dmvnorm(data[j,], moyenne[k,], variance[,,k], log = TRUE)
# }))
#}
#log_likelihood[i] <- sum(log_k)


# for(k in 1:K) {
#   temp <- 0
#   for(n in 1:N) {
#     temp <- temp + tk[n, k] * t(data[n,] - moyenne[k,]) * (data[n,] - moyenne[k,])}
#   print(temp)
#   variance[,,k] <- temp/nk[k]
# }


# for(k in 1:K) {
#   temp <- replace(tk[, k], tk[, k]==0, 10^-8)
#   variance[,,k] = cov.wt(data, wt = temp , cor = TRUE, method = "ML")$cov
#   print(paste("det", k))
#   print(det(variance[,,k]))
# }