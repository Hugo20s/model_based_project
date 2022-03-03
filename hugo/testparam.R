K <- 3
f <- 2
N <- 4
#n = 3
#f = 2
#k = 3


d <- cbind(c(1,2,3,4), c(4,5,6,7))    #n;f
t <- cbind(c(0.1,0.2,0.3, 0.6), c(0.7,0.7,0.4, 0.1), c(0.2, 0.1, 0.3, 0.2))         #n;k
st <- colSums(t) #k
moyenne <- matrix(0, K, f)
variance <- array(0, c(f,f, K))
for(k in 1:K) {
  moyenne[k,] <- (1/sum(t[,k])) * colSums(apply(d,2, function(x) x*t[, k])  )
  variance[,,k] <- cov.wt(d, wt = t[,k], method="ML")$cov
}

dens_emd <- list()
dens_all<- list()
dens_mv<- list()
for (k in 1:K) {
  dens_emd[[k]] <- emdbook::dmvnorm(d, moyenne[k,], variance[,,k] )  
  #dens_all[[k]] <- dmultinorm_all(data, moyenne[k,], variance[,,k], f)
  dens_mv[[k]] <- mvtnorm::dmvnorm(d, moyenne[k,], variance[,,k] )  
}

M = matrix(c(1,0.5,0.5,0.5,1,0.5,0.5,0.5,1),nrow=3)
mu <- c(1,2,3)
x <- matrix(1:6,nrow=2)
res1 <- emdbook::dmvnorm(x, mu, M)
res2 <- mvtnorm::dmvnorm(x , mu, M)
res3 <- dmultinorm_all(x, mu, M)


apply(M, 1, which.max)


log(res3)

i <- 1
sweep(x[i,], 2, mu, FUN = '-', check.margin=FALSE)

(x[i,] - mu) * 2


param_sauv <- parameters
tk_sauv <- tk 
list_densite <- list()
list_densite_all <- list()
for (k in 1:K){
  list_densite[[k]] <- mvtnorm::dmvnorm(data, moyenne[k,], variance[,,k])
  #list_densite_all[[k]] <- dmultinorm_all(data, moyenne[k,], variance[,,k], p)
}
# for (k in 1:K){
#   print("emd")
#   print(emdbook::dmvnorm(data, moyenne[k,], variance[,,k]))
#   print("all")
#   #print(dmultinorm_all(data, moyenne[k,], variance[,,k], p))
#   print("mv")
#   print(mvtnorm::dmvnorm(data, moyenne[k,], variance[,,k]))
# }
i <- i + 1
#return()
