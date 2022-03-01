library(mvtnorm)
Estep <- function(data,mix.pi,mu,sigma,K){
  result <- apply(data, 1, function(xt){
    gammar_i <- sapply(1:K,function(k) {
      mix.pi[k] * dmvnorm(xt, mu[,k], sigma[,,k])
    });
    gammar_i / sum(gammar_i);
  });
  t(result);
}

Mstep <- function(gamma,data,K,N,degree) {
  gamma.sum <- apply(gamma,2,sum)
  new.mix <- gamma.sum/N;
  new.mu <- t(t(t(data) %*% gamma) / gamma.sum);
  new.sigma <- array(0, dim=c(degree, degree, K));
  for(k in 1:K) {
    sig <- matrix(0,degree,degree);
    for(n in 1:N) {
      sig <- sig + gamma[n, k] * (data[n,] %*% t(data[n,]));
    }
    new.sigma[,,k] <- sig / gamma.sum[k] - new.mu[,k] %*% t(new.mu[,k])
  }  
  list(new.mix, new.mu, new.sigma);
}

set.seed(1010);

K <- 3;  #number of the clusters
degree <- 2;  #degree of the data(2D data this time)
N <- 1000; #number of samples

mix.pi <- runif(K);
mix.pi <- runif(K)/sum(mix.pi);

mu <- 0.1*matrix(runif(K*degree), nrow=degree, ncol=K);

sigma <- array(0, dim=c(degree, degree, K));
for(k in 1:K){
  sigma[,,k] <- 0.1*diag(runif(degree));
}

data <-  t(apply(rmultinom(N,1,mix.pi),2,function(num) {
  maxindex <- which.max(num)
  rmvnorm(1, mu[,maxindex], sigma[,,maxindex])
}));

n.trial <- 1000

count <- 0
errorlist <-  rep(0,n.trial)
for(i in 1:n.trial){
  count <- count +1
  gamma <- Estep(data,mix.pi,mu,sigma,K);
  result <- Mstep(gamma,data,K,N,degree);
  new.mix.pi <- result[[1]]; new.mu <- result[[2]]; new.sigma <- result[[3]];
  error <-sum((new.mix.pi-mix.pi)^2) + sum((new.mu-mu)^2) + sum((new.sigma-sigma)^2)
  errorlist[i] <- error
  mix.pi <- new.mix.pi;
  mu <- new.mu;
  sigma <- new.sigma;
  if (error < 1e-8) break;
}

plot(errorlist[1:count],log="y")
