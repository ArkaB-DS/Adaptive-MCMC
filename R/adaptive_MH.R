### R Code for Rosenthal's Adaptive MH
set.seed(1)
library(mvtnorm)
library(clusterGeneration)

d <- 5 # dimension
N <- 5e4 # length of MC
X <- matrix(0,nrow=N+1,ncol=d)


# generating a Sigma for the target distribution
eig <- seq(1e-2,1e3,length=d) 
Sigma <- genPositiveDefMat(dim=d,covMethod="eigen",eigenvalue=eig)$Sigma 

# target - N_10(mean=0, var= Sigma)

beta <- .05
p=0 # for calculating acceotance probability

# implementing the adaptive MH
for(i in 2:(N+1))
{
 print(i)
#Sigma_n<- 1/(N+1-1)*tcrossprod(t(X)-apply(X,2,mean),t(X)-apply(X,2,mean))
Sigma_n<-cov(X)
# Y <-  (1-beta)*rmvnorm(n=1,mean=X[i-1,],sigma=2.38^2/d*Sigma_n)+ beta*rmvnorm(n=1,mean=X[i-1,],sigma=0.1^2/d*diag(d))
 if(rbinom(n=1,size=1,prob=beta)){ 
 Y<- rmvnorm(n=1,mean=X[i-1,],sigma=0.1^2/d*diag(d))
 } else {
 Y<-rmvnorm(n=1,mean=X[i-1,],sigma=2.38^2/d*Sigma_n)  }
 Y <- as.vector(Y)
 if (runif(1)<dmvnorm(Y,mean=rep(0,d),sigma=Sigma)/dmvnorm(X[i-1,],mean=rep(0,d),sigma=Sigma))
 {
  X[i,] <- Y
  p <- p+1
 } else{
  X[i,] <- X[i-1,]
       }
}

p/N # acceptance probabaility ; =0.851
# the code was run for i=34391
ts.plot(X[1:2e4,1])
acf(X[,1])
plot(density(X[,5]))
