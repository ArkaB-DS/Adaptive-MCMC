set.seed(1)
library(mvtnorm)
library(clusterGeneration)

d <- 4 # dimension
N <- 5e4 # length of MC
X <- matrix(10,nrow=N+1,ncol=d)


# generating a Sigma for the target distribution
eig <- seq(1,1e1,length=d) 
Sigma <- genPositiveDefMat(dim=d,covMethod="eigen",eigenvalue=eig)$Sigma 

# target - N_10(mean=0, var= Sigma)

beta <- .05
p=0 # for calculating acceotance probability

# implementing the adaptive MH
for(i in 2:(N+1))
{
 print(i)
 Sigma_n<-cov(X)
 if(rbinom(n=1,size=1,prob=beta)){ 
 Y<- rmvnorm(n=1,mean=X[i-1,],sigma=0.1^2/d*diag(d))
 } else {
 Y<-rmvnorm(n=1,mean=X[i-1,],sigma=2.38^2/d*Sigma_n)  }
 Y <- as.vector(Y)
 U <- runif(1)
 r <- dmvnorm(Y,mean=rep(0,d),sigma=Sigma)/dmvnorm(X[i-1,],mean=rep(0,d),sigma=Sigma)  
 if (U<r)
 {
  X[i,] <- Y
  p <- p+1
 } else{
  X[i,] <- X[i-1,]
       }
}

p/N
par(mfrow=c(2,2),bg="pink") 
for (i in 1:4 )ts.plot(X[1:3e4,i],ylab=i,col="brown")
for (i in 1:4 )acf(X[,i],ylab=i,col="blue",main="")
#for (i in 1:4 )plot(density(X[,i]),ylab=i,col="blue",main="")
