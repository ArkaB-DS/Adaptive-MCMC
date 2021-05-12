set.seed(1) # for reproducibility
library(mvtnorm) # to deal with multivariate normal distribution
library(mvnfast)

d <- 10 # dimension
N <- 1e5 # length of MC
X <- matrix(0,nrow=N+1,ncol=d)


# generating a Sigma for the target distribution
M <- matrix(rnorm(d*d),nrow=d)
Sigma_inv <- solve(M%*%t(M))

beta <- .05
p=0 # for calculating acceptance probability

# implementing the adaptive MH
for(i in 2:(N+1))
{
 if (i%%1e4==0) print(i)
 # proposing
 if (i<=2*d) 
     {
     Y<- rnorm( n = d , mean = X[i-1,] , sd = 0.1/sqrt(d) )
     } else {    
 if( rbinom( n = 1 , size = 1 , prob = beta) ){ 
 Y<- rnorm( n = d , mean = X[i-1,] , sd = 0.1/sqrt(d) )
 } else {
 Sigma_n <- cov(X[1:i,])
 Y<- rmvnorm( n = 1 , mean = X[i-1,] , sigma = 2.38^2/d*Sigma_n )  
        }
            }
 Y <- as.vector(Y)
 r <-  (crossprod(X[i-1, ], Sigma_inv)%*%X[i-1, ] - crossprod(Y, Sigma_inv)%*%Y)/2
 if ( log(runif(1)) < r )
 {
  X[i,] <- Y
  p <- p+1
 } else{
  X[i,] <- X[i-1,]
       }
}
cat("The acceptance probability is: ",p/N,"\n")
#par(mfrow=c(2,1),bg="pink") 
ts.plot(X[1:1e4,1],ylab=i,col="brown")
# acf(X[,3],ylab=i,col="blue",main="")
