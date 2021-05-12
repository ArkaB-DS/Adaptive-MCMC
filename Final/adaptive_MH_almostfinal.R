set.seed(1) # for reproducibility
library(mvtnorm) # to deal with multivariate normal distribution
library(tictoc)

d <- 20 # dimension
N <- 2e5 # length of MC
X <- matrix(10,nrow=N+1,ncol=d)

# generating a Sigma for the target distributions
M <- matrix(rnorm(d*d),nrow=d)
Sigma_inv <- solve(M%*%t(M))

beta <- .05
p=0 # for calculating acceptance probability

# implementing the adaptive MH
tic()
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
toc()
cat("The acceptance probability is: ",p/N,"\n")
pdf("Final_Plots.pdf")
for (i in 1:d) ts.plot(X[,i],ylab=i)
for (i in 1:d) acf(X[,i],ylab=i,col="blue",main="")
dev.off()
R=range(Sigma_n-M%*%t(M))
prob=p/N
FNorm=sum((Sigma_n-M%*%t(M))^2)
save(X,prob,R,FNorm,file="AM_d10_1e5_start10.Rdata")
