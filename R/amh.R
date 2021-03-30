set.seed(1)
library(mvtnorm)
library(clusterGeneration)

d <- 4 # dimension
N <- 5e5 # length of MC
X <- matrix(10, nrow = N, ncol=d)


# generating a Sigma for the target distribution
eig <- seq(0.01, 100, length = d) 
Sigma <- genPositiveDefMat(dim=d, covMethod="eigen", eigenvalue=eig)$Sigma 
E <- solve(Sigma)


beta <- .05
p = 0 # for calculating acceotance probability
n0 = 0

# implementing the adaptive MH

for(i in 2:N){
	# Feedback
	if(i%%10000 == 0) print(i)

	# Proposing the next step. For first few steps, run a simple MH sampler. 
	# After that, run the adaptive version
	if(p <= d){
		Y <- rnorm(n = d, mean=X[i-1,], sd = 2.38/sqrt(d))
		n0 = i
	}else{
		if(rbinom(n = 1, size = 1, prob = beta)){ 
			Y <- rnorm(n = d, mean=X[i-1,], sd = sqrt(0.1/d))
		}else{
			Sigma_n <- cov(X[1:i-1, ])
			Y <- rmvnorm(n = 1, mean=X[i-1,], sigma = 2.38^2/d*Sigma_n)  
		}
	}
  
	# The accept-reject step
	Y <- as.vector(Y)
	r <- (crossprod(X[i-1, ], E)%*%X[i-1, ] - crossprod(Y, E)%*%Y)/2
    
	if (log(runif(1)) < r){
		X[i,] <- Y
		p <- p + 1
	}else{
		X[i,] <- X[i-1,]
	}
}

p/(N-1)
par(mfrow=c(2,2),bg="pink") 
for (i in 1:4 )ts.plot(X[,i],ylab=i,col="brown")
for (i in 1:4 )acf(X[,i],ylab=i,col="blue",main="")