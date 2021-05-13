#########################################################################
#####
#####  Implementing Adaptive MH Algorithm from Roberts & Rosenthal, 2009
#####		
#########################################################################

set.seed(1) 	 # for reproducibility
library(mvtnorm) # for multivariate normal 
library(tictoc)  # to keep the track of time



###################    TARGET
###    Target is a d-dimensional Gaussian with mean 0 and 
###    variance-covariance matrix MM'. Entries of M are 
###    iid N(0, 1) variates

d <- 20             					# dimension
M <- matrix(rnorm(d*d),nrow=d)			# M matrix for target covariance
Sigma_inv <- solve(M%*%t(M))			# Inverse of target covariance



###################    SAMPLER

N <- 2e5    			# number of samples
p <- 0  				# records acceptances
beta <- 0.05	

### Declaration and (bad) Initialization of the chain
X <- matrix(10, nrow=N+1, ncol=d)


### Sampling begins 

tic()

for(i in 2:(N+1)){
 	
 	if (i%%1e4==0) print(i) 	# For feedback
	
	# Proposal Step
	if (i<=2*d){
	     Y<- rnorm( n = d , mean = X[i-1,] , sd = 0.1/sqrt(d) )
	} 
	else{    
	 	if( rbinom( n = 1 , size = 1 , prob = beta) ){ 
	 		Y<- rnorm( n = d , mean = X[i-1,] , sd = 0.1/sqrt(d) )
	 	} 
	 	else{
	 		Sigma_n <- cov(X[1:i,])
	 		Y <- rmvnorm( n = 1 , mean = X[i-1,] , sigma = 2.38^2/d*Sigma_n )  
	    }
    }
	 
	# Accept - Reject Step
	r <-  (crossprod(X[i-1, ], Sigma_inv)%*%X[i-1, ] - crossprod(Y, Sigma_inv)%*%Y)/2
	
	if ( log(runif(1)) < r ){
	  	X[i,] <- Y
	  	p <- p+1
	} else{
	  	X[i,] <- X[i-1,]
    }
}

toc()

# Saving the plots

#pdf("am_final.pdf")
#for (i in 1:d) ts.plot(X[,i],ylab=i)
#for (i in 1:d) acf(X[,i],ylab=i,col="blue",main="")
#dev.off()

# Comparing Sample covariance with Target covariance

R <- range(Sigma_n-M%*%t(M))
prob <- p/N
FNorm <- sum((Sigma_n-M%*%t(M))^2)

save(X,prob,R,FNorm, file="am_final.Rdata")
