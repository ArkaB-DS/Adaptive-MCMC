#########################################################################
#####
#####  Implementing Adaptive MH Algorithm from Roberts & Rosenthal, 2009
#####		
#########################################################################

set.seed(1) 	 # for reproducibility
library(mvtnorm) # for multivariate normal 
library(tictoc)  # to keep the track of time



###################    TARGET
###    Target is a “banana-shaped” distribution, as proposed by 
###    Haario et al. (1999, 2001)

d <- 20          # dimension
B <- 0.1		 # rotation parameter



###################    SAMPLER

N <- 2e5    			# number of samples
p <- 0  				# records acceptances
beta <- 0.05						


### Declaration and (bad) Initialization of the chain

X <- matrix(0, nrow = N, ncol = d)
X[1,] <- c(20, 30, rep(2.5, d-2)) 


### Sampling begins 

tic()

for(i in 2:N){
	
  	if(i%%10000 == 0) print(i)   # For feedback

  	# Proposal step
  	if( i <= 2*d){
		Y <- rnorm(d, mean = X[i-1, ], sd = 0.1/sqrt(d))
	}
	else{
		if(rbinom(1, 1, beta)){
			Y <- rnorm(d, mean = X[i-1, ], sd = 0.1/sqrt(d))
		}
		else{
			S <- cov(X[1:i-1, ])
			Y <- rmvnorm(1, mean = X[i-1, ], sigma = 2.38^2/d*S)
		}
	}

	# Accept - Reject Step
	r <- -(Y[1]^2)/200 - ((Y[2] + B*(Y[1]^2) - 100*B)^2)/2 - sum(Y[-c(1, 2)]^2)/2 +
				(X[i-1,1]^2)/200 + ((X[i-1,2] + B*(X[i-1,1]^2) - 100*B)^2)/2 + sum(X[i-1,-c(1, 2)]^2)/2

	if (log(runif(1)) < r){
		 X[i,] <- Y
  		 p <- p+1

	}else{
		X[i,] <- X[i-1,]
	}
}

toc()

# Saving the plots

#pdf("plots_banana_final.pdf")
#for(k in 1:d)  ts.plot(X[, k], type = "l", main = "Trace plot for AM")
#dev.off()

prob <- p/N
save(X, prob, file = "banana_final.Rdata")
