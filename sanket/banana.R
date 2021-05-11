#########################################################################
#####
#####  Implementing Adaptive MH Algorithm from Roberts & Rosenthal, 2009
#####		
#########################################################################


set.seed(1) # for reproducibility
#library(mvtnorm) # for multivariate normal 


###################    TARGET
###    Target is a “banana-shaped” distribution, as proposed by 
###    Haario et al. (1999, 2001)

d <- 20
B <- 0.1

logf <- function(x){
	z <- -(x[1]^2)/200 - ((x[2] + B*(x[1]^2) - 100*B)^2)/2 - sum(x[-c(1, 2)]^2)/2
	return(z)
}



###################    PROPOSAL
###     Proposal is a mixture of Gaussian 

propose <- function(n, X){
	# n is the number of current step
	# X the generated process
	
	# For n <= 2d, sample from a fixed kernel. For n > 2d, sample from a mixture kernel.

	if( n <= 2*d){
		Y <- rnorm(d, mean = X[n-1, ], sd = 0.1/sqrt(d))
	}
	else{
		if(rbinom(1, 1, beta)){
			Y <- rnorm(d, mean = X[n-1, ], sd = 0.1/sqrt(d))
		}
		else{
			S <- cov(X[1:n-1, ])
			Z <- rnorm(d, mean = X[n-1, ], sd = 2.38/sqrt(d))
			Y <- chol(S)%*%Z
		}
	}
	return(as.vector(Y))
}


###################    ACCEPTANCE RULE

alpha <- function(x, y){
	r <- logf(y) - logf(x)
	
	if (log(runif(1)) < r){
		return(c(y, 1))
	}else{
		return(c(x, 0))
	}
}
 
 
###################    SAMPLER

N <- 1e5    						# number of samples
beta <- 0.05
p <- 0  							# records acceptances

X <- matrix(0, nrow = N, ncol = d)
X[1, ] <- rnorm(d, mean = 0, sd = 1)
current <- X[1, ]


for(i in 2:N){
	
  if(i%%10000 == 0) print(i)
  
	prop <- propose(i, X)
	move <- alpha(X[i-1, ], prop)
	
	X[i, ] <- move[1:d]
	p <- p + move[d+1]
}

par(mfrow = c(1, 1))
pdf("plots_banana.pdf")
for(k in 1:d){
  ts.plot(X[, k], type = "l", col = "blue", main = "Trace plot for AM")
}
dev.off()

p/N   #0.0479


