#########################################################################
#####
#####  Implementing Adaptive MH Algorithm from Roberts & Rosenthal, 2009
#####		
#########################################################################


set.seed(1) # for reproducibility
#library(mvtnorm) # for multivariate normal 

#"%^%" <- function(x, n) 
#  with(eigen(x), vectors %*% (values^n * t(vectors)))


###################    TARGET
###    Target is a d-dimensional Gaussian distribution N(0, MM') 
###    where entries of M are iid N(0, 1)

d <- 10


M <- matrix(rnorm(d^2, 0, 1), nrow = d)
E <- M%*%t(M)
E_inv <- solve(M%*%t(M))



###################    PROPOSAL
###     Proposal is a mixture of Gaussian 

propose <- function(p, x, S){
	# n is the number of current step
	# X the generated process
	
	# For n <= 2d, sample from a fixed kernel. For n > 2d, sample from a mixture kernel.

	if( p <= d){
		Y <- rnorm(d, mean = x, sd = 0.1/sqrt(d))
	}
	else{
		if(rbinom(1, 1, beta)){
			Y <- rnorm(d, mean = x, sd = 0.1/sqrt(d))
		}
		else{
			Z <- rnorm(d, mean = x, sd = 2.38/sqrt(d))
			Y <- chol(S)%*%Z
		}
	}
	return(as.vector(Y))
}


###################    ACCEPTANCE RULE

alpha <- function(x, y){
	r <- crossprod(x, E_inv)%*%x - crossprod(y, E_inv)%*%y
	
	if (log(runif(1)) < r){
		return(c(y, 1))
	}else{
		return(c(x, 0))
	}
}
 
 
###################    SAMPLER

N <- 1e6    						# number of samples
beta <- 0.05
p <- 0  							# records acceptances

X <- matrix(0, nrow = N, ncol = d)
X[1, ] <- rep(0, d)
#X[1, ] <- rnorm(d, mean = 0, sd = diag(M%*%t(M)))
current <- X[1, ]
S <- diag(rep(1, d))
#lmabda <- matrix(0, nrow = N, ncol = d)


for(i in 2:N){
	
  if(i%%10000 == 0) print(i)
  
	prop <- propose(p, X[i-1, ], S)
	move <- alpha(X[i-1, ], prop)
	
	X[i, ] <- move[1:d]
	p <- p + move[d+1]
	
	if(p > d) S <- cov(X[1:i, ])
#	lambda[i, ] <- eigen((S %^% (0.5)) %*% (E %^% (-0.5)))$values
}

par(mfrow = c(1, 1))
pdf("plots_adaptiveMH.pdf")
for(k in 1:d){
  ts.plot(X[, k], type = "l", col = "black", main = "Trace plot for AM")
}
dev.off()

p/N   #0.11487


