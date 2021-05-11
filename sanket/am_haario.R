#########################################################################
#####
#####  Implementing AM Algorithm from Haario et al (2001) but with 
#####	 unbounded support.
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
e <- 0.1



###################    PROPOSAL
###     Proposal is a mixture of Gaussian 

propose <- function(n, x, S){
	# n is the number of current step
	# X the generated process
	
	# For n <= 2d, sample from a fixed kernel. For n > 2d, sample from a mixture kernel.

	if( n <= 2*d){
		Y <- rnorm(d, mean = x, sd = 2.38/sqrt(d))
	}
	else{
		Z <- rnorm(d, mean = x, sd = 2.38/sqrt(d))
		Y <- chol(S + diag(rep(e, d)))%*%Z
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

N <- 1e5    						# number of samples
beta <- 0.05
p <- 0  							# records acceptances

X <- matrix(0, nrow = N, ncol = d)
X[1, ] <- rep(0, d)
#X[1, ] <- rnorm(d, mean = 0, sd = diag(M%*%t(M)))
current <- X[1, ]
S <- diag(rep(1, d))


for(i in 2:N){
	
  if(i%%10000 == 0) print(i)
  
	prop <- propose(i, X[i-1, ], S)
	move <- alpha(X[i-1, ], prop)
	
	X[i, ] <- move[1:d]
	p <- p + move[d+1]
	
	S <- cov(X[1:i, ])
}

par(mfrow = c(1, 1))
pdf("plots_haario.pdf")
for(k in 1:d){
  ts.plot(X[, k], type = "l", col = "black", main = "Trace plot for AM")
}
dev.off()

p/N


