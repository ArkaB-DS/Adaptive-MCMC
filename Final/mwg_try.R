
K <- 10        # no of theta parameters
r <- seq(5, 95, by = 10)
Y <- numeric(0)

for(i in 1:K){
	Y <- c(Y, rnorm(r[i], mean = i - 1, sd = 10))
}


B <- 2e3 

###########  Declaration and Initialization of parameter vectors

A <- rep(0, 50*B + 1)
A[1] <- 0.5

u <- rep(0, 50*B + 1)

theta <- matrix(0, nrow = 50*B + 1, ncol = K)

V <- rep(0, 50*B + 1)
V[1] <- (sum((Y - rep(theta[1, ], r))^2)/2 + 1)/(sum(r)/2 + 1)


#################################################################


ls <- rep(0, K + 3)       #####  ls ~ [A, u, theta1, theta2, ... , thetaK, V]
p <- rep(0, K + 3)
sig <- exp(ls)

for(n in 1:B){

	for(i in 2:51){

		j = (n - 1)*50 + i
		y <- rnorm(K+3, mean = c(A[j-1], u[j-1], theta[j-1, ], V[j-1]), sd = sig)
		z <- log(runif(K+3))


		### Updating A

		r <- 1/A[j-1] - 1/y[1] + 2*(log(A[j-1]) - log(y[1]))

		if(z[1] < 1){
			A[j] <- y[1]
			p[1] <- p[1] + 1
		}else{
			A[j] <- A[j-1]
		}


		### Updating u

		r <- (u[j-1]^2 - y[2]^2)/2

		if(z[2] < 1){
			u[j] <- y[2]
			p[2] <- p[2] + 1
		}else{
			u[j] <- u[j-1]
		}


		### Updating theta

		for(k in 1:K){
 
			r <- log(A[j]^2 + (theta[j-1, k] - u[j])^2) - log(A[j]^2 + (y[2+k] - u[j])^2)

			if(z[2+k] < 1){
				theta[j, k] <- y[2+k]
				p[2+k] <- p[2+k] + 1
			}else{
				theta[j, k] <- theta[j-1, k]
			}
		}

		### Updating V

		r <- (sum((Y - rep(theta[j, ], r))^2)/2 + 1)*(1/V[j-1] - 1/y[K+3]) + (sum(r)/2 + 2)*(log(V[j-1]) - log(y[K+3]))

		if(z[K+3] < 1){
			V[j] <- y[K+3]
			p[K+3] <- p[K+3] + 1
		}else{
			V[j] <- V[j-1]
		}

	}

	####  Updating ls and sigma

	t <- as.numeric(p/(n*50) <= 0.44)
	ls <- ls + ((-1)^(t))*min(0.01, n^(-1/2))
	sig <- exp(ls)
}

print(p)
print(ls)