library(tictoc)

K <- 10       # no of theta parameters
R <- round(seq(5, 500, length.out = K))
sumR <- c(0, cumsum(R))

Y <- numeric(0)

for(i in 1:K) Y <- c(Y, rnorm(R[i], mean = i - 1, sd = 10))

B <- 1e5

###########  Declaration and Initialization of parameter vectors

A <- rep(0, 50*B + 1)
#A[1] <- 0.5
A[1] <- 1

u <- rep(0, 50*B + 1)

theta <- matrix(0, nrow = 50*B + 1, ncol = K)

V <- rep(0, 50*B + 1)
#V[1] <- (sum((Y - rep(theta[1, ], R))^2)/2 + 1)/(sum(R)/2 + 1)
V[1] <- 1

#################################################################


ls <- matrix(0, nrow = B, ncol = K + 2)       #####  ls ~ [A, u, theta1, theta2, ... , thetaK, V]
p <- rep(0, K + 2)
sig <- exp(ls[1, ])

tic()

for(n in 1:B){

	if(n%%1000 == 0) print(n)

	for(i in 2:51){

		j = (n - 1)*50 + i
		y <- rnorm(K+2, mean = c(A[j-1], u[j-1], theta[j-1, ]), sd = sig)
		z <- log(runif(K+2))
		
		### Updating V

		temp <- rgamma(1, shape = sum(R)/2 + 1, rate = sum((Y - rep(theta[j-1, ], R))^2)/2 + 1)
		V[j] <- 1/temp

		### Updating A

		if(y[1] <= 0){
			A[j] <- A[j-1]
		}
		else{
			r <- (2*K - 2)*(log(y[1]) - log(A[j-1])) + (1/A[j-1] - 1/y[1]) + 
					sum(log(A[j-1]^2 + (theta[j-1, ] - u[j-1])^2)) -
					 sum(log(y[1]^2 + (theta[j-1, ] - u[j-1])^2))
			
			if(z[1] < r){
				A[j] <- y[1]
				p[1] <- p[1] + 1
			}else{
				A[j] <- A[j-1]
			}
		}
		
		### Updating V

##		temp <- rgamma(1, shape = sum(R)/2 + 1, rate = sum((Y - rep(theta[j, ], R))^2)/2 + 1)
##		V[j] <- 1/temp


		### Updating u

		r <- (u[j-1]^2 - y[2]^2)/2 + sum(log(A[j]^2 +
			 (theta[j-1, ] - u[j-1])^2)) - sum(log(A[j]^2 + (theta[j-1, ] - y[2])^2))

		if(z[2] < r){
			u[j] <- y[2]
			p[2] <- p[2] + 1
		}else{
			u[j] <- u[j-1]
		}


		### Updating theta

		for(k in 1:K){
 
			r <- (sum((Y[(sumR[k]+1):sumR[k+1]] - theta[j-1, k])^2) - sum((Y[(sumR[k]+1):sumR[k+1]] - y[2+k])^2))/(2*V[j])
					+ log(A[j]^2 + (theta[j-1, k] - u[j])^2) - log(A[j]^2 + (y[2+k] - u[j])^2)

			if(z[2+k] < r){
				theta[j, k] <- y[2+k]
				p[2+k] <- p[2+k] + 1
			}else{
				theta[j, k] <- theta[j-1, k]
			}
		}

	}

	####  Updating ls and sigma

	t <- as.numeric(p/(n*50) <= 0.44)
	ls[n + 1, ] <- ls[n, ] + ((-1)^(t))*min(0.01, n^(-1/2))
	sig <- exp(ls[n + 1, ])
}

toc()

print(p/(50*B))

pdf("MWG_try.pdf")
for(k in 1:K)  ts.plot(ls[, k], type = "l")
dev.off()

pdf("MWG_params.pdf")
ts.plot(A[1:1e5])
ts.plot(V[1:1e5])
ts.plot(u[1:1e5])
for(k in 1:K)  ts.plot(theta[1:1e5, k], type = "l")
dev.off()



