set.seed(1)

### portkey 2-coin algorithm
p2c <- function(beta=0.75,x,y,k=10)
{
 while(1){
 S <- rbinom(n = 1, size = 1, prob = beta)
 if (S==0) return(0)
 if (S==1)
 {
  C1 <- rbinom(n = 1, size = 1, prob = x/(x+y) )
  if (C1==1)
  {
   C2 <- ifelse(runif(1)<(dweibull(y, scale=rgamma(1,10,100),shape=k)/(k/exp(1)*y)),1,0)
   if (C2==1) {
    return(1) } else (next)
  }
 }
  else
  {
   C2 <- ifelse(runif(1)<(dweibull(x, scale=rgamma(1,10,100),shape=k)/(k/exp(1)*x)),1,0)
   if (C2==1) {
	return(0) } else (next)
  }
         }
}

## MC function for adaptive portkey
# X for adaptive

# beta---0.05 (for adaptive MCMC)

MCA <- function(N=1e5,b=0.75){
X <- numeric(N) # creating  the MC
X[1] <- .1

p=0 # adaptive 
for (i in 2:N)
{ print(i)
 Sigma_n <- (N-1)/N*var(X)
 Y <- ifelse(rbinom(n=1,size=1,prob=0.05),rnorm(n=1,mean=X[i-1],sd=0.1),
             rnorm(n=1,mean=X[i-1],sd=sqrt(2.38^2*Sigma_n)))
 if (Y<0) {
 X[i] <- X[i-1]
 next }
 if(p2c(x=X[i-1],y=Y,beta=b))
 {
  X[i] <- Y
  p <- p+1
 } else{
 X[i] <- X[i-1]
       } 
}
 return (list("chain"=X,"acc.prob" =p/N))
}

## MC function for non-adaptive portkey

MCNA <- function(N=1e5,b=0.75){
Z <- numeric(N) # creating  the MC
Z[1] <- .1

p=0 # adaptive 
for (i in 2:N)
{ print(i)
 Y <- rnorm(1,mean=Z[i-1],sd=0.13)
 if (Y<0) {
 Z[i] <- Z[i-1]
 next }
 if(p2c(x=Z[i-1],y=Y,beta=b))
 {
  Z[i] <- Y
  p <- p+1
 } else{
 Z[i] <- Z[i-1]
       } 
}
 return (list("chain"=Z,"acc.prob"=p/N))
}

A <- MCA(b=0.75)
B <- MCNA(b=0.75)


#par(mfrow=c(2,1))
print(A$acc.prob)
ts.plot(A$chain[.99e5:1e5],lwd=2,
xlab="Iteration",col="red",ylab="",
main=expression(paste(beta," = 0.75")))
print(B$acc.prob)
lines(B$chain[.99e5:1e5],
xlab="Iteration",col="blue",ylab="",
main=expression(paste(beta," = 0.75")))

par(mfrow=c(2,1))
acf(A$chain,
xlab="Iteration",col="red",ylab="",
main=expression(paste(beta," = 0.75")))
acf(B$chain,
xlab="Iteration",col="seagreen",ylab="",
main=expression(paste(beta," = 0.75")))






