### optimal scaling nd motivation for adaptive mcmc
set.seed(1011)

# target - N(0,1)
# proposal - N(x,h) where h is to be chosen

N <- 1e5  # length of the chain

rwmh <- function(N=1e3,h){
X <- numeric(N) # creating MC
X <- -1 # starting value
p <- 0 # for calculating acc. prob.
# implementing RWMH
for(i in 2:N)
{
 Y <- rnorm(1,mean=X[i-1],sd=sqrt(h))
 if(runif(1)<dnorm(Y)/dnorm(X[i-1])) 
 {
  X[i] <- Y
  p <- p+1
  } else {
  X[i] <- X[i-1]
  }
} 
return( list(acc.prob=p/N, chain=X) )
}

# plots

#par(mfrow=c(2,3))
layout(matrix(c(1,2,3,4,5,6,7,7,7), ncol=3, byrow=TRUE), heights=c(4,4,1))
plot.ts(rwmh(h=1000)$chain,col="seagreen",
xlab="(a) Proposal variance too large",ylab="d1")
ts.plot(rwmh(h=0.0001)$chain,col="seagreen",
xlab="(b) Proposal variance too small",ylab="d2")
ts.plot(rwmh(h=2)$chain,col="seagreen",
xlab="(c) Proposal variance approximately optimized",ylab="d3")
acf(rwmh(h=1000)$chain,col="seagreen",main="")
acf(rwmh(h=0.0001)$chain,col="seagreen",main="")
acf(rwmh(h=2)$chain,col="seagreen",main="")

title("Fig. 2. Simple Metropolis algorithm with (a) too-large variance (left plots), (b) too-small variance (middle) and (c) appropriate variance
(right). Trace plots (top) and autocorrelation plots (below) are shown for each case.",
line = -45, outer = TRUE,font.main=3)
