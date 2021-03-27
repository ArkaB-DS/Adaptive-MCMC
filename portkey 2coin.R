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

## MC function

MC <- function(N=1e5,b=0.75){
X <- numeric(N) # creating  the MC
X[1] <- .1
p=0

for (i in 2:N)
{
# print(i)
 Y <- rnorm(1, mean = X[i-1], sd = 0.13)
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
 return (list("chain"=X, "acc.prob." =p/N))
}


par(mfrow=c(2,2))

print(MC(b=0.99)$acc.prob.)
ts.plot(MC(b=0.99)$chain[.99e5:1e5],
xlab="Iteration",col="seagreen",ylab="",
main=expression(paste(beta," = 0.99")))

print(MC(b=0.90)$acc.prob.)
ts.plot(MC(b=0.90)$chain[.99e5:1e5],
xlab="Iteration",col="seagreen",ylab="",
main=expression(paste(beta," = 0.90")))

print(MC(b=0.75)$acc.prob.)
ts.plot(MC(b=0.75)$chain[.99e5:1e5],
xlab="Iteration",col="seagreen",ylab="",
main=expression(paste(beta," = 0.75")))

