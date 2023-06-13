require(SDEtools)
library(fields)
## Loading required package: SDEtools
lambda <- 1
xi <- 1
gamma <- 1
# Define model.
f = function(x) lambda*(xi-x);
g = function(x) gamma*sqrt(x);
# Advection-diffusion form
D = function(x) 0.5*gamma^2*x;
u = function(x) f(x) - 0.5*gamma^2;
## Generator
## Define grid
h = 0.1
xii = seq(0,5,h) # Cell Interfaces
dx <- diff(xii)
xc = 0.5*(tail(xii,-1) + head(xii,-1)) # Cell centers
nc <- length(xc)
G = fvade(u,D,xii,'r'); # Generator


P <- as.matrix(expm(G*h)) # Probability transition matrix in the discritized markov chain

contour(xc,xc,P/rep(dx,rep(nc,nc)),xlab="From",ylab="To")
abline(0,1,lty="dashed")

image.plot(xc,xc,P/rep(dx,rep(nc,nc)))

mu <- P %*% xc
dmudt <- (mu - xc) / tsample
V <- P %*% (xc^2) - mu^2
dVdt <- V/tsample
par(mfrow=c(2,1))
plot(xc,dmudt,type="l",lwd=3,xlab="x",ylab="Drift")
lines(xc,f(xc))
plot(xc,dVdt,type="l",lwd=3,xlab="x",ylab="Growth rate in variance")
lines(xc,g(xc)^2)

#ex2
par(mfrow=c(1,1))
obs <- read.table("hmm-obs.txt",header=TRUE)
plot(Y ~t, data=obs)

#ex3
nu <- 0.5
likelihood <- function(x,y) dpois(y,lambda=nu*x)
ys <- seq(0,5, by =1)
dltab = outer(xc,ys,likelihood);
matplot(xc,dltab,type="l")
legend("topright",lty=1+ys,legend=ys,col=1+ys)

#ex4
ltab <- dltab[,obs$Y+1]
ltab2 = outer(xc,obs$Y,likelihood);
image.plot(xc,1:length(obs$Y),ltab)
image.plot(xc,1:length(obs$Y),ltab2)

#ex5
#P =expm(G* dt)



# Ex6
sim <- read.table("hmm-states.txt",header=TRUE)
time = obs$t
estimated_distribution = est$psi
state = xc
image.plot(time,state,estimated_distribution, zlim = c(0,sort(est$psi, decreasing = TRUE)[2]),col=grey(seq(0, 1, length = 256)))
lines(sim$t,sim$X,lwd=1, col = 'green')
lines(obs$t,apply(est$psi,1,function(p)sum(p*xc)),lwd = 2, col = 'red')

# Ex7

hmmfilter <- function(G)
{
  phi <- Matrix(array(0,c(length(obs$t),length(xc))))
  psi <- phi
  ## Initialize with the stationary distribution
  mu <- StationaryDistribution(G)
  mu <- mu/sum(mu)
  ## Compute transition probabilities
  P <- expm(G*tsample)
  const <- numeric(length(obs$t))
  phi[1,] <- mu
  ## Include the first data update
  psi[1,] <- phi[1,] * ltab[,1]
  const[1] <- sum(psi[,1] )
  psi[,1] <- psi[,1] / const[1]
  ## Main time loop over the remaining time steps
  for(i in 2:length(obs$t))
  {
    phi[i,] = psi[i-1,] %*% P # Time update
    psi[i,] = phi[i,] * ltab[,i] # Data update
    const[i] <- sum(psi[i,]) # Normalization
    psi[i,] = psi[i,] / const[i]
  }
  return(list(c=const,phi=as,matrix(phi),psi=as.matrix(psi),loglik=sum(log(const))))
}

## Likelihood estimation of xi
xis <- seq(0,2,0.05) # Try these values of xi
loglik <- function(xi) # Likelihood function of one value of xi
{
  f <- function(x) lambda*(xi-x);
  D <- function(x) 0.5*gamma^2*x;
  u <- function(x) f(x) - 0.5*gamma^2;
  G = fvade(u,D,xii,'r');
  return(hmmfilter(G)$loglik)
}
## Tabulate likelihood function
ls <- sapply(xis,loglik)
plot(xis,ls,type="b")
grid()
abline(v=xi) ## Include the true value used in the simulation
abline(h=max(ls)-0.5*qchisq(0.95,df=1)) ## Include confidence interval



