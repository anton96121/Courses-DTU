library(SDEtools)
library(fields)
library(Matrix)
library(Bessel)
source("C:\\Users\\anton\\OneDrive - Københavns Universitet\\Uni\\Uni\\11. semester\\SDE\\week9\\fvade.R")

#brow
rBM <- function(tvec){
  return(cumsum(rnorm(length(tvec),mean=0,sd=sqrt(diff(c(0,tvec))))))
}

# Von Misses distribution
VonMises  = function(x,kappa){
  I0 = BesselI(kappa, 0, expon.scaled = F)
  pdf_probs = exp(kappa*cos(x))/(2*3.141593*I0)
  return(pdf_probs)
  
}
# Define dynamics
sigma = 1
## Drift and intensity
f <- function(x) -sin(x)
g <- function(x) sigma
## Diffusivity and its spatial derivative
D <- function(x) 1/2*sigma^2
Dp <- function(x) 0
## Advective flow field
u <- function(x) f(x) - Dp(x)
## Parameters
xmin <- -3.141593
## Spatial grid
xmax <- 3.141593
xv <- seq(xmin,xmax,length=301)
dx <- diff(xv)
xc <- xv[-1] - 0.5*dx
G <- fvade(u,D,xv,'r')


Tend = 100
times = seq(0,Tend,length=100*Tend+1)
x0 = 0
path <- heun(f,g,times,x0)

par(mfrow=c(1,1))
plot(path$times,path$X)


pi <- StationaryDistribution(G)
phi <- pi / dx

hist(path$X, density = T, freq = F, main = "Von Mises SDE")
lines(xc,phi,type="l")
lines(xc,VonMises(xc, kappa = 2/sigma^2), col = 'red')
legend("topleft",
       c("Simuation", "Numerical PDF", "Analytical PDF"), fill=c("gray","black", "red"))


# Ex 4
acf(cos(path$X), lag.max = 500)

# Ex 5
Gs = as(G, "sparseMatrix")
expLt = expm(Gs*1)
# mu cal in maple int(1/(2*Pi*2.279585)*exp(2/1*cos(x))*cos(x), x= -Pi..Pi)
mu = 0.6977747505
h = function(x) return(cos(x))
rho = function(x) return(VonMises(x,2))
rho_sum = rho(xc)/sum(rho(xc))
hvec = h(xc)


sum((rho_sum%*%(hvec-mu))[1]*(expLt%*%hvec))

sum((rho_sum*(hvec-mu))*(expLt%*%hvec))



