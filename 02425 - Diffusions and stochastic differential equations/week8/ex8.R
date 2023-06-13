library(SDEtools)
library(fields)
library(Matrix)
source("C:\\Users\\anton\\OneDrive - Københavns Universitet\\Uni\\Uni\\11. semester\\SDE\\week8\\fvade.R")

# Define dynamics
## Drift and intensity
f <- function(x) lambda * (xi - x)
g <- function(x) gamma*sqrt(abs(x))
## Diffusivity and its spatial derivative
D <- function(x) 0.5*gamma^2*x
Dp <- function(x) 0.5*gamma^2
## Advective flow field
u <- function(x) f(x) - Dp(x)
## Parameters
xi <- 2
gamma <- 1
lambda <- 1/2
## Spatial grid
xmax <- 10
xv <- seq(0,xmax,length=301)
dx <- diff(xv)
xc <- xv[-1] - 0.5*dx
## Discretize the generator
G <- fvade(u,D,xv,'r')
## Loading required package: Matrix
##
## Attaching package: 'Matrix'
## The following object is masked from 'package:spam':
##
## det
## Compute the stationary density
pi <- StationaryDistribution(G)
## user system elapsed
## 0.008 0.000 0.008
phi <- pi / dx
plot(xc,phi,type="l")
plot(function(x)dgamma(x,rate=2*lambda/gamma^2,shape=2*lambda*xi/gamma^2),
     from=0,to=xmax,add=TRUE,col="red")

# PDF

## Transients
## Initial condition for the SDE
x0 <- xi/2
## Initial condition for the FKE
phi0 <- numeric(length(xc))
phi0[sum(xc<x0)] <- 1
## Time grid
tv <- seq(0,10,0.03333333)
## Solve the FKE
Gs = as(G, "sparseMatrix")
PHI <- sapply(tv,function(t) as.numeric(phi0 %*% expm(Gs*t)))/dx

gg = rev(sort(PHI))
ind_max = which(PHI == max(PHI), arr.ind = TRUE)
PHI[ind_max] = gg[2]+10*(gg[2]-gg[3])

image.plot(tv,xc,t(PHI))
text(1,9,"text")


plot(xc,PHI[,1])
plot(xc,PHI[,2])
plot(xc,PHI[,3])
plot(xc,PHI[,4])
plot(xc,PHI[,5])
plot(xc,PHI[,6])
plot(xc,PHI[,7])
plot(xc,PHI[,8])
plot(xc,PHI[,9])
plot(xc,PHI[,10])
plot(xc,PHI[,20])
plot(xc,PHI[,30])
plot(xc,PHI[,40])
plot(xc,PHI[,50])
plot(xc,PHI[,60])
plot(xc,PHI[,70])
plot(xc,PHI[,80])
plot(xc,PHI[,90])
plot(xc,PHI[,length(tv)])
lines(xc,phi)

#CDF
CDF <- apply(PHI*dx,2,cumsum)
image.plot(tv,xc,t(CDF))
EX <- apply(PHI*dx*xc,2,sum)
lines(tv,EX,lwd=2)
EX2 <- apply(PHI*dx*xc^2,2,sum)
VX <- EX2 - EX^2
sX <- sqrt(VX)
lines(tv,EX+sX,lty="dashed")
lines(tv,EX-sX,lty="dashed")

