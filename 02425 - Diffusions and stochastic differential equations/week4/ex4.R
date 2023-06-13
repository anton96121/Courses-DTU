
#### EX2 ####
## Helper function for simulation of Brownian motion, as per week 3
rBM <- function(tvec){
  return(cumsum(rnorm(length(tvec),mean=0,sd=sqrt(diff(c(0,tvec))))))
}

## System parameters
m <- 1 # [kg]
k <- 0.5 # [N/m]
c <- 0.2 # [N*s/m]
sigma <- 10 # [N sqrt(s)]

## Q 1.1
A <- array(c(0,-k/m,1,-c/m),c(2,2))
G <- array(c(0,sigma/m),c(2,1))

## Simulation parameters
T <- 1000
dt <- 0.01

## Setup arrays
tvec <- seq(0,T,dt)
P <- array(0,c(2,2,length(tvec))) # This also computes the solution of the Lyapunov equation
X <- array(0,c(2,length(tvec)))

## Simulate sample path of Brownian motion
B <- rBM(tvec)
dB <- diff(B)

## Main time loop, Euler stepping the SDE and the Lyapunov equation
for(i in 1:(length(tvec)-1)){
  X[,i+1] <- X[,i] + A %*% X[,i] * dt + G * dB[i]
  P[,,i+1] <- P[,,i] + (A %*% P[,,i] + P[,,i] %*% t(A) + G %*% t(G)) * dt
}

## Q 2.2: Plot the sample path; first part only for clarity
par(mfrow=c(2,1))
plot(tvec,X[1,],type="l",xlim=c(0,100),xlab="Time",ylab="Position")
plot(tvec,X[2,],type="l",xlim=c(0,100),xlab="Time",ylab="Velocity")


#### EX3 ####
lyap <- function(A,Q)
{
  A <- as.matrix(A)
  I <- diag(rep(1,nrow(A)))
  P <- kronecker(I,A)+kronecker(A,I)
  X <- -solve(P,as.numeric(Q))
  return(matrix(X,nrow=nrow(A)))
}
Pinf <- P[,,length(tvec)]
Pinf2 <- lyap(A,G%*%t(G))


## Print the empirical covariance and compare with the solution of the Lyapunov equation
print(cov(t(X)))
print(Pinf)
print(Pinf2)

#### EX5 ####
require(expm)

ivec <- seq(0,50/dt,1)
tvec <- ivec * dt
rhovec <- sapply(tvec,function(t) (Pinf2 %*% expm(t(A)*t) ))

par(mfrow=c(1,1))
acf(X[1,],lag.max=50/dt)
lines(ivec,rhovec[1,]/rhovec[1,1],col="red",lwd=3)
acf(X[2,],lag.max=50/dt)
lines(ivec,rhovec[4,]/rhovec[4,1],col="red",lwd=3)










#### EX6 ####
I <- diag(2)
H <- function(w) (solve(1i*w*I-A) %*% G)
ws = seq(0.2,18,0.01)
Hs <- sapply(ws,H)

par(mfrow=c(3,2))
min_y = min(min(abs(Hs[1,])),min(abs(Hs[2,])))
max_y = max(max(abs(Hs[1,])),max(abs(Hs[2,])))
plot(ws,abs(Hs[1,]),log="xy",type="l",ylab="Abs(H)", main = "Position", ylim = c(min_y,max_y))
plot(ws,abs(Hs[2,]),log="xy",type="l",ylab="Abs(H)", main = "Velocity", ylim = c(min_y,max_y))

min_y = min(min(Arg(Hs[1,])),min(Arg(Hs[2,])))
max_y = max(max(Arg(Hs[1,])),max(Arg(Hs[2,])))
plot(ws,Arg(Hs[1,]),log="x",type="l",ylab="Arg(H)", ylim = c(-pi,pi))
plot(ws,Arg(Hs[2,]),log="x",type="l",ylab="Arg(H)", ylim = c(-pi,pi))

min_y = min(min(abs(Hs[1,])^2*sigma^2),min(abs(Hs[2,])^2*sigma^2))
max_y = max(max(abs(Hs[1,])^2*sigma^2),max(abs(Hs[2,])^2*sigma^2))
plot(ws,abs(Hs[1,])^2*sigma^2,type="l",log="xy",ylab=expression(S[Q](omega)), ylim = c(min_y,max_y))
plot(ws,abs(Hs[2,])^2*sigma^2,type="l",log="xy",ylab=expression(S[V](omega)), ylim = c(min_y,max_y))


