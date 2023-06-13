## Main simulation loop of Np particles, simulation time Tstop,
## time step dt,
## with floating velocity v and diffusivity D
EulerSim <- function(v=0,D=1,H=1,Tstop=10,dt=1e-3,Np=1,Z0=H/2,u=function(z)0)
{
  ## Simulation duration and time step
  Nt <- ceiling(Tstop/dt)+1
  ## Initialize horiz. and vert. position
  Z <- X <- array(NA,c(Nt,Np))
  X[1,] <- 0
  Z[1,] <- Z0
  ## Useful shorthands
  s <- sqrt(2*D*dt)
  vdt <- v*dt
  ## Main time loop
  for(i in 2:Nt)
  {
    ## Update vertical position; reflect at bottom
    Z[i,] <- abs(Z[i-1,] + vdt + rnorm(Np,sd=s))
    ## Reflect at surface
    Z[i,] <- H - abs(H-Z[i,])
    ## Update horizontal position
    X[i,] <- X[i-1,] + u(Z[i,])*dt
  }
  return(list(Z=Z,X=X))
}

#### Question 1a ####
Tlong <- 100
dt <- 1e-3
v <- 2
H <- 1
D <- 1

# Positive drift, 2
VeryLongSimUp <- EulerSim(v=v,D=D,Tstop=Tlong,dt=dt,Np=1)
# Negative drift, -2
VeryLongSimDown <- EulerSim(v=-v,D=D,Tstop=Tlong,dt=dt,Np=1)

tv <- seq(0,Tlong,dt)
plot(tv,VeryLongSimUp$Z,type="l",col="green",xlab="Time",ylab="z",
     xlim=c(0,2))
lines(tv,VeryLongSimDown$Z,type="l",col="red")

breaks <- seq(0,H,length=21)
hist(VeryLongSimUp$Z,breaks=breaks,freq=FALSE,col=rgb(0,1,0,0.25),
     main="",xlab="z")
hist(VeryLongSimDown$Z,breaks=breaks,freq=FALSE,col=rgb(1,0,0,0.25),add=TRUE)
# The densities are from question 2 which can be found in the folder week1
plot(function(z)exp(z*v/D)*v/D/(exp(v*H/D)-1),from=0,to=1,add=TRUE)
plot(function(z)-exp(z*-v/D)*v/D/(exp(-v*H/D)-1),from=0,to=1,add=TRUE)


#### Question 1b ####
# Simulate a short run; compare with Gaussian
Tshort <- 0.01
Np <- 1000
ShortSim <- EulerSim(v=v,D=D,H=H,Tstop=Tshort,dt=1e-4,Np=Np)
hist(tail(ShortSim$Z,1),freq=FALSE,breaks=21,xlab="z",main="")

# Density is from page 18 and given as N(x0+v*t,2*D*t).
plot(function(z)dnorm(z,mean=H/2+v*Tshort,sd=sqrt(2*D*Tshort)),
     from=0,to=1,add=TRUE)



#### Question 3 ####
# See pdf in week1 for explanation.
Cs <- function(Pe)function(zp)exp(Pe*zp)*Pe/(exp(Pe)-1)
plot(Cs(10),from=0,to=1,ylim=c(0,10),xlab="z/H",ylab="C")
plot(Cs(1),from=0,to=1,add=TRUE,col="red")
plot(Cs(0.1),from=0,to=1,add=TRUE,col="green")


#### Question 4 ####
Np <- 1000
## Horiztonal flow field
u <- function(z) log(1+3*z)
FloatSim <- EulerSim(v=v,D=D,Tstop=10,Np=Np,u=u)
SinkSim <- EulerSim(v=-v,D=D,Tstop=10,Np=Np,u=u)
hist(tail(FloatSim$X,1),breaks=seq(0,15,0.25),freq=FALSE,col="green",xlab="X",main="")
hist(tail(SinkSim$X,1),freq=FALSE,breaks=seq(0,15,0.25),add=TRUE,col="red")
legend(x="topleft",legend=c('Light particles','Heavy particles'),fill=c("green","red"))

meanFloat <- mean(tail(FloatSim$X,1))
varFloat <- var(as.numeric(tail(FloatSim$X,1)))
meanSink <- mean(tail(SinkSim$X,1))
varSink <- var(as.numeric(tail(SinkSim$X,1)))
tab <- array(c(meanFloat,meanSink,varFloat,varSink),c(2,2))
colnames(tab) <- c("Mean","Var")
rownames(tab) <- c("Float","Sink")
print(tab)

