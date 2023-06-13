library(MASS)
library(matlib)
library(extrafont)
loadfonts()
par(family = "Latin Modern Roman")

############### Opg 2.e #################3

expRt = function(t,a,x0){
  x = matrix(c(exp(a*t)*cos(sqrt(-a^2 + 1)*t),exp(a*t)*sin(sqrt(-a^2 + 1)*t), -exp(a*t)*sin(sqrt(-a^2 + 1)*t), exp(a*t)*cos(sqrt(-a^2 + 1)*t)),2,2)
  y = x%*%x0
  return(y)
}

y = 2

aus = 0.4
tus = seq(0,pi/sqrt(1-aus^2),length.out = 1000)

acen = 0
tcen = seq(0,pi/sqrt(1-acen^2),length.out = 1000)

as = -0.4
ts = seq(0,pi/sqrt(1-as^2),length.out = 1000)

DataUnstable = matrix(0,1000,2)
DataCenter = matrix(0,1000,2)
DataStable = matrix(0,1000,2)
DataEnd = matrix(0,999,2)
ys = matrix(0,1000,2)
xs = matrix(0,1000,2)

yus = matrix(0,1000,2)
xus = matrix(0,1000,2)

ycen = matrix(0,1000,2)
xcen = matrix(0,1000,2)

axis = seq(-20,20,length.out = 100)

for(i in 1:1000){
  ys[i,] = c(-(as*axis[i])/(sqrt(-as^2+1)),axis[i])
  xs[i,] = c(axis[i]/sqrt(-as^2 + 1), 0)
  
  yus[i,] = c(-(aus*axis[i])/(sqrt(-aus^2+1)),axis[i])
  xus[i,] = c(axis[i]/sqrt(-aus^2 + 1), 0)
  
  ycen[i,] = c(-(acen*axis[i])/(sqrt(-acen^2+1)),axis[i])
  xcen[i,] = c(axis[i]/sqrt(-acen^2 + 1), 0)
  
}



for(i in 1:1000){
  DataUnstable[i,] = expRt(tus[i],aus, c(-aus*y/sqrt(-aus^2 + 1),y))
  DataCenter[i,] = expRt(tcen[i],acen, c(-acen*y/sqrt(-acen^2 + 1),y))
  DataStable[i,] = expRt(ts[i],as, c(-as*y/sqrt(-as^2 + 1),y))
  
  if(i!=1000){
    aend = i/500-1
    DataEnd[i,] = expRt(pi/sqrt(1-aend^2),aend, c(-aend*y/sqrt(-aend^2 + 1),y))
  }
  
}






eqscplot(DataUnstable[1,1],DataUnstable[1,2], col = 'red',cex = 1.25, pch = 19, xlim = c(-5,5), ylim = c(-9,3),ylab = "y(t)", xlab = 'x(t)')
lines(DataUnstable, col = 'red')
lines(yus, col = 'red', lty = 'dashed')
lines(DataStable, col = 'blue')
lines(ys, col = 'blue', lty = 'dashed')
lines(DataCenter, col = 'green')
lines(ycen, col = 'green', lty = 'dashed')
lines(xcen)
points(DataCenter[1,1],DataCenter[1,2], col="green", pch=19)
points(DataStable[1,1],DataStable[1,2], col="blue", pch=19)
points(DataCenter[1000,1],DataCenter[1000,2], col="green", pch=19)
points(DataStable[1000,1],DataStable[1000,2], col="blue", pch=19)
points(DataUnstable[1000,1],DataUnstable[1000,2], col="red", pch=19)
lines(DataEnd, lty = 'dashed')
lines(c(-15, 15), c(2,2), lty = 'dotted')
legend("left", legend = c("a = -0.4", "a = 0", "a = 0.4","v-axis for a=-0.4","v-axis for a=0","v-axis for a=0.4", 't = 0', expression("t" == frac(pi,sqrt(1-a^2)))),
       lwd = 1, col = c('blue','green', 'red','blue','green', 'red','black','black'), lty = c('solid','solid','solid','dashed','dashed','dashed', 'dotted','dashed'), cex=0.8, bty="n")




## 

expRtAB = function(t,a,b,x0){
  expA = matrix(c(exp(a*t)*cos(sqrt(-a^2 + 1)*t),exp(a*t)*sin(sqrt(-a^2 + 1)*t), -exp(a*t)*sin(sqrt(-a^2 + 1)*t), exp(a*t)*cos(sqrt(-a^2 + 1)*t)),2,2)
  expB = matrix(c(exp(b*t)*cos(sqrt(-b^2 + 1)*t),exp(b*t)*sin(sqrt(-b^2 + 1)*t), -exp(b*t)*sin(sqrt(-b^2 + 1)*t), exp(b*t)*cos(sqrt(-b^2 + 1)*t)),2,2)
  Pa = matrix(c(sqrt(1-a^2),0,a,1),2,2)
  Pb = matrix(c(sqrt(1-b^2),0,b,1),2,2)
  x = Pa%*%expA%*%inv(Pa)%*%(Pb%*%expB%*%inv(Pb)%*%x0)
  return(x)
}


y = 2

aus = 0.55
bus = -0.45
tus = seq(0,pi/sqrt(1-aus^2),length.out = 1000)

acen = 0
bcen = 0
tcen = seq(0,pi/sqrt(1-acen^2),length.out = 1000)

as = 0
bs = 0.3
ts = seq(0,pi/sqrt(1-as^2),length.out = 1000)

DataUnstable = matrix(0,1000,2)
DataCenter = matrix(0,1000,2)
DataStable = matrix(0,1000,2)


for(i in 1:1000){
  DataUnstable[i,] = expRtAB(tus[i],aus,bus, c(0,y))
  DataCenter[i,] = expRtAB(tcen[i],acen,bcen, c(0,y))
  DataStable[i,] = expRtAB(ts[i],as,bs, c(0,y))
  
}

eqscplot(DataUnstable,type = 'l', col = 'red',cex = 1.25, pch = 19, xlim = c(-5,5), ylim = c(-9,22),ylab = "y(t)", xlab = 'x(t)')
lines(DataCenter, col = 'blue')
lines(DataStable, col = 'green')


