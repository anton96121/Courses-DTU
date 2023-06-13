library(matlib)
library(MASS)
saddelmux = c(-4,-3,-2,-1,0)
saddelmuxp = c(0,1,2,3,4)
saddelmuy = c(-8,-6,-4,-2,0)
saddelmuyp = c(0,0,0,0,0)

gg = c(-10,-9,-8,-7,-6,-5,-4)
ggy = c(0,0,0,0,0,0,0)

gg2 = c(10,9,8,7,6,5,4)
ggy = c(0,0,0,0,0,0,0)

centermux = c(-4,-3,-2,-1,0)
centermupx = c(0,1,2,3,4)
centermuy = c(0,0,0,0,0)
centermuyp = c(0,2,4,6,8)

center_col = 'steelblue'
saddle_col = 'coral'

eqscplot(saddelmux,saddelmuy,type = 'l', xlim = c(-4,4), ylim = c(-8,8), lty = 2, lwd = 2,col = saddle_col, main = 'Bifurcation Diagram',xlab = expression(mu), ylab = expression(x[1]), cex.lab = 1.4)
lines(saddelmuxp,saddelmuyp, lty = 2, col = saddle_col, lwd = 2)
lines(gg,ggy, lty = 1, col = center_col, lwd = 2)
lines(centermux,centermuy, lty = 1, col = center_col, lwd = 2)
lines(centermupx,centermuyp, lty = 1, col = center_col, lwd = 2)
lines(gg2,ggy, lty = 2, col = saddle_col, lwd = 2)
grid()
legend("topright", legend = c("Center (Stable)", "Saddle (Unstable)"),
       lwd = 2, col = c(center_col,saddle_col), lty = c('solid','dashed'), cex=1, bty="n")

