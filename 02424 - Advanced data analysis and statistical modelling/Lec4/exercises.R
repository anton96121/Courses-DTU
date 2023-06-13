##### EX 3.4 #####

x <- rep(1:3,2)
method <- c(rep(0,3),rep(1,3))
y <- c(0.22,0.38,0.72,0.31,0.66,0.99)

plot(y ~ x, col=method+1, pch = method+1)

# covariance structure
v <- x^((method<=1)+(method==0))
Sigma <- diag(v)
SigmaInv <- solve(Sigma)

# design matrix with mixed effect
X0 <- cbind(1, method, x, method * x)

## q1
# estimated beta, eq 3.27, p. 49
beta <- solve(t(X0)%*%SigmaInv%*%X0)%*%t(X0)%*%SigmaInv%*%y

## q2
# estimated variance, thm 3.5, p. 53
sigma2 <- t(y-X0%*%beta)%*%SigmaInv%*%(y-X0%*%beta)/(length(x)-length(beta))

## q3
M0 = lm(y~1+factor(method)*x,weights = 1/v)
M1 = lm(y~1+factor(method):x,weights = 1/v)
M2 = lm(y~1+x,weights = 1/v)
M3 = lm(y~1,weights = 1/v)

anova(M1,M0)
anova(M2,M0)
anova(M3,M0) # we see removing the slope is too agressive. Then common slope(M2)

##### EX 3.7 #####
data(anscombe)
X <- anscombe
par(mfrow=c(2,2))
plot(X$y1 ~ X$x1, col=1, pch = 1, xlim = c(0,20), ylim = c(0,14))
plot(X$y2 ~ X$x2, col=2, pch = 2, xlim = c(0,20), ylim = c(0,14))
plot(X$y3 ~ X$x3, col=3, pch = 3, xlim = c(0,20), ylim = c(0,14))
plot(X$y4 ~ X$x4, col=4, pch = 4, xlim = c(0,20), ylim = c(0,14))

M1 <- lm(X$y1 ~ X$x1)
summary(M1)

M2 <- lm(X$y2 ~ X$x2)
summary(M2)

M3 <- lm(X$y3 ~ X$x3)
summary(M3)

M4 <- lm(X$y4 ~ X$x4)
summary(M4)


