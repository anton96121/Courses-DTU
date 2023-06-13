
#### Q3 ####
N <- 1000000
X <- runif(N)
Y <- runif(N)
Z <- X+Y

Z_means = rep(NA, 10)

for (i in 1:10){
  Z_means[i] = mean(Z[which(X<0.1*i & 0.1*(i-1)<=X)])
}

plot(seq(0,0.9,0.1)+0.05, Z_means, xlab = "X", ylab = "E[Z|X]", main = "Conditional mean", xlim = c(0,1), pch = 19, lwd = 1)
for (i in 1:11){
  abline(v=(i-1)/10)
}


#### Q4 ####
lines(seq(0,0.9,0.1)+0.05, seq(0,0.9,0.1)+0.05+1/2)

#### Q5 ####
X_means = rep(NA, 10)

for (i in 1:10){
  X_means[i] = mean(X[which(Z<0.2*i & 0.2*(i-1)<=Z)])
}
  
plot(seq(0,1.8,0.2)+0.1, X_means, xlab = "Z", ylab = "E[X|Z]", main = "Conditional mean", xlim = c(0,2), pch = 19, lwd = 1)
for (i in 1:11){
  abline(v=(2*i-2)/10)
}

lines(seq(0,1.8,0.2)+0.1, (seq(0,1.8,0.2)+0.1)/2)

#### Q6 ####

## Based on analytical expressions
mean(X)
mean(Z/2)

##
mean(X_means) 
mean(X)
