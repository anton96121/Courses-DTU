# ex2
x = seq(-1,8,length.out = 10)
y = c(1.4,4.7,5.1,8.3,9.0,14.5,14.0,13.4,19.2,18)
fit = lm(y~x)
plot(fit)


# ex3
x = c(13,5,28,28,15,4,13,4,10,17,11,13,12,17,3)


fun <- function(theta, x){
  - sum(dnbinom(x, theta[1], theta[2], log = TRUE))
}

# Warning is given because NaNs are produced due to boundary is visited
nlminb(c(1,0.5),fun,x=x,lower=c(0.0000001,0.0000001), upper = c(Inf,0.9999999))

# No problems because boundary is excluded
nlminb(c(1,0.5),fun,x=x,lower=c(0.0000001,0.0000001), upper = c(Inf,0.9999999))
