
#### EX1 ####
brown_sim = function(t) cumsum(rnorm(length(t),mean=0,sd=sqrt(diff(c(0,t)))))
t = c(0,0.5,1.5,2)
N = 10000
B = sapply(1:N,function(i)brown_sim(t))

# We see that the mean vector is a vector of zeros as 
# the theory says on page 70 in the book.
print(apply(B,1,mean))
# Similarly we also see that the covariance matrix corresponds
# to the theory on page 70 in the book.
print(var(t(B)))


#### EX2 ####
# 2.1
t = seq(0,1,0.001)
N = 1000
B = sapply(1:N,function(i)brown_sim(t))

# 2.2
S1 = apply(B,2,max)
S1_pdf = function(x) 2*dnorm(-x/sqrt(1)) # thm 4.3.2
S1_cdf = function(x) 1-2*pnorm(-x/sqrt(1)) # thm 4.3.2

par(mfrow=c(1,2))
hist(S1, main = "PDF, S1 = max{Bt : 0 <= t <= 1}", freq = F)
curve(S1_pdf,from=0,to=max(S1), add = TRUE,col="red", lwd = 2)

plot(ecdf(S1), main = "CDF, S1 = max{Bt : 0 <= t <= 1}", freq = F)
curve(S1_cdf,from=0,to=max(S1), add = TRUE,col="red", lwd = 2)

# 2.3
b=0.5
tau = apply(B,2,function(x)min(which(b<=x),N))/N
tau_pdf = function(t) b*t^(-3/2)*dnorm(b/sqrt(t)) # thm 4.3.3
tau_cdf = function(t) 2*pnorm(-abs(b)/sqrt(t)) # thm 4.3.3

par(mfrow=c(1,2))
hist(tau, main = "PDF, tau = min{t : Bt >= b}", freq = F)
curve(tau_pdf,from=0,to=max(tau),add=TRUE,col="red", lwd = 2)

plot(ecdf(tau), main = "CDF, tau = min{t : Bt >= b}")
curve(tau_cdf,from=0,to=max(tau),add=TRUE,col="red", lwd = 2)

#### EX3 ####
set.seed(1)
# 3.1
delta = 2^(-24)
t = seq(0,1,delta)
N = 1
B = sapply(1:N,function(i)brown_sim(t))

# 3.2
dB = abs(diff(B))

(V_d = sum(dB)) # Def 4.3.1
(V2_d = sum(dB^2)) # Def 4.3.2

# 3.3 
B2 = B[seq(1,length(t),2)]

# 3.4
halfs = 20
v_d_vec = rep(0,halfs)
v2_d_vec = rep(0,halfs)

B_half = B

dB = abs(diff(B_half))
v_d_vec[1] = sum(dB)
v2_d_vec[1] = sum(dB^2)
for (i in 2:halfs){
  B_half = B_half[seq(1,length(B_half),2)]
  
  dB = abs(diff(B_half))
  v_d_vec[i] = sum(dB)
  v2_d_vec[i] = sum(dB^2)
}


deltas <- delta*2^(1:halfs)/2

par(mfrow=c(1,2))
plot(log10(deltas),v_d_vec,xlab='log10(delta)',ylab="Discretized total variation V(B)",pch=16,log="y")
lines(log10(deltas),sqrt(2/pi/deltas))
grid()
plot(log10(deltas),v2_d_vec,xlab='log10(delta)',ylim = c(0.7,1.3),ylab=expression("Discretized quadratic variation [B]"[1]), log = "y")
lines(log10(deltas),v2_d_vec)
grid()
#abline(h=1)

# We see that the total variation goes towards 1 (because t = 1) when delta goes towards 0.
# This which agrees with thm 4.3.1.
  
