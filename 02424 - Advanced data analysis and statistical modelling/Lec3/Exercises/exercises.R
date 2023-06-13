setwd('C:\\Users\\anton\\OneDrive - Københavns Universitet\\Uni\\Uni\\10. semester\\02424 - adv stat\\Lec3\\Exercises')
############################### Q1 ##########################
# ex1
Y = c(4.4,3.4,3.3,2.5,7.3,4.9,4.8,4.4)
plot(Y, ylim = c(0,8), xlim = c(0,8))
grid()

# We assume a common slope, a intercept for the roof and a constant effect of being at the ground

X =matrix(c(1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,1,2,3,4,5,6,7,8),8,3)

#ex2
beta = solve(t(X)%*%X)%*%t(X)%*%log(Y)

#ex3
exp(t(c(1,1,9))%*%beta)


############################## Q2 ###########################
data = read.csv("SO2.csv", header = T, sep = "")

## ex1 ##
# Correlation panel
panel.cor <- function(x, y){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=2)
  txt <- paste0("corr = ", r)
  cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt)
}
# Customize upper panel
upper.panel<-function(x, y){
  points(x,y, pch = 19)
}
# Create the plots
pairs(data, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)

## ex2 ##
model1 = lm(data$Pollution ~ data$Industry + 1)
summary(model1)

model2 = lm(data$Pollution ~ data$Industry + data$Population + 1)
summary(model2)

anova(model1,model2)
# We see population is actually significant. Can maybe be a result due to the 
#  very large city from data point 22.

model3 = lm(data$Pollution ~ data$Industry + data$Population + data$Temp + 1)
summary(model3)

# we conclude at model 2.


## We now try to remove city 22
data2 = data[-22,]

pairs(data2, 
      lower.panel = panel.cor,
      upper.panel = upper.panel)

model1 = lm(data2$Pollution ~ data2$Industry + 1)
summary(model1)

model2 = lm(data2$Pollution ~ data2$Industry + data2$Population + 1)
summary(model2)

anova(model1,model2)

# Still significant

### ex3 ###
# see solution