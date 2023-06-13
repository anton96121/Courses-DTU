library("rstudioapi") 
path_file_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path_file_dir)
dat<-read.table("challenger_data.txt", head=TRUE)
dat
par(mfrow=c(1,2))
plot(dat$temp, dat$failed, xlab='Temperature', ylab='No damaged (out of 6)')
plot(dat$pres, dat$failed, xlab='Pressure', ylab='No damaged (out of 6)')

dat$resp<-cbind(dat$failed,dat$n-dat$failed)
fit0<-glm(formula = resp ~ temp+pres,
          family = binomial(link = logit), 
          data = dat)
drop1(fit0, test='Chisq')


fit1<-glm(formula = resp ~ temp,
          family = binomial(link = logit), 
          data = dat)
drop1(fit1, test='Chisq')
## Done with reduction
summary(fit1)


temp <- seq(31,85)

pred <- predict(fit1,newdata = data.frame(temp=temp),type="response")

plot(temp,pred,type="l")

## Odds
exp(-10*coef(fit1)[2])

## Done.. 
