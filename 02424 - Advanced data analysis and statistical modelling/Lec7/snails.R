
##################################################
# Snails
##################################################
library("rstudioapi")
path_file_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path_file_dir)

snails <- read.table("snails.txt",header=TRUE)

head(snails)
snails$p <- snails$death/snails$n
library(scatterplot3d)
col <- numeric(dim(snails)[1])
pch <- numeric(dim(snails)[1])

col[snails$species=="A"] <- 1
col[snails$species=="B"] <- 2
pch[snails$exposure==1] <- 1
pch[snails$exposure==2] <- 2
pch[snails$exposure==3] <- 3
pch[snails$exposure==4] <- 4


scatterplot3d(snails[ ,c("humidity","temp","p")],color=col,pch=pch)

plot(p~humidity,color=col,pch=pch,data=snails)
plot(p~humidity,col=col,pch=pch,data=snails)
plot(p~temp,col=col,pch=pch,data=snails)
plot(p~exposure,col=col,pch=pch,data=snails)

snails$resp <- cbind(snails$death,snails$n-snails$death)

fit0 <- glm(resp~factor(humidity)*factor(temp)*factor(exposure)*species,data=snails,family=binomial)
fit1 <- update(fit0,.~.-factor(humidity):factor(temp):factor(exposure):species)
fit2 <- update(fit1,.~.-factor(humidity):factor(temp):factor(exposure)-
                 factor(humidity):factor(temp):species-
                 factor(humidity):factor(exposure):species-
                 factor(temp):factor(exposure):species)
fit3 <- update(fit2,.~.-factor(humidity):factor(exposure)- factor(humidity):factor(temp) -
                 factor(humidity):species)

fit4 <- update(fit3,.~.-factor(temp):species-factor(exposure):species )
fit5 <- update(fit4,.~.-factor(temp):factor(exposure))

fit6 <- glm(formula = resp ~ humidity + factor(temp) + factor(exposure) +
              species, family = binomial, data = snails)


fit7 <- glm(formula = resp ~ humidity + temp + factor(exposure) +
              species, family = binomial, data = snails)




anova(fit7,test="Chisq")

summary(fit1)

summary(fit2)

summary(fit3)

summary(fit4)

summary(fit5)

summary(fit7)

plot(fit7)

anova(fit1,fit0,test="Chisq")
anova(fit4,fit0,test="Chisq")
anova(fit5,fit0,test="Chisq")
anova(fit6,fit0,test="Chisq")
anova(fit7,fit6,test="Chisq")

anova(fit2,test="Chisq")
drop1(fit2,test="Chisq")

drop1(fit3,test="Chisq")

drop1(fit4,test="Chisq")

drop1(fit5,test="Chisq")

drop1(fit6,test="Chisq")


summary(fit0)

glm(resp~factor(humidity)+factor(temp)+factor(exposure)+species, family=binomial,data=snails)


glm(resp~factor(humidity),family=binomial,data=snails)



head(snails)


fit8 <- glm(formula = resp ~ humidity + temp + factor(exposure) +
              species, family = binomial, data = snails[snails$exposure!=1, ])


fit9 <- glm(formula = resp ~ humidity + temp + I(exposure-1) + I((exposure-1)^2) +
              species, family = binomial, data = snails)

summary(fit8)
summary(fit9)

anova(fit9,fit7,test="Chisq")

summary(fit7)

4.10485*1-0.70029*1^2
4.10485*2-0.70029*2^2
4.10485*3-0.70029*3^2


fit10 <- glm(formula = resp ~ humidity+ temp  + species+ I(exposure-1) + I((exposure-1)^2),
             family = binomial, data = snails)

summary(fit10)
drop1(fit10,test="Chisq")

anova(fit10,fit9,test="Chisq")

hum <- seq(40,90,by=1)
temp <- seq(5,30,by=1)
pred1 <- predict(fit10, newdata = data.frame(humidity=hum[1], temp=temp,
                                             species = "B", exposure = 4))

s3d <- scatterplot3d(cbind(hum[1],temp,pred1),type="l",color=gray(0.5),
                     xlim=c(40,90),ylim=c(5,30),zlim=c(-3.5,4.5))

for(i in 2:length(hum)){
  pred1 <- predict(fit10, newdata = data.frame(humidity=hum[i],
                                               temp=temp, species = "B", exposure = 4))
  s3d$points3d(cbind(hum[i],temp,pred1),
               pch=19,type="l",col=gray(0.5))
  
}

for(i in 1:length(temp)){
  pred1 <- predict(fit10, newdata = data.frame(humidity=hum,
                                               temp=temp[i], species = "B", exposure = 4))
  s3d$points3d(cbind(hum,temp[i],pred1),
               pch=19,type="l",col=gray(0.5))
  
}



pred1 <- predict(fit10, newdata = data.frame(humidity=hum[1], temp=temp,
                                             species = "B", exposure = 4),
                 type="response")

s3d <- scatterplot3d(cbind(hum[1],temp,pred1),type="l",color=gray(0.5),
                     xlim=c(40,90),ylim=c(5,30),zlim=c(0,1))

I <- snails$species=="B" & snails$exposure==4
s3d$points3d(snails$humidity[I],snails$temp[I],snails$p[I],
             pch=19,type="p",col=2)


for(i in 2:length(hum)){
  pred1 <- predict(fit10, newdata = data.frame(humidity=hum[i],
                                               temp=temp, species = "B", exposure = 4),
                   type="response")
  s3d$points3d(cbind(hum[i],temp,pred1),
               pch=19,type="l",col=gray(0.5))
  
}

for(i in 1:length(temp)){
  pred1 <- predict(fit10, newdata = data.frame(humidity=hum,
                                               temp=temp[i], species = "B", exposure = 4),
                   type="response")
  s3d$points3d(cbind(hum,temp[i],pred1),
               pch=19,type="l",col=gray(0.5))
  
}



plot(fit10)


fit11 <- glm(formula = resp ~ exposure + I(exposure*exposure) +
               humidity+ temp  + species,
             family = binomial(link=logit), data = snails)

summary(fit11)

