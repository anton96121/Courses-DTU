setwd('C:\\Users\\anton\\OneDrive - Københavns Universitet\\Uni\\Uni\\10. semester\\02424 - adv stat\\Lec5')
data = read.csv("AIDS.csv", header = T, sep = ",")
plot(data$cases)


cases.glm.M0 <- glm(formula = cases ~ year*quarter, family = poisson(link = 'log'), data = data)
drop1(cases.glm.M0,test = 'Chisq')

cases.glm.M1 <- update(cases.glm.M0,.~.-year:quarter)
drop1(cases.glm.M1,test = 'Chisq')

# Done