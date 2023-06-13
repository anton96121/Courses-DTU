library("rstudioapi")
library('lubridate')
library("xtable")
library("car")
library("lme4")
library(latex2exp)
library(zoo)
library(nlme)

confint_sigma2 <- function(object, level=0.95) {
  if (class(object) != "lm") stop("The object is not a lm object")
  alpha <- 1 - level
  n <- length(object$residuals)
  p <- length(coef(object))
  s2 <- sigma(object) ^ 2
  den_lower <- qchisq(p=alpha/2, df=n-p, lower.tail=FALSE)
  den_upper <- qchisq(p=1-alpha/2, df=n-p, lower.tail=FALSE)
  lower <- (n-p) * s2 / den_lower
  upper <- (n-p) * s2 / den_upper
  int <- matrix(c(lower, upper), ncol=2)
  colnames(int) <- paste0(100*c(alpha/2, 1-alpha/2), " %")
  rownames(int) <- "Sigma2"
  return(int)
}

path_file_dir = dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path_file_dir)


data = read.csv("../data/data_project.csv")
data$date = as.Date(data$date)
data <- data[data$date <="2008-12-31",]
data$weekday <- weekdays(data$date) 
data$weekday <- substr(data$weekday,1,3)
data$days_from_start <- data$date - data$date[1]
data$buss_day = 'buss_day'
data[data$weekday %in% c('Saturday','Sunday'),'buss_day'] = 'weekend'


#m1_future = read.csv("../data/data_future_found/data_m1_2002_to_2005.csv");
#m2_future = read.csv("../data/data_future_found/data_m2_2002_to_2005.csv");
#m3_future = read.csv("../data/data_future_found/data_m3_2002_to_2005.csv");
#m4_future = read.csv("../data/data_future_found/data_m4_2002_to_2005.csv");

all_futures_data = read.csv("../data/futures_data.csv")
all_futures_data[ , 1:3] <- lapply(all_futures_data[ , 2:3], as.Date)
all_futures_data = as.data.frame(all_futures_data)

date_futures_start <- min(all_futures_data[
  all_futures_data$forward_month=='M4',"delivery_start"])
# we want to years
date_spot_stop <- date_futures_start %m+% years(2)

#plot(all_futures_data[all_futures_data$forward_month=='M4',"price"])


# make shorter series for analysis
data_short <- data[data$date <date_spot_stop,]
par(mfrow=c(1,1))
plot(data_short$price)

#weekday_vec = c("Monday","Tuesday", "Wednesday", "Thursday", "Friday", "Saturday", "Sunday")
weekday_vec = c("Mon","Tue", "Wed", "Thu", "Fri", "Sat", "Sun")


#### Analysis for the large data set #####
### The very simple model 
par(mfrow=c(1,1))
m0 <- lm(log(price) ~ days_from_start+factor(weekday)-1, data = data)
summary(m0)    
par(mfrow=c(2,2))
plot(m0)

### This could be the simple model. Very noticeably, the residuals exhibit
# lots of variance that seems very systematic and need more variables to 
# explain
par(mfrow=c(1,1))
plot(resid(m0))



#### Analysis for the short data set ####
### The very simple model
png(file="../figures_report/application/short_period/simple_model_fit/short_lm_no_log_diag.png"
    , width=800
    , height=600)
par(mfrow=c(1,2))
m0_simple <- lm(price ~ days_from_start+factor(weekday)-1
                , data = data_short)
summary(m0_simple)    

table_out <- data.frame(confint(m0_simple)[,1])
table_out[,2] <- coef(m0_simple)
table_out[,3] <- confint(m0_simple)[,2]
xtable(table_out,digits=4)
sigma(m0_simple)^2
confint_sigma2(m0_simple)



par(mfrow=c(2,2))
plot(m0_simple)
dev.off()
### This could be the simple model. Very noticeably, the residuals exhibit
# lots of variance that seems very systematic
par(mfrow=c(1,1))

png(file="../figures_report/application/short_period/simple_model_fit/short_lm_no_log_fitted.png"
    , width=800
    , height=600)
plot(fitted(m0_simple), rstandard(m0_simple)
     , xlab='Fitted values'
     , ylab='Standardized residuals')
dev.off()

png(file="../figures_report/application/short_period/simple_model_fit/short_lm_no_log_residuals_time_series.png"
    , width=800
    , height=600)
plot(data_short$date,rstandard(m0_simple)
     , ylab='standardized residuals'
     , xlab='date')
dev.off()

### introduce log transformation
m0_simple_log <- lm(log(price) ~ days_from_start+factor(weekday)-1
                    , data = data_short)
summary(m0_simple_log)    

png(file="../figures_report/application/short_period/simple_model_fit/short_lm_w_log_diag.png"
    , width=800
    , height=600)
par(mfrow=c(2,2))
plot(m0_simple_log)
dev.off()


png(file="../figures_report/application/short_period/simple_model_fit/short_lm_log_fitted.png"
    , width=800
    , height=600)
plot(fitted(m0_simple_log), rstandard(m0_simple_log)
     , xlab='Fitted values'
     , ylab='Standardized residuals')
dev.off()

png(file="../figures_report/application/short_period/simple_model_fit/short_lm_log_residuals_time_series.png"
    , width=800
    , height=600)
plot(data_short$date,rstandard(m0_simple_log)
     , ylab='standardized residuals'
     , xlab='date')
dev.off()


par(mfrow=c(1,1))
png(file="../figures_report/application/short_period/simple_model_fit/short_lm_log_day_of_week_full.png"
    , width=800
    , height=600)
plot(factor(data_short$weekday,levels=weekday_vec)
     , rstandard(m0_simple_log)
     , ylab=TeX("standardized residuals")
     , xlab='')
dev.off()

png(file="../figures_report/application/short_period/simple_model_fit/short_lm_log_day_of_week_trunc.png"
    , width=800
    , height=600)
plot(factor(data_short$weekday,levels=weekday_vec)
     , rstandard(m0_simple_log)
     , ylab='standardized residuals'
     , ylim=c(-3,3)
     , xlab='')
dev.off()


png(file="../figures_report/application/short_period/simple_model_fit/short_lm_log_acf.png"
    , width=800
    , height=600)
acf(rstandard(m0_simple_log)
    , main='')
dev.off()


table_out <- data.frame(confint(m0_simple_log)[,1])
table_out[,2] <- coef(m0_simple_log)
table_out[,3] <- confint(m0_simple_log)[,2]
xtable(table_out,digits=3)
sigma(m0_simple_log)^2
confint_sigma2(m0_simple_log)



### A more elaborate model...

# no improvements in the residuals nor the 
m0_larger_log <- lm(log(price) ~ days_from_start*factor(weekday)-1
                    , data = data_short)
summary(m0_larger_log)    

table_out <- data.frame(confint(m0_larger_log)[,1])
table_out[,2] <- coef(m0_larger_log)
table_out[,3] <- confint(m0_larger_log)[,2]
xtable(table_out,digits=6)
sigma(m0_larger_log)^2
confint_sigma2(m0_larger_log)

png(file="../figures_report/application/short_period/simple_model_fit/short_lm_log_intesection_diag.png"
    , width=800
    , height=600)
par(mfrow=c(2,2))
plot(m0_larger_log)
dev.off()

AIC(m0_larger_log)
drop1(m0_larger_log,test="F")
m1_larger <- update(m0,.~.-days_from_start:factor(weekday))
anova(m0_simple_log, m0_larger_log, test='F')




# a finder characterizations of the distributions.
#par(mfrow=c(2,2))
#for (weekday_i in weekday_vec){
#  data_subset = residuals(m0_simple_log)[data_short$weekday==weekday_i]
#  hist(data_subset[(data_subset>-1 & data_subset<1)]
#            , breaks=seq(-1,1,by=0.1)
#            , main=weekday_i
#            , xlab='')
#}


#### Remove the spikes
par(mfrow=c(1,1))
mask_outliers = which(abs(rstandard(m0_simple_log)) > 2.5)

# Indicate the spikes
cols = rep("steelblue",dim(data_short)[1])
cols[mask_outliers] <- "coral"
  png(file="../figures_report/application/short_period/simple_model_fit/short_lm_log_spikes_identified.png"
      , width=800
      , height=600)
  plot(data_short$date, data_short$price
       , ylab='EUR/MWh'
       , xlab='Date'
       , col=cols)
  legend('topright'
         , c('spikes', 'non spikes')
         , col=c('coral','steelblue')
         , pch=c(1,1)
         , bty='n')
  dev.off()
  length(mask_outliers)/dim(data_short)[1]
  
  legend(x = "topright",
         legend = c("EPMF", "Pois","Binom"), bty = 'n',
         col=c("coral", "darkGreen", "brown"),  pch=c(15,9,16),cex=1.3)
  # make plot of the spikes
  data_short$spike_vec <- 0
  data_short[mask_outliers,'spike_vec'] <- data_short[mask_outliers, 'price']
  
  
  png(file="../figures_report/application/short_period/simple_model_fit/short_lm_log_spikes_identified_bar.png"
      , width=800
      , height=600)
  barplot(data_short$spike_vec
          , ylab='EUR/MWh'
          , xlab='Date'
          , col=cols)
  dev.off()
  
  # re-estiamte model without the spikes
  m0_simple_log_rm_outlier <- lm(log(price) ~ days_from_start+factor(weekday)-1
                                 , data = data_short[-mask_outliers,])
  png(file="../figures_report/application/short_period/simple_model_fit/short_lm_w_log_rm_out_diag.png"
      , width=800
      , height=600)
  par(mfrow=c(2,2))
  plot(m0_simple_log_rm_outlier)
  dev.off()
  
  summary(m0_simple_log_rm_outlier)    
  AIC(m0_simple_log_rm_outlier)
  
  
  par(mfrow=c(1,1))
  png(file="../figures_report/application/short_period/simple_model_fit/short_lm_w_log_rm_out_res_fitted.png"
      , width=800
      , height=600)
  plot(fitted(m0_simple_log_rm_outlier), rstandard(m0_simple_log_rm_outlier)
       , ylab='standardized residuals'
       , xlab='date')
  dev.off()
  
  par(mfrow=c(1,1))
  png(file="../figures_report/application/short_period/simple_model_fit/short_lm_w_log_rm_out_res_ts.png"
      , width=800
      , height=600)
  plot(data_short[-mask_outliers,'date'], rstandard(m0_simple_log_rm_outlier)
       , ylab='standardized residuals'
       , xlab='date')
  dev.off()
  
  png(file="../figures_report/application/short_period/simple_model_fit/short_lm_w_log_rm_out_day_of_week.png"
      , width=800
      , height=600)
  plot(factor(data_short[-mask_outliers,'weekday'],levels=weekday_vec)
       ,rstandard(m0_simple_log_rm_outlier)
       , ylab='standardized residuals'
       , xlab='')
  dev.off()
  
  table_out <- data.frame(confint(m0_simple_log_rm_outlier)[,1])
  table_out[,2] <- coef(m0_simple_log_rm_outlier)
  table_out[,3] <- confint(m0_simple_log_rm_outlier)[,2]
  xtable(table_out,digits=3)
  sigma(m0_simple_log_rm_outlier)^2
  confint_sigma2(m0_simple_log_rm_outlier)
  
  #mask_outliers_second = which(abs(rstandard(m0_simple_log_rm_outlier)) > 3.3)
  #cols = rep("steelblue",dim(data_short)[1])
  #cols[mask_outliers_second] <- "coral"
  #plot(data_short$price, col=cols)
  
  m0_simple_log_rm_outlier_gls <- gls(log(price) ~ days_from_start+factor(weekday)-1
                                      , data = data_short[-mask_outliers,]
                                      , method='REML')
  summary(m0_simple_log_rm_outlier_gls)
  AIC(m0_simple_log_rm_outlier_gls)
  plot(m0_simple_log_rm_outlier_gls)
  plot(ACF(m0_simple_log_rm_outlier_gls),alpha=0.01)
  
  ## test if the same as before
  #m0_simple_log_rm_outlier_gls <- gls(log(price) ~ days_from_start+factor(weekday)-1
  #                                    , data = data_short[-mask_outliers,]
  #                                    , method='ML')
  #AIC()
  
  
  #data_short$weekday = factor(data_short$weekday)
  m0_gls_log_rm_outlier_varId <- gls(log(price) ~ days_from_start+factor(weekday)-1
                                     , weights = varIdent(form =~ 1|factor(weekday))
                                     , data = data_short[-mask_outliers,])
  
  summary(m0_gls_log_rm_outlier_varId)
  plot(m0_gls_log_rm_outlier_varId)
  par(mfrow=c(1,1))
  plot(fitted(m0_gls_log_rm_outlier_varId)
       , residuals(m0_gls_log_rm_outlier_varId, 'pearson'))
  png(file="../figures_report/application/short_period/simple_model_fit/short_lm_w_log_rm_out_varIndent_day_of_week.png"
      , width=800
      , height=600)
  plot(factor(data_short[-mask_outliers,'weekday'],levels=weekday_vec)
       ,residuals(m0_gls_log_rm_outlier_varId,'normalized')
       , ylab='standardized residuals'
       , xlab='')
  dev.off()
  intervals(m0_gls_log_rm_outlier_varId)
  
  anova(m0_simple_log_rm_outlier_gls, m0_gls_log_rm_outlier_varId)
  
  
  
  table_out <- data.frame(confint(m0_gls_log_rm_outlier_varId)[,1])
  table_out[,2] <- coef(m0_gls_log_rm_outlier_varId)
  table_out[,3] <- confint(m0_gls_log_rm_outlier_varId)[,2]
  xtable(table_out,digits=6)
  sigma(m0_gls_log_rm_outlier_varId)^2
  confint_sigma2(m0_gls_log_rm_outlier_varId)
  summary(m0_gls_log_rm_outlier_varId)
  
  
  
  m0_gls_log_rm_outlier_ar1 <- gls(log(price) ~ days_from_start+factor(weekday)-1
                                   #, weights = varIdent(form =~ 1|factor(weekday))
                                   , correlation = corAR1(form =~ as.numeric(days_from_start))
                                   , data = data_short[-mask_outliers,])
  
  summary(m0_gls_log_rm_outlier_ar1)
  
  
  table_out <- data.frame(confint(m0_gls_log_rm_outlier_ar1)[,1])
  table_out[,2] <- coef(m0_gls_log_rm_outlier_ar1)
  table_out[,3] <- confint(m0_gls_log_rm_outlier_ar1)[,2]
  xtable(table_out,digits=3)
  sigma(m0_gls_log_rm_outlier_ar1)^2
  m0_gls_log_rm_outlier_ar1$varBeta
  intervals(m0_gls_log_rm_outlier_ar1)
  
  confint_sigma2(m0_gls_log_rm_outlier_varId)
  
  
  AIC(m0_gls_log_rm_outlier_ar1)
  plot(m0_gls_log_rm_outlier_ar1)
  plot(ACF(m0_gls_log_rm_outlier_ar1, resType = 'n'),alpha=0.01)
  
  par(mfrow=c(1,1))
  png(file="../figures_report/application/short_period/simple_model_fit/short_lm_w_log_rm_out_corar1_acf.png"
      , width=800
      , height=600)
  acf(residuals(m0_gls_log_rm_outlier_ar1, 'normalized')
      , main='')
  dev.off()
  
  png(file="../figures_report/application/short_period/simple_model_fit/short_lm_w_log_rm_out_corar1_qqplot_norm_issue.png"
      , width=800
      , height=600)
  qqPlot(residuals(m0_gls_log_rm_outlier_ar1, 'normalized'))
  dev.off()
  
  
  png(file="../figures_report/application/short_period/simple_model_fit/short_lm_w_log_rm_out_corar1_diag.png"
      , width=800
      , height=600)
  par(mfrow=c(2,2))
  plot(fitted(m0_gls_log_rm_outlier_ar1)
       ,residuals(m0_gls_log_rm_outlier_ar1, 'normalized')
       , ylab='normalized residuals'
       , xlab='fitted' )
  abline(h=0,col='red')
  qqPlot(residuals(m0_gls_log_rm_outlier_ar1, 'normalized')
         , ylab="sample quantiles" )
  plot(factor(data_short[-mask_outliers,'weekday'],levels=weekday_vec)
       , residuals(m0_gls_log_rm_outlier_ar1, 'normalized')
       , ylab='normalized residuals'
       , xlab='date')
  plot(data_short[-mask_outliers,'date']
       , residuals(m0_gls_log_rm_outlier_ar1, 'normalized')
       , ylab='normalized residuals'
       , xlab='date')
  dev.off()
  
  
  par(mfrow=c(1,1))
  plot(data_short[-mask_outliers,'date']
       , residuals(m0_gls_log_rm_outlier_ar1, 'response')
       , ylab='normalized residuals'
       , xlab='date')
  
  
  
  
  qqnorm(residuals(m0_gls_log_rm_outlier_ar1, 'normalized'))
  plot(residuals(m0_gls_log_rm_outlier_ar1, 'normalized'))
  
  plot(factor(data_short[-mask_outliers,'weekday'],levels=weekday_vec)
       , residuals(m0_gls_log_rm_outlier_varId_ar1, 'normalized')
       , ylab='standardized residuals'
       , xlab='date')
  anova(m0_simple_log_rm_outlier_gls,m0_gls_log_rm_outlier_ar1)
  
  
  
  
  
  m0_gls_log_rm_outlier_varId_ar1 <- gls(log(price) ~ days_from_start+factor(weekday)-1
                                         , weights = varIdent(form =~ 1|factor(weekday))
                                         , correlation = corAR1(form =~ as.numeric(days_from_start))
                                         , data = data_short[-mask_outliers,])
  
  summary(m0_gls_log_rm_outlier_varId_ar1)
  AIC(m0_gls_log_rm_outlier_varId_ar1)
  
  anova(m0_gls_log_rm_outlier_ar1,m0_gls_log_rm_outlier_varId_ar1)
  
  
  mask_outliers_new = which(abs(residuals(m0_gls_log_rm_outlier_varId_ar1, 'normalized')) > 3)
  
  # Indicate the spikes
  cols = rep("steelblue",dim(data_short)[1])
  cols[mask_outliers] <- "coral"
    
  plot(data_short[-mask_outliers,'date'], data_short[-mask_outliers,'price']
       , ylab='EUR/MWh'
       , xlab='Date'
       , col=cols)
  legend('topright'
         , c('spikes', 'non spikes')
         , col=c('coral','steelblue')
         , pch=c(1,1)
         , bty='n')
  
  
  
  
  
  ###### export the relevant series 
  
  
  plot(data_short[-mask_outliers,'date']
       , exp(residuals(m0_gls_log_rm_outlier_varId_ar1, 'response'))
       , ylab='normalized residuals'
       , xlab='date'
       , type = "l"
       , col='steelblue')
  
  
  plot(data_short[-mask_outliers,'date']
       , data_short[-mask_outliers,'price'] - exp(fitted(m0_gls_log_rm_outlier_ar1))
       , ylab='residuals in original domain'
       , xlab='date'
       , type = "l"
       , col='steelblue')
  
  
  ##### export the residuals to fit  #####
  data_short$ou_estim_vec <- NA
  data_short[-mask_outliers,'ou_estim_vec'] <- data_short[-mask_outliers,'price'] - exp(fitted(m0_gls_log_rm_outlier_ar1))
  data_short$ou_estim_vec <- na.locf(na.locf(data_short[,'ou_estim_vec']), fromLast = TRUE)
  data_short <- data_short[-c(7)]
  
  data_short$fitted_vals <- exp(predict(m0_gls_log_rm_outlier_ar1,data_short))
  data_short$spike_vec_w_sign <- as.numeric(data_short$spike_vec)
  data_short[mask_outliers, 'spike_vec_w_sign'] <- data_short[mask_outliers, 'spike_vec_w_sign'] - data_short[mask_outliers, 'fitted_vals']
  
  
  write.csv(data_short, "../data/data_for_parm_estimation_short_series.csv", row.names=FALSE)
  
  
  
  #### Estimating Jump Diffusion Parameters #### 
  png(file="../figures_report/application/short_period/simple_model_fit/short_spikes_with_sign.png"
      , width=800
      , height=600)
  barplot(data_short[,'spike_vec_w_sign']
          , ylab='EUR/MWh'
          , xlab='Date'
          , col=cols
          , ylim=c(-50,180))
  grid()
  dev.off()
  