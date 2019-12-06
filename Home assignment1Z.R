rm(list=ls(all=TRUE))
graphics.off()
# First clear the workspace

library(psych)
library(lm.beta)
library(rgl)
library(tidyverse)
library(gridExtra)
library(lsr)
library(MVA)
library(ggplot2)
library(car)
# Load packages

data = read.csv("https://raw.githubusercontent.com/kekecsz/PSYP13_Data_analysis_class-2019/master/home_sample_1.csv ")
# Open dataset

data1 = data
View(data)
describe(data)
summary(data)

ggplot() +	
  aes(x = data$pain) +	
  geom_histogram( bins = 50)

ggplot() +	
  aes(x = data$age) +	
  geom_histogram( bins = 40)

ggplot() +	
  aes(x = data$STAI_trait) +	
  geom_histogram( bins = 40)	

ggplot() +	
  aes(x = data$pain_cat) +	
  geom_histogram( bins = 40)

ggplot() +	
  aes(x = data$mindfulness) +	
  geom_histogram( bins = 50)

ggplot() +	
  aes(x = data$cortisol_serum) +	
  geom_histogram( bins = 50)

ggplot() +	
  aes(x = data$cortisol_saliva) +	
  geom_histogram( bins = 50)
# Histograms for variables pain, age, STAI trait, pain catastrophizing, mindfulness and cortisol measures
# Checking data for any abnormalities or deviations

data1$STAI_trait[data1$STAI_trait==3.5]=35.0
# Correcting STAI trait value of participant number 18 from "3.5" to "35.0", assumed typo

data1$sex = as.numeric(data1$sex)
# Changing variable sex from character to numeric

data2 = data1[,c("age","sex","pain","STAI_trait","pain_cat","mindfulness","cortisol_serum","cortisol_saliva")]

###########################################
  
model1 <- lm(pain ~ age + sex, data = data2)
# Linear model 1 with multiple predictors: age and sex

summary(model1)
# View linear model 1

model2 <- lm(formula = pain ~ age + sex + STAI_trait + pain_cat + mindfulness + cortisol_serum + cortisol_saliva, data = data2)
# Linear model 2 with multiple predictors: age, sex, STAI trait, pain catastrophizing and cortisol measures

summary(model2)
# View linear model 2

###########################################

lev = hat(model.matrix(model1))
  plot(lev)
  data2[lev > .05,]
# Checking for multivariate outliers in linear model 1 using leverage

N= nrow(data2)
  mahad=(N-1)*(lev-1 / N)
  tail(sort(mahad),5)
  order(mahad,decreasing=T)[c(5,4,3,2,1)]
# Checking for multivariate outliers in linear model 1 using Mahalanobis distance, looking at the five most extreme scores
# Looking at Chi-square table and a p value that equals .001 and a degrees of freedom of 2, the cutoff score is 13.82 - no participant is considered a multivariate outlier

lev = hat(model.matrix(model2))
  plot(lev)
  data2[lev > .10,]
# Checking for multivariate outliers in linear model 2 using leverage

N= nrow(data2)
  mahad=(N-1)*(lev-1 / N)
  tail(sort(mahad),5)
  order(mahad,decreasing=T)[c(5,4,3,2,1)]
# Looking at Chi-square table and a p value that equals .001 and a degrees of freedom of 7, the cutoff score is 24.32 - participant number 59 is considered a multivariate outlier and will thus be excluded from the dataset

data2 <- data2[-59,]
# Removing participant number 59

model1 = update(model1)
model2 = update(model2)
# Updating the models after the participant is removed

###########################################

lev = hat(model.matrix(model1))
plot(lev)
data2[lev > .05,]

N= nrow(data2)
mahad=(N-1)*(lev-1 / N)
tail(sort(mahad),5)
order(mahad,decreasing=T)[c(5,4,3,2,1)]
# Looking at Chi-square table and a p value that equals .001 and a degrees of freedom of 2, the cutoff score is 13.82 - no participant is considered a multivariate outlier

lev = hat(model.matrix(model2))
plot(lev)
data2[lev > .10,]

N= nrow(data2)
mahad=(N-1)*(lev-1 / N)
tail(sort(mahad),5)
order(mahad,decreasing=T)[c(5,4,3,2,1)]
# Looking at Chi-square table and a p value that equals .001 and a degrees of freedom of 7, the cutoff score is 24.32 - no participant is considered a multivariate outlier
# Rerunning tests for multivariate outliers - no outliers found

###########################################
  
plot(x = model1, which = 4)
plot(x = model2, which = 4)
# Checking regression model 1 and 2 for influental outliers using Cook's distance
# No outliers detected

resid <- residuals(model1)
  hist(resid)
  describe(resid)
  qqnorm(resid)
  shapiro.test(resid)
  
resid <- residuals(model2)
  hist(resid)
  describe(resid)	
  qqnorm(resid)
  shapiro.test(resid)
# Testing regression model 1 and 2 for the normality of residuals using the Shapiro-Wilk test
# Tests not showing any abnormalities

yhat.2 <- fitted.values(object = model1)
  plot( x = yhat.2,
  y = data2$pain,
  xlab = "Fitted Values",
  ylab = "Observed Values")

yhat.2 <- fitted.values(object = model2)
  plot( x = yhat.2,
  y = data2$pain,
  xlab = "Fitted Values",
  ylab = "Observed Values")
# Checking the linearity of the relationship between predictors and outcomes in regression model 1 and 2
# Tests not showing any abnormalities

plot(x = model1, which = 3)
plot(x = model2, which = 3)
ncvTest(model1)
ncvTest(model2)
# Checking the homogeneity of variance for regression model 1 and 2
# Tests not showing any abnormalities

vif(mod = model1)
vif(mod = model2)
# Checking for multicollinearity in regression model 1 and 2
# Regression model 2 showing two fairly large correlations between predictor variables cotrisol serum and cortisol saliva
# Variable cortisol saliva has a higher value than cortisol serum and will thus be excluded from a third regression model

###########################################

model3 <- lm(formula = pain ~ age + sex + STAI_trait + pain_cat + mindfulness + cortisol_serum, data = data2)
# Creating an updated second model without the variable salivary cortisol

#####################################################

lev = hat(model.matrix(model3))
  plot(lev)
  data2[lev > .10,]

N= nrow(data2)
  mahad=(N-1)*(lev-1 / N)
  tail(sort(mahad),5)
  order(mahad,decreasing=T)[c(5,4,3,2,1)]

plot(x = model3, which = 4)

resid <- residuals(model3)
  hist(resid)
  describe(resid)
  qqnorm(resid)
  shapiro.test(resid)

yhat.2 <- fitted.values(object = model3)
  plot( x = yhat.2,
  y = data2$pain,
  xlab = "Fitted Values",
  ylab = "Observed Values")

plot(x = model3, which = 3)
ncvTest(model3)

vif(model3, which = 3)
# Rerunning model diagnostics for the updated model
# No outliers, abnormalities, or deviations found

#####################################################

summary(model1)
summary(model3)
# Model test statistics for model 1 and the updated model 2
# When looking at adjusted R-squared, the updaed model 2 explains more of the variance in the data than model 1

AIC(model1)	
AIC(model3)
# Looking at the two models fit using the AIC function
# The difference between the AIC for model 1 and the updated model 2 are bigger than 2
# The updated model 2 has a smaller AIC value, meaning that it fits the data better than model 1

anova(model1, model3)	
# Looking at the two models fit, based on residual error and degrees of freedom, using the anova function
# The updated model 2 is significantly better at predicting pain than model 1

coef_table = function(model){
  mod_sum = summary(model3)
  mod_sum_p_values = as.character(round(mod_sum$coefficients[,4], 3))	
  mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"] = substr(mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"], 2, nchar(mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"]))	
  mod_sum_p_values[mod_sum_p_values == "0"] = "<.001"
  mod_sum_table = cbind(as.data.frame(round(cbind(coef(model), confint(model), c(0, lm.beta(model)$standardized.coefficients[c(2:length(model$coefficients))])), 2)), mod_sum_p_values)	
  names(mod_sum_table) = c("b", "95%CI lb", "95%CI ub", "Std.Beta", "p-value")	
  mod_sum_table["(Intercept)","Std.Beta"] = "0"	
  return(mod_sum_table)
}

sm_table = coef_table(model3)
sm_table
# Creating a coefficient table for the updated second model
# All variables except for STAI trait are significant when it comes to predicting pain in this dataset

coef_table = function(model){
  mod_sum = summary(model1)
  mod_sum_p_values = as.character(round(mod_sum$coefficients[,4], 3))	
  mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"] = substr(mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"], 2, nchar(mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"]))	
  mod_sum_p_values[mod_sum_p_values == "0"] = "<.001"
  mod_sum_table = cbind(as.data.frame(round(cbind(coef(model), confint(model), c(0, lm.beta(model)$standardized.coefficients[c(2:length(model$coefficients))])), 2)), mod_sum_p_values)	
  names(mod_sum_table) = c("b", "95%CI lb", "95%CI ub", "Std.Beta", "p-value")	
  mod_sum_table["(Intercept)","Std.Beta"] = "0"	
  return(mod_sum_table)
}

sm_table = coef_table(model1)
sm_table
# Creating a coefficient table for the first model
# All variables are significant when it comes to predicing pain in this dataset