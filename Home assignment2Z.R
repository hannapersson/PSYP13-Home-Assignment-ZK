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
# Checking data for abnormalities or deviations

data1$STAI_trait[data1$STAI_trait==3.5]=35.0
# Correcting STAI trait value of participant number 18 from "3.5" to "35.0" - assumed typo

data1$sex = as.numeric(data1$sex)
# Changing variable sex from character to numeric

data1$Participant <- c(1:160)
# Adding variable "Participant" as a "dummy variable" to predict in linear model

data1 <- data1[-59,]
# Removing participant number 59 since this participant was identified as an outlier in the first assignment

data2 = data1[,c("age","sex","pain","STAI_trait","pain_cat","mindfulness","cortisol_serum","weight","IQ","household_income")]

fullmodel <- lm(formula = pain ~ age + sex + STAI_trait + pain_cat + mindfulness + cortisol_serum + weight + IQ + household_income, data = data2)
# Linear model 1 with multiple predictors: age, sex, STAI trait, pain catastrophizing and cortisol serum

###########################################

lev = hat(model.matrix(fullmodel))
  plot(lev)
  data2[lev > .12,]
# Checking for multivariate outliers using leverage

N= nrow(data2)
  mahad=(N-1)*(lev-1 / N)
  tail(sort(mahad),5)
  order(mahad,decreasing=T)[c(5,4,3,2,1)]
# Checking for multivariate outliers using Mahalanobis distance, looking at the five most extreme scores
# Looking at Chi-square table and a p value that equals .001 and a degrees of freedom of 9, the cutoff score is 27.88 - no participant is considered a multivariate outlier

plot(x = fullmodel, which = 4)
# Checking for influental outliers using Cook's distance
# No outliers detected

resid <- residuals(fullmodel)
  hist(resid)
  describe(resid)
  qqnorm(resid)
  shapiro.test(resid)
# Testing the normality of residuals using the Shapiro-Wilk test
# Tests not showing any abnormalities
  
yhat.2 <- fitted.values(object = fullmodel)
  plot( x = yhat.2,
  y = data2$pain,
  xlab = "Fitted Values",
  ylab = "Observed Values")
# Checking the linearity of the relationship between predictors and outcomes
# Tests not showing any abnormalities
  
plot(x = fullmodel, which = 3)
ncvTest(fullmodel)
# Checking the homogeneity of variance
# Tests not showing any abnormalities

vif(mod = fullmodel)
# Checking for multicollinearity in regression model

###########################################

step(object = fullmodel,
  direction = "backward")
# start at the full model, then removing predictor variables one by one until the model produces the lowest AIC score

backwardmod <- lm(formula = pain ~ age + sex + pain_cat + mindfulness + cortisol_serum + weight, data = data2)
# Creating the backward model from the backward regression

###########################################

lev = hat(model.matrix(backwardmod))
  plot(lev)
  data2[lev > .10,]

N= nrow(data2)
  mahad=(N-1)*(lev-1 / N)
  tail(sort(mahad),5)
  order(mahad,decreasing=T)[c(5,4,3,2,1)]
# Looking at Chi-square table and a p value that equals .001 and a degrees of freedom of 6, the cutoff score is 22.46 - no participant is considered a multivariate outlier

plot(x = backwardmod, which = 4)

resid <- residuals(backwardmod)
  describe(resid)
  hist(resid)
  qqnorm(resid)
  shapiro.test(resid)

yhat.2 <- fitted.values(object = backwardmod)
plot( x = yhat.2,
      y = data2$pain,
      xlab = "Fitted Values",
      ylab = "Observed Values")

plot(x = backwardmod, which = 3)
ncvTest(backwardmod)

vif(backwardmod, which = 3)
# Rerunning data diagnostics for the backward model
# No outliers, abnormalities, or deviations found

###########################################

theorybasedmod <- lm(formula = pain ~ age + sex + STAI_trait + pain_cat + mindfulness + cortisol_serum, data = data2)
# Creating the theory-based model from assignment 1

summary(theorybasedmod)
summary(backwardmod)
# Model test statistics for the theory-based model and the backward model
# When looking at adjusted R-squared, the backward model explains more of the variance in the data than theory-based model (0.418 vs. 0.397)

AIC(theorybasedmod)	
AIC(backwardmod)
# Looking at the two models fit using the AIC function
# The difference between the AIC for theory-based model and backward model is bigger than 2
# The backward model has a smaller AIC value, meaning that it fits the data better than the theory-based model

coef_table = function(model){
  mod_sum = summary(backwardmod)
  mod_sum_p_values = as.character(round(mod_sum$coefficients[,4], 3))	
  mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"] = substr(mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"], 2, nchar(mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"]))	
  mod_sum_p_values[mod_sum_p_values == "0"] = "<.001"
  mod_sum_table = cbind(as.data.frame(round(cbind(coef(model), confint(model), c(0, lm.beta(model)$standardized.coefficients[c(2:length(model$coefficients))])), 2)), mod_sum_p_values)	
  names(mod_sum_table) = c("b", "95%CI lb", "95%CI ub", "Std.Beta", "p-value")	
  mod_sum_table["(Intercept)","Std.Beta"] = "0"	
  return(mod_sum_table)
}

sm_table = coef_table(backwardmod)
sm_table
# Creating a coefficient table for the backward model
# All variables except for STAI trait are significant when it comes to predicting pain

coef_table = function(model){
  mod_sum = summary(theorybasedmod)
  mod_sum_p_values = as.character(round(mod_sum$coefficients[,4], 3))	
  mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"] = substr(mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"], 2, nchar(mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values != "1"]))	
  mod_sum_p_values[mod_sum_p_values == "0"] = "<.001"
  mod_sum_table = cbind(as.data.frame(round(cbind(coef(model), confint(model), c(0, lm.beta(model)$standardized.coefficients[c(2:length(model$coefficients))])), 2)), mod_sum_p_values)	
  names(mod_sum_table) = c("b", "95%CI lb", "95%CI ub", "Std.Beta", "p-value")	
  mod_sum_table["(Intercept)","Std.Beta"] = "0"	
  return(mod_sum_table)
}

sm_table = coef_table(theorybasedmod)
sm_table

data3 = read.csv("https://raw.githubusercontent.com/kekecsz/PSYP13_Data_analysis_class-2019/master/home_sample_2.csv")

summary(data3)

data3$sex = as.numeric(data3$sex)
# Changing variable sex from character to numeric

data3$Participant <- c(1:160)
# Adding variable participant as a "dummy variable" to predict in linear model

Mod1 = lm(formula = Participant ~ age + sex + pain + STAI_trait + pain_cat + mindfulness + cortisol_serum + weight, data = data3)

###########################################

lev = hat(model.matrix(Mod1))
  plot(lev)
  data2[lev > .12,]
# Checking Mod 1 for multivariate outliers using leverage

N= nrow(data3)
  mahad=(N-1)*(lev-1 / N)
  tail(sort(mahad),5)
  order(mahad,decreasing=T)[c(5,4,3,2,1)]
# Checking Mod 1 for multivariate outliers using Mahalanobis distance, looking at the five most extreme scores
# Looking at Chi-square table and a p value that equals .001 and a degrees of freedom of 8, the cutoff score is 26.13 - no participant is considered a multivariate outlier

###########################################

pred_theorybasedmod <- predict(theorybasedmod, data3)	
pred_backwardmod <- predict(backwardmod, data3)	
# Creating predicted values for theory-based model and backward model

RSS_theorybasedmod = sum((data3[,"pain"] - pred_theorybasedmod)^2)	
RSS_backwardmod = sum((data3[,"pain"] - pred_backwardmod)^2)	
RSS_theorybasedmod	
RSS_backwardmod
# Checking the sum of squared residuals for the theory-based model and the backward model
# The backward has more residual error than the theory-based model (231.78 compared to 227.63), meaning that the theory-based model is a better predictor for pain based on the new data set