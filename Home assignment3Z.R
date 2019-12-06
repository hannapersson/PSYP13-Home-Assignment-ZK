rm(list=ls(all=TRUE))
graphics.off()
# First clear the workspace

library(psych)
library(tidyverse)	
library(cAIC4)
library(r2glmm)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggplot2)
library(car)
library(influence.ME)
library(lm.beta)
# Load packages

data = read.csv("https://raw.githubusercontent.com/kekecsz/PSYP13_Data_analysis_class-2019/master/home_sample_3.csv")
# Open dataset

data1=data

View(data1)
describe(data1)
summary(data1)

ggplot() +	
  aes(x = data1$pain) +	
  geom_histogram( bins = 50)

ggplot() +	
  aes(x = data1$age) +	
  geom_histogram( bins = 40)

ggplot() +	
  aes(x = data1$STAI_trait) +	
  geom_histogram( bins = 40)	

ggplot() +	
  aes(x = data1$pain_cat) +	
  geom_histogram( bins = 40)

ggplot() +	
  aes(x = data1$mindfulness) +	
  geom_histogram( bins = 50)

ggplot() +	
  aes(x = data1$cortisol_serum) +	
  geom_histogram( bins = 50)
# Data diagnostics

(data1$sex = droplevels(replace(data1$sex, data1$sex == "Female", "female")))
# Change sex variable "Female" to "female" - assumed typo

data1 <- data1[-195,]
# Removing participant 195 since her mindful attention awareness scale ("mindfulness") score was higher than the maximum score (ranging from 1-6)

stdCoef.merMod <- function(object) {	
  sdy <- sd(getME(object,"y"))	
  sdx <- apply(getME(object,"X"), 2, sd)	
  sc <- fixef(object)*sdx/sdy	
  se.fixef <- coef(summary(object))[,"Std. Error"]	
  se <- se.fixef*sdx/sdy	
  return(data.frame(stdcoef=sc, stdse=se))	
}
# Function to extract standardized beta coefficients from linear mixed models 

(data1$hospital = factor(data1$hospital))
# Asign hospital as a grouping factor

mod_rnd_inter = lmer(pain ~ age + sex + STAI_trait + pain_cat + mindfulness + cortisol_serum + (1|hospital), data = data1)
# Creating linear random intercept model

###########################################

summary(mod_rnd_inter)
  
influence_observation = influence(mod_rnd_inter, obs = T)$alt.fixed
  
view(influence_observation)

data_plot_influence = as_tibble(influence_observation) %>% 	
  gather(colnames(influence_observation), value = coefficient, key = predictor)	

data_plot_influence %>% 	
  ggplot() +	
  aes(x = 1, y = coefficient, group = predictor) +	
  geom_violin() +	
  facet_wrap( ~ predictor, scales = "free")
# Looking for influential outliers using violin plot
# No outliers found
  
resid <- residuals(mod_rnd_inter)
hist(resid)
describe(resid)
qqnorm(resid)
shapiro.test(resid)
# Testing the normality of residuals using the Shapiro-Wilk test
# Tests not showing any abnormalities

yhat.2 <- fitted.values(object = mod_rnd_inter)
plot( x = yhat.2,
      y = data1$pain,
      xlab = "Fitted Values",
      ylab = "Observed Values")
# Checking the linearity of the relationship between predictors and outcome
# Test not showing any abnormalities

ggplot() +	
  aes(x = data1$age, y = resid) +	
  geom_point()
ggplot() +	
  aes(x = data1$sex, y = resid) +	
  geom_point()

ggplot() +	
  aes(x = data1$STAI_trait, y = resid) +	
  geom_point()	

ggplot() +	
  aes(x = data1$pain_cat, y = resid) +	
  geom_point()	

ggplot() +	
  aes(x = data1$cortisol_serum, y = resid) +	
  geom_point()	

ggplot() +	
  aes(x = data1$mindfulness, y = resid) +	
  geom_point()
# Scatterplot of the residuals and the fixed predictors separately
# Scatterplots showing a nonlinear relationship between residuals and each predictor separately

plot(mod_rnd_inter, arg = "pearson")	

homosced_mod = lm(resid^2 ~ hospital, data = data1)	
summary(homosced_mod)
# Checking the homogeneity of variances by viewing a plot with the stardardized residuals and predicted values
# Also checking homoscedasticity across clusters and running a significance test
# The plot did not show heteroscedasticity (a funnel shape would indicate this) 
# The model turned out non significant

pairs.panels(data1[,c("age", "STAI_trait", "pain_cat", "cortisol_serum","mindfulness")], col = "red", lm = T)	
# Checking the multicollinearity of the fixed predictors
# No multicollinearity found

###########################################

cAIC(mod_rnd_inter)$caic
# Calculating the random intercept model's conditional AIC

summary(mod_rnd_inter)
# Looking at the model's coefficients

r2beta(mod_rnd_inter, method = "nsj", data = data1)	
# Calculating the marginal R squared value based on the the fixed factors, not taking into account random effect terms (hospital)
# The 95% CI contains 0 in three variables: mindfulness, sexmale and STAI_trait. These variables are not significant predictors for this data

r.squaredGLMM(mod_rnd_inter)
# Calculating the marginal and the conditional R squared value

sm = summary(mod_rnd_inter)		
  sm_p_values = as.character(round(sm$coefficients[,"Pr(>|t|)"], 3))		
  sm_p_values[sm_p_values != "0" & sm_p_values != "1"] = substr(sm_p_values[sm_p_values != "0" & sm_p_values != "1"], 2, nchar(sm_p_values[sm_p_values != "0" & sm_p_values != "1"]))		
  sm_p_values[sm_p_values == "0"] = "<.001"		

coef_CI = suppressWarnings(confint(mod_rnd_inter))		
coef_CI
# Calculating 95% CI for each fixed variable

stdCoef.merMod(mod_rnd_inter)	
# Calculating standardised betas for each predictor

data2 = read.csv("https://raw.githubusercontent.com/kekecsz/PSYP13_Data_analysis_class-2019/master/home_sample_4.csv")
# Opening the second dataset

summary(data2)

ggplot() +	
  aes(x = data2$pain) +	
  geom_histogram( bins = 50)

ggplot() +	
  aes(x = data2$age) +	
  geom_histogram( bins = 40)

ggplot() +	
  aes(x = data2$STAI_trait) +	
  geom_histogram( bins = 40)	

ggplot() +	
  aes(x = data2$pain_cat) +	
  geom_histogram( bins = 40)

ggplot() +	
  aes(x = data2$mindfulness) +	
  geom_histogram( bins = 50)

ggplot() +	
  aes(x = data2$cortisol_serum) +	
  geom_histogram( bins = 50)
# Data diagnostics

influence_observation2 = influence(mod_rnd_inter, obs = T)$alt.fixed

view(influence_observation)

data_plot_influence = as_tibble(influence_observation2) %>% 	
  gather(colnames(influence_observation2), value = coefficient, key = predictor)	

data_plot_influence %>% 	
  ggplot() +	
  aes(x = 1, y = coefficient, group = predictor) +	
  geom_violin() +	
  facet_wrap( ~ predictor, scales = "free")
# Looking for influential outliers using violin plot
# No outliers found

(data2$hospital = factor(data2$hospital))
# Asign hospital as a grouping factor

pred_mod_rnd_inter <- predict(object = mod_rnd_inter, newdata = data2, allow.new.levels=TRUE)
# Creating predicted values for linear random intercept model

RSS_pred_mod_rnd_inter = sum((data2[,"pain"] - pred_mod_rnd_inter)^2)	
RSS_pred_mod_rnd_inter
# Calculating the model's residual sum of squares

mod_mean <- lm(pain ~ 1, data = data2)	

error_plotter <- function(mod, col = "black", x_var = NULL){	
  mod_vars = as.character(mod$call[2])	
  data = as.data.frame(eval(parse(text = as.character(mod$call[3]))))	
  y = substr(mod_vars, 1, as.numeric(gregexpr(pattern ='~',mod_vars))-2)	
  x = substr(mod_vars, as.numeric(gregexpr(pattern ='~',mod_vars))+2, nchar(mod_vars))	
  
  data$pred = predict(mod)	
  
  if(x == "1" & is.null(x_var)){x = "response_ID"	
  data$response_ID = 1:nrow(data)} else if(x == "1"){x = x_var}	
  
  plot(data[,y] ~ data[,x], ylab = y, xlab = x)	
  abline(mod)	
  
  for(i in 1:nrow(data)){	
    clip(min(data[,x]), max(data[,x]), min(data[i,c(y,"pred")]), max(data[i,c(y,"pred")]))	
    abline(v = data[i,x], lty = 2, col = col)	
  }	
  
}

error_plotter(mod_mean, col = "red", x_var = "pain")	

TSS = sum((data2$pain - predict(mod_mean))^2)	
TSS	
# Calculating the total sum of squared differences

R2 = 1-(RSS_pred_mod_rnd_inter/TSS)	
R2
# Calculating the model's R squared
# The model explains 0.307 of the variance in data 2
# Comparing the R squared for data 2 with the marginal R squared and conditional R squared for data 1, the model is better at explaining data 1 than data 2

r2beta(mod_rnd_inter, method = "nsj", data = data1)
# Calculating standardised beta for each predictor to see which predictor that is the most influential

best_rnd_slope = lmer(pain ~ cortisol_serum + (cortisol_serum|hospital), data = data1)
# Creating a random slope model that allows for both random intercept and random slope

pred_best_rnd_slope <- predict(object = best_rnd_slope, newdata = data1, allow.new.levels=TRUE)
# Creating predicted values for the new linear mixed effects model

best_rnd_slope %>% 		
  ggplot() +		
  aes(y = pain, x = cortisol_serum, group = hospital)+		
  geom_point(aes(color = hospital), size = 4) +		
  geom_line(color='red', aes(y=pred_best_rnd_slope, x=cortisol_serum))+		
  facet_wrap( ~ hospital, ncol = 2)
# Plotting the model for each hospital separately

Var_Fix_effect_age <- var(predict(lm(pain~age,data=data1)))
  Var_Fix_effect_age

Var_Fix_effect_sex <- var(predict(lm(pain~sex,data=data1)))
  Var_Fix_effect_sex

Var_Fix_effect_STAI <- var(predict(lm(pain~STAI_trait,data=data1)))
  Var_Fix_effect_STAI

Var_Fix_effect_pain <- var(predict(lm(pain~pain_cat,data=data1)))
  Var_Fix_effect_pain

Var_Fix_effect_mind <- var(predict(lm(pain~mindfulness,data=data1)))
  Var_Fix_effect_mind

Var_Fix_effect_cort <- var(predict(lm(pain~cortisol_serum,data=data1)))
  Var_Fix_effect_cort
# Calculating variance components for each fixed effect predictor in the random intercept model