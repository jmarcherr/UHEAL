library(lattice)# graphical display (xyplot)
library(lme4)     # subject specific models (GLMM)
library(psych)    # OPTIONAL: scatterplots with effects
library(ggplot2)  # to make other plots
library(dplyr)
library(LMMstar)
library(plotrix)               # Load plotrix
#library(lfe)

uheal_data <- read.csv("uheal_data.csv", header=TRUE, stringsAsFactors=FALSE)
uheal_data <- filter(uheal_data, CP_new == 0)

# Make categorical variables into factors:
uheal_data$subid <- factor(uheal_data$subid)
# gender factor
uheal_data$gender_fact <-factor(uheal_data$gender)
levels(uheal_data$gender_fact) <- c("Female", "Male")
# translate to sex
uheal_data$sex <-uheal_data$gender
uheal_data$sex_fact <-factor(uheal_data$sex)
levels(uheal_data$sex_fact) <- c("Female", "Male")
# NH/HI
uheal_data$CP <- factor(uheal_data$CP_new)

summary(uheal_data)
str(uheal_data)

# AP amplitude
fm.ap <- lm(AP_amp_pm ~-1+PTA_hf*PTA_lf_new*Age*sex, data=uheal_data)

summary(fm.ap)
xyplot(AP_amp_pm~PTA_hf|gender_fact, group=gender, data=uheal_data)
xyplot(AP_amp_pm ~ Age| gender_fact, data=uheal_data,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.lmline(x, y, col = "red")
       })
plot(fitted(fm.ap),residuals(fm.ap))
abline(h=0)
qqnorm(residuals(fm.ap))
qqline(residuals(fm.ap), col = "red") 
hist(residuals(fm.ap))

# Extract estimates and confidence intervals
cf <-confint(fm.ap,level = 0.99)


####### FFR_sig ##############
uheal_data_FFR = filter(uheal_data,FFR_sig == 1)
# FFR SNR
fm2 <- lm(FFR_SNR~ -1 + PTA_hf*PTA_lf_new*Age*sex,uheal_data_FFR)

summary(fm2)
xyplot(FFR_SNR~Age|gender_fact, group=gender, data=uheal_data_FFR)
xyplot(FFR_SNR~ Age| gender_fact, data=uheal_data_FFR,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.lmline(x, y, col = "red")
       })
plot(fitted(fm.ap),residuals(fm.ap))
abline(h=0)
qqnorm(residuals(fm.ap))
qqline(residuals(fm.ap), col = "red") 
qqnorm(residuals(fm2))
hist(residuals(fm2))
# Extract estimates and confidence intervals
confint(fm2,level = 0.99)

# memr
fm3 <- lm(memr~-1+ PTA_hf*PTA_lf_new*Age*sex,uheal_data)

summary(fm3)
qqnorm(residuals(fm3))
hist(residuals(fm3))
# Extract estimates and confidence intervals
confint(fm3,level = 0.99)


#### fixed effects between measures
# AP ~ FFR
fm4 <- lm(AP_amp_pm~FFR_SNR,uheal_data_FFR)
summary(fm4)

# AP ~ memr
fm5 <- lm(AP_amp_pm~memr,uheal_data)
summary(fm5)

# PTA
fm6 <- lm(PTA_hf~AP_amp_pm*Age,uheal_data)
summary(fm6)

############ sørens approach

library(car)
# simulate some data, we can discuss details if you like
# set.seed(4)
# n=120
# age = rnorm(n)*12+50
# gender = round(runif(n))
# 
# E = mvtnorm::rmvnorm(n=n,mean=c(0,0.1,0.5,0.05),sigma=toeplitz(seq(1,0.1,l=4)))
# X = cbind(Age,sex,data = uheal_data)
# 
# 
# B = cbind(c(-0.1,0),c(-0.2,0),c(0,0),c(0,1))
# 
# Y = X%*%B + E
# 
# dat = data.frame(EFR=Y[,1],
#                  ABR=Y[,2],
#                  MEMR=Y[,3],
#                  AP=Y[,4],
#                  age=X[,1],
#                  sex=X[,2])
uheal_data$FFR_SNR[uheal_data$FFR_sig==0]<-NaN
#png(file="scatterplot.png")
scatterplotMatrix(~ FFR_SNR + AP_amp_pm + memr
                  + PTA_hf + Age | sex_fact,
                  data=uheal_data, smooth=FALSE, regLine=FALSE, ellipse=FALSE, diagonal=TRUE, legend=list(coords="bottomleft"))
#dev.off()

  # Approach 1: univariate, with separate models, populate p value data frame
p_val_age = NULL
for (c in c('FFR_SNR','memr','AP_amp_pm','WV_amp_pm','SP_amp','PTA_lf','PTA_hf')){
  print(summary(lm(paste(c,"~ Age+sex"), data=uheal_data)))
  p_val_age = cbind(p_val_age,summary(lm(paste(c,"~ Age+sex"), data=uheal_data))$coefficients[2,4])
}
print(p.adjust(p_val_age,method='bonferroni'))


# Approach 2:_ 
mod.age = lm(cbind(FFR_SNR, AP_amp_pm,WV_amp_pm, memr, PTA_hf)
             ~ Age+sex, data=uheal_data)
print(mod.age)

manova.age = Anova(mod.age)
summary(manova.age)

# AP~ FFR 
p_val_FFR = NULL
for (c in c('AP_amp_pm','AP_amp_pm+Age')){
  print(summary(lm(paste(c,"~ FFR_SNR"), data=uheal_data)))
  p_val_FFR = cbind(p_val_FFR,summary(lm(paste(c,"~ FFR_SNR"), data=uheal_data))$coefficients[2,4])
}
print(p.adjust(p_val_FFR,method='bonferroni'))


## Try FFR ~ AP + Age

lmAP=lm('FFR_SNR~PTA_hf+AP_amp_pm+Age',data=uheal_data)
summary(lmAP)

res <- c(residuals(lmAP))
FFR <- uheal_data$FFR_SNR[!is.nan(uheal_data$FFR_SNR) & !is.nan(uheal_data$AP_amp_pm)]
AP <- uheal_data$AP_amp_pm[!is.nan(uheal_data$AP_amp_pm) & !is.nan(uheal_data$FFR_SNR)]

