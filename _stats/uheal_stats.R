library(lattice)  # graphical display (xyplot)
library(lme4)     # subject specific models (GLMM)
library(psych)    # OPTIONAL: scatterplots with effects
library(ggplot2)  # to make other plots
library(dplyr)
library(LMMstar)
uheal_data <- read.csv("uheal_data2.csv", header=TRUE, stringsAsFactors=FALSE)


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
uheal_data$CP <- factor(uheal_data$CP)

summary(uheal_data)
str(uheal_data)

# AP amplitude
fm.ap <- lm(AP_amp ~-1+PTA_hf*Age*sex, data=uheal_data)

summary(fm.ap)
xyplot(AP_amp~PTA_hf|gender_fact, group=gender, data=uheal_data)
xyplot(AP_amp ~ Age| gender_fact, data=uheal_data,
       panel = function(x, y, ...) {
         panel.xyplot(x, y, ...)
         panel.lmline(x, y, col = "red")
       })
plot(fitted(fm.ap),residuals(fm.ap))
abline(h=0)
qqnorm(residuals(fm.ap))
qqline(residuals(fm.ap), col = "red") 

# Extract estimates and confidence intervals
confint(fm.ap)


# FFR SNR
fm2 <- lm(FFR_SNR_sig~ -1 + PTA_hf*Age*gender,uheal_data)

summary(fm2)
xyplot(FFR_SNR_sig~Age|gender_fact, group=gender, data=uheal_data)
xyplot(FFR_SNR_sig~ Age| gender_fact, data=uheal_data,
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

# memr
fm3 <- lm(memr~-1+ Age*gender*PTA_hf*PTA_lf,uheal_data)

summary(fm3)
qqnorm(residuals(fm3))
hist(residuals(fm3))
