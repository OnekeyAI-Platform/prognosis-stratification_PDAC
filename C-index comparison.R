#C index comparison
getwd()
setwd("F:/data/results")
getwd()
rm(list = ls())
mydata <- read.csv("AUC_train.csv")
#install.packages("compareC")
library(compareC)
colnames(mydata)
# linear prediction
cox_fit1 <- coxph(Surv(duration, event) ~ Radiomics, data = mydata, x = TRUE, y = TRUE)
cox_fit2 <- coxph(Surv(duration, event) ~ DL, data = mydata, x = TRUE, y = TRUE)
lp1 <- predict(cox_fit1, type = "lp")  # Radiomics model
lp2 <- predict(cox_fit2, type = "lp")  # DL model
# C-index comprison
result_lp <- compareC(
  timeX = mydata$duration,
  statusX = mydata$event,
  scoreY = lp1,
  scoreZ = lp2
)








