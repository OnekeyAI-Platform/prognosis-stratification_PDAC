##DCA curve and calibration
getwd()
setwd("C:/Users/lenovo/Desktop/data/results")
getwd()
library(timeROC)
library(ggDCA)
library(rms)
library(survival)
library(ggsci)
mydata <- read.csv("All_model_prediction.csv")
colnames(mydata)
dd <- datadist(mydata)
options(datadist = "dd")
#-----------------------------------DCA curve-----------------------------------
f1<-coxph(Surv(duration,event)~Radiomics,mydata,x=TRUE,y=TRUE)
f2<-coxph(Surv(duration,event)~DL,mydata,x=TRUE,y=TRUE)
f3<-coxph(Surv(duration,event)~DL25D,mydata,x=TRUE,y=TRUE)
f4<-coxph(Surv(duration,event)~Pre,mydata,x=TRUE,y=TRUE)
f5<-coxph(Surv(duration,event)~Post,mydata,x=TRUE,y=TRUE)
f6<-coxph(Surv(duration,event)~All_Post,mydata,x=TRUE,y=TRUE)
df1 <- ggDCA::dca(f1, f2, f3, f4, f5, f6,
                  times = 12 
)
ggplot(df1,linetype = F)+
  scale_color_lancet(name="Model Type",labels=c("Radiomics","2D DL","2.5D DL","DL_Rad_FB","DL_Rad_DB","Combined"))+
  theme_bw(base_size = 15)+
  theme(legend.position = c(0.8,0.75),
        legend.background = element_blank(),
        legend.title = element_text(size = 15),  
        legend.text = element_text(size = 13),   
        axis.title = element_text(size = 16),    
        axis.text = element_text(size = 15)      
  )
#-------------------------------------------------------------------------------
#-----------------------------------calibration---------------------------------
library(dplyr)
library(tidyr)
mydata <- read.csv("All_model_AUC_prediction.csv")
dd <- datadist(mydata)
options(datadist = "dd")
f1 <- cph(Surv(duration, event) ~ Radiomics,
               data = mydata, x=T,y=T,surv = T,
               time.inc = 12)
f2 <- cph(Surv(duration, event) ~ DL,
          data = mydata, x=T,y=T,surv = T,
          time.inc = 12)
f3 <- cph(Surv(duration, event) ~ DL25D,
          data = mydata, x=T,y=T,surv = T,
          time.inc = 12)
f4 <- cph(Surv(duration, event) ~ Pre,
          data = mydata, x=T,y=T,surv = T,
          time.inc = 12)
f5 <- cph(Surv(duration, event) ~ Post,
          data = mydata, x=T,y=T,surv = T,
          time.inc = 12)
f6 <- cph(Surv(duration, event) ~ All_Post,
          data = mydata, x=T,y=T,surv = T,
          time.inc = 12)
cal1 <- calibrate(f1, cmethod="KM", method="boot",u=12,m=9,B=500)
cal2 <- calibrate(f2, cmethod="KM", method="boot",u=12,m=9,B=500)
cal3 <- calibrate(f3, cmethod="KM", method="boot",u=12,m=9,B=500)
cal4 <- calibrate(f4, cmethod="KM", method="boot",u=12,m=9,B=500)
cal5 <- calibrate(f5, cmethod="KM", method="boot",u=12,m=9,B=500)
cal6 <- calibrate(f6, cmethod="KM", method="boot",u=12,m=9,B=500)

plot(cal1,lwd = 2,lty = 0,errbar.col = c("#008B45FF"),
     xlim = c(0,1),ylim= c(0,1),
     xlab = "Nomogram-prediced OS (%)",ylab = "Observed OS (%)",
     col = c("#008B45FF"),
     cex.lab=1.6,cex.axis=1.6, cex.main=1.2, cex.sub=0.6)
lines(cal1[,c('mean.predicted',"KM")],
      type = 'b', lwd = 2, col = c("#008B45FF"), pch = 16)

plot(cal2,lwd = 2,lty = 0,errbar.col = c("#F1A340FF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#F1A340FF"),add = T)
lines(cal2[,c('mean.predicted',"KM")],
      type = 'b', lwd = 2, col = c("#F1A340FF"), pch = 16) #F1A340FF金棕

plot(cal3,lwd = 2,lty = 0,errbar.col = c("#B2182B"),
     xlim = c(0,1),ylim= c(0,1),col = c("#B2182B"),add = T)
lines(cal3[,c('mean.predicted',"KM")],
      type = 'b', lwd = 2, col = c("#B2182B"), pch = 16)

plot(cal4,lwd = 2,lty = 0,errbar.col = c("#8E44ADFF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#8E44ADFF"),add = T)
lines(cal4[,c('mean.predicted',"KM")],
      type = 'b', lwd = 2, col = c("#8E44ADFF"), pch = 16)

plot(cal5,lwd = 2,lty = 1,errbar.col = c("#2166AC"),
     xlim = c(0,1),ylim= c(0,1),col = c("#2166AC"),add = T)
lines(cal5[,c('mean.predicted',"KM")],
      type = 'b', lwd = 2, col = c("#2166AC"), pch = 16)

plot(cal6,lwd = 2,lty = 0,errbar.col = c("#E64B35FF"),
     xlim = c(0,1),ylim= c(0,1),col = c("#E64B35FF"),add = T)
lines(cal6[,c('mean.predicted',"KM")],
      type = 'b', lwd = 2, col = c("#E64B35FF"), pch = 16)
box(lwd = 2) 
abline(0,1,lty = 3, 
       lwd = 2, 
       col = "grey70" 
)



