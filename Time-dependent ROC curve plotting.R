#Time-dependent ROC curve
getwd()
setwd("C:/Users/lenovo/Desktop/results_ROc")
getwd()
mydata <- read.csv("All_model_AUC_train.csv")
library(survival)
library(survminer)
library(riskRegression)
names(mydata)
f1<-coxph(Surv(duration,event)~Radiomics,mydata,x=TRUE,y=TRUE)
f2<-coxph(Surv(duration,event)~DL,mydata,x=TRUE,y=TRUE)
f3<-coxph(Surv(duration,event)~DL25D,mydata,x=TRUE,y=TRUE)
f4<-coxph(Surv(duration,event)~Pre,mydata,x=TRUE,y=TRUE)
f5<-coxph(Surv(duration,event)~Post,mydata,x=TRUE,y=TRUE)
f6<-coxph(Surv(duration,event)~All_Post,mydata,x=TRUE,y=TRUE)
A1<-riskRegression::Score(list("Radiomics"=f1,"2D DL"=f2,"2.5D DL"=f3,"DL_Rad_FB"=f4,"DL_Ras_DB"=f5,"Combined"=f6),
                          formula=Surv(duration,event)~1,
                          data=mydata,
                          metrics="auc",
                          null.model=F, 
                          times=seq(5,60,2))

plotAUC(A1)
auc<-plotAUC(A1)
ggplot() +
  geom_line(
    data = auc,
    aes(times,AUC,group=model,col=model),
    linewidth = 0.9 
  ) +
  geom_hline(yintercept=1, linetype=2, linewidth=1) + 
  theme_classic() +
  labs(title = "Time-Dependent ROC Curves", x="Months", y="AUC",color="Model") +
  ylim(0.5, 1)+
  scale_x_continuous(breaks = seq(0, 60, by = 12)) +  
  theme(
    plot.title = element_text(size = 15),  
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 14),    
    axis.text.x = element_text(size = 15),
    axis.text.y = element_text(size = 15) 
  )

