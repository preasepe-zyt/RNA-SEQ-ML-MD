library(survival)
library(survminer)
library(rms)
library(rmda)
library(ggDCA)
library(ggplot2)
library(caret)
as.data.frame(fa4) -> data
ddist <- datadist(as.data.frame(roc_randomforest))
options(datadist='ddist')
f <- lrm(group_list ~ ALOX5+C5AR1+VDR+HMGCR, data=roc_randomforest,x=TRUE,y=TRUE) 
summary(f)   ##也能用此函数看具体模型情况，模型的系数，置信区间等
###  nomogram

nomogram <- nomogram(f,fun=function(x)1/(1+exp(-x)), ##逻辑回归计算公式
                     fun.at = c(0.01,0.6,0.99),#风险轴刻度
                     funlabel = "Risk of Depression", #风险轴便签
                     lp=F,  ##是否显示系数轴
                     conf.int = F, ##每个得分的置信度区间，用横线表示，横线越长置信度越
                     abbrev = F#是否用简称代表因子变量,
                     
)
pdf("nomogram.pdf",width = 10,height = 6,family = "Arial")
plot(nomogram,xfrac=.45)
dev.off()
#dca决策曲线
pdf("DCA.pdf",width = 10,height = 8,family = "Arial")
dca_lrm <- dca(f ,model.names = "Nomogram")
ggplot(dca_lrm, lwd = 2)
dev.off()
#Calibration curves of nomogram
cal2 <- calibrate(f, method='boot', B=5000)
plot(cal2,xlim=c(0,1.0),ylim=c(0,1.0), lwd = 2)
