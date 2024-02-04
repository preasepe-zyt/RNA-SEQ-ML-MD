library(pROC)
library(ggplot2)
#建立曲线
roc -> roc_auc
ppi_auc <- c()
for (i in 1:3) {
    ppi_auc[i] <- roc(roc_auc[,1],roc_auc[,(i+1)]) %>% auc()
}

#计算FULL auc
auc(rocobj1)
auc(rocobj2)
auc(rocobj3)
#绘制曲线
plot(rocobj1)
#其他参数美化
plot(rocobj1, print.auc=TRUE, auc.polygon=TRUE, 
     grid=c(0.1, 0.2), 
     grid.col=c("green", "red"), 
     max.auc.polygon=TRUE,
     auc.polygon.col="skyblue",
     print.thres=TRUE)
#计算partial AUC选择关注(sp或se)一定范围的数据
plot(rocobj1, print.auc=TRUE, auc.polygon=TRUE, 
     partial.auc=c(1, 0.8),partial.auc.focus="sp", 
     grid=c(0.1, 0.2),grid.col=c("green", "red"),
     max.auc.polygon=TRUE, 
     auc.polygon.col="skyblue",
     print.thres=TRUE, reuse.auc=FALSE)
#方法2
library(ROCR)
data(ROCR.simple)
df <- data.frame(ROCR.simple)
#计算tpr，fpr
pred<-prediction(df$predictions, df$labels)
perf<- performance(pred,"tpr","fpr")
perf
plot(perf,colorize=TRUE)









