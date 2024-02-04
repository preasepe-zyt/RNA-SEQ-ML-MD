
#引用包
library(pROC)
library(extrafont)
library(tidyverse)
loadfonts(device = "pdf")
#导入字体进入pdf
#showtext包可给定字体文件，加载到 R环境中，生成新的字体家族名字，后期调用这个名字设定字体，并且支持中文写入pdf不乱码
library(showtext)
showtext_auto(enable=TRUE)
#train
ids_exprs %>% mutate(gene=rownames(ids_exprs)) -> data
data[intersect(data$gene,final_genes$x),] %>% select(c(1:22)) %>% t() -> train
cbind(group_list,train) %>% as.data.frame() -> train
write.csv(train,"train.txt")
train$group_list -> group_list

#pls-da
data[intersect(data$gene,c("MMP9","ALOX5","C5AR1","MAPK14")),] %>% select(c(1:22)) %>% t() -> pls_da
cbind(group_list,pls_da) %>% as.data.frame() -> roc_pls_da
pls_da_model <- glm(group_list~ MMP9+ALOX5+C5AR1+MAPK14, roc_pls_da,family = "binomial",control=list(maxit=100))

#lasso
data[intersect(data$gene,gene_index_lasso$geneids),] %>% select(c(1:22)) %>% t() -> lasso
cbind(group_list,lasso) %>% as.data.frame() -> roc_lasso
colnames(roc_lasso)
lasso_model <- glm(group_list~ ALOX5+NR1H4+MALT1, roc_lasso,family = "binomial",control=list(maxit=100))
lasso_model <- glm(group_list~ NR1H4+MALT1, roc_lasso,family = "binomial",control=list(maxit=100))
#SVM_RFE
data[intersect(data$gene,SVM_RFE_gene$VDR),] %>% select(c(1:22)) %>% t() -> SVM_RFE
cbind(group_list,SVM_RFE) %>% as.data.frame() -> roc_SVM_RFE
colnames(roc_SVM_RFE)
SVM_RFE_model <- glm(group_list~ ALOX5+C5AR1+HMGCR, roc_SVM_RFE,family = "binomial")

#randomforest
data[intersect(data$gene,gene_index_randomforest$Feature),] %>% select(c(1:22)) %>% t() -> randomforest
cbind(group_list,randomforest) %>% as.data.frame() -> roc_randomforest
colnames(roc_randomforest)
randomforest_model <- glm(group_list~ ALOX5+C5AR1+VDR+HMGCR, roc_randomforest,family = "binomial")

#GBM
data[intersect(data$gene,GBM$var),] %>% select(c(1:22)) %>% t() -> GBM
cbind(group_list,GBM) %>% as.data.frame() -> roc_GBM
colnames(roc_GBM)
GBM_model <- glm(group_list~ ALOX5+C5AR1+VDR+NR1H4+HMGCR+MMP9+MALT1+MAPK14+ALOX5AP,roc_GBM,family = "binomial")
GBM_model <- glm(group_list~ C5AR1+VDR+NR1H4+HMGCR+MMP9+MALT1+MAPK14+ALOX5AP,roc_GBM,family = "binomial",control=list(maxit=100))
#xboost
data[intersect(data$gene,xboost$Feature),] %>% select(c(1:22)) %>% t() -> xboost
cbind(group_list,xboost) %>% as.data.frame() -> roc_xboost
colnames(roc_xboost)
xboost_model <- glm(group_list~ ALOX5+C5AR1+VDR+NR1H4+HMGCR+MAPK14,roc_xboost,family = "binomial",control=list(maxit=100))
xboost_model <- glm(group_list~ C5AR1+VDR+NR1H4+HMGCR+MAPK14,roc_xboost,family = "binomial",control=list(maxit=100))
#ridge
data[intersect(data$gene,gene_index_ridge$V1),] %>% select(c(1:22)) %>% t() -> ridge
cbind(group_list,ridge) %>% as.data.frame() -> roc_ridge
colnames(roc_ridge)
ridge_model <- glm(group_list~ ALOX5+C5AR1+VDR+NR1H4+HMGCR+MMP9+MALT1+MAPK14+ALOX5AP,roc_ridge,family = "binomial",control=list(maxit=100))
ridge_model <- glm(group_list~ C5AR1+VDR+NR1H4+HMGCR+MMP9+MALT1+MAPK14+ALOX5AP,roc_ridge,family = "binomial",control=list(maxit=100))

#ridge
data[intersect(data$gene,gene_index_ridge$V1),] %>% select(c(1:22)) %>% t() -> ridge
cbind(group_list,ridge) %>% as.data.frame() -> roc_ridge
colnames(roc_ridge)
ridge_model <- glm(group_list~ ALOX5+C5AR1+VDR+NR1H4+HMGCR+MMP9+MALT1+MAPK14+ALOX5AP,roc_ridge,family = "binomial",control=list(maxit=100))
ridge_model <- glm(group_list~ C5AR1+VDR+NR1H4+HMGCR+MMP9+MALT1+MAPK14+ALOX5AP,roc_ridge,family = "binomial",control=list(maxit=100))
#ELASTIC NET 	
data[intersect(data$gene,c("MMP9","C5AR1","MAPK14","ALOX5")),] %>% select(c(1:22)) %>% t() -> ELASTIC_NET 
cbind(group_list,ELASTIC_NET ) %>% as.data.frame() -> roc_ELASTIC_NET 
colnames(roc_ELASTIC_NET)
ELASTIC_NET_model <- glm(group_list~ MMP9+ALOX5+MAPK14+C5AR1,roc_ELASTIC_NET,family = "binomial",control=list(maxit=100))

library(performance)
library(ggplot2)
pdf("ROC_组合.pdf",width = 8,height = 8,family="serif")
output <- performance_roc(GBM_model,ridge_model,lasso_model,xboost_model,
                          pls_da_model,SVM_RFE_model,ELASTIC_NET_model,randomforest_model)
#output <- output %>%
    # group_by(Model,Sensitivity) %>%
    # summarise(mean_sp = mean(Specificity))
output$Model <- factor(output$Model,levels=c("GBM_model","ridge_model","lasso_model","xboost_model",
                                             "pls_da_model","SVM_RFE_model","ELASTIC_NET_model","randomforest_model"))

ggplot(output,aes(x=Specificity, y=Sensitivity, fill=Model,color=Model))+
    geom_line(size=1)+
    mytheme+
    scale_y_continuous(limits = c(-0.1,1.2),expand = c(0,0))+
    labs(x="1 - Specificity (False Positive Rate)",y="Sensitivity (True Positive Rate)",title = "")+
    scale_color_discrete(labels = c("GBM model: 95.00%",
                                    "Ridge model: 93.37%",
                                    "Lasso model: 95.83%",
                                    "Xboostmodel: 93.33%",
                                    "PLS Logistic model: 96.67%",
                                    "SVM-RFE_model: 95.83%",
                                    "Elastic Net model: 96.67%",
                                    "Randomforest model: 97.50%"))
    
    

dev.off()


# 3.3 出图
if(T){mytheme <- theme(text = element_text(family  = "serif"),
                       plot.title = element_text(size = 20,color="black",hjust = 0.5,family="serif",face = "bold"),
                       axis.title = element_text(size = 20,color ="black", family="serif",face = "bold"), 
                       axis.text = element_text(size= 20,color = "black", family="serif",face = "bold"),
                       #panel.grid.minor.y = element_blank(),
                       #panel.grid.minor.x = element_blank(),
                       axis.text.x = element_text(hjust = 1,family="serif",face = "bold" ),
                       legend.text = element_text(size= 20, family="serif",face = "bold"),
                       legend.title= element_text(size= 20, family="serif",face = "bold"),
                       axis.text.y = element_text(margin = margin(0,0,0,0.2,'cm')),
                       axis.ticks.length = unit(.40, "cm"),
                       legend.justification=c(2,0),
                       legend.position=c(1.6,0.1),
                       axis.line = element_line(size = 1.5),
                       axis.ticks = element_line(size = 0.8),
                       plot.margin = margin(t = 2),
                       axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0),size= 30, family="serif",face = "bold"),
                       axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size= 30, family="serif",face = "bold"),
                       panel.border = element_rect(linewidth = 2)
                       #legend.spacing.y = unit(-100, 'cm')
                       #axis.ticks = element_line(linewidth = 2,size = 2)
) }


#引用包
library(timeROC)
#导入初始数据：

roc_randomforest -> roc

y=colnames(roc)[1]
#定义颜色
bioCol=c("#35A585","#CC78A6",'#0796E0',"#C71585")
if(ncol(roc)>4){
    bioCol=rainbow(ncol(roc))}

#绘制
pdf("roc曲线.pdf",width=6,height=6,family = "serif")
par(family="serif", font = 2)
roc1=roc(roc[,y], as.vector(roc[,2]), levels = c( 1,0), direction = "auto")
aucText=c( paste0(colnames(roc)[2],", AUC=",sprintf("%0.3f",auc(roc1))))
plot(roc1, col=bioCol[1], xlab="1 - Specificity (False Positive Rate)", ylab="Sensitivity (True Positive Rate)",
     cex.axis = 1.5, cex.lab = 1.6, lwd = 4, font = 2)
for(i in 3:ncol(roc)){
    roc1=roc(roc[,y], as.vector(roc[,i]))
    lines(roc1, col=bioCol[i-1], lwd = 4)
    aucText=c(aucText, paste0(colnames(roc)[i],", AUC=",sprintf("%0.3f",auc(roc1))) )
}

legend("bottomright", aucText,lwd=2,bty="n",col=bioCol[1:(ncol(roc)-1)],cex=1)
dev.off()


