rm(list=ls())  #清空环境内变量
options(stringsAsFactors = F)  #避免自动将字符串转换为R语言因子
options(digits = 2)#保留两位小数
suppressMessages()#抑制一些信息
library(GEOquery)
library(stringr)
library(ggplot2)
library(reshape2)
library(limma)
library(extrafont)
library(tidyverse)
loadfonts(device = "pdf")
#导入字体进入pdf
#showtext包可给定字体文件，加载到 R环境中，生成新的字体家族名字，后期调用这个名字设定字体，并且支持中文写入pdf不乱码

library(showtext)
showtext_auto(enable=TRUE)
#表达矩阵
gset <- getGEO("GSE98793", destdir=".", AnnotGPL = F, getGPL = F)  #设置后两个为T比较耗时，而且作用不大
exp<-exprs(gset[[1]])  #exp即为表达矩阵
pdata<- pData(gset[[1]]) #使用pData()函数提取临床信息
GPL <- gset[[1]]@annotation 
#pdata %>% filter(`treatment:ch1` == "none") -> pdata
GSE98793 <- exp3[which(rownames(exp3)=="ALOX5"),]



#临床信息中哪一列提供了分组信息需要自己去鉴别
#采用字符串分隔函数按空格进行分隔取第三位
group_list<- pdata$`subject group:ch1`
group_list<- ifelse(str_detect(pdata$`subject group:ch1`,"CNTL; healthy control"),"Control","Major Depressive Disorder")
#设置样本分组
#group_list <- factor(group_list,levels = c("Control","Major Depressive Disorder"))  #设置为因子变量
table(group_list)
GSE98793 <- rbind(group=group_list,value=GSE98793,dataset=rep("GSE98793",nrow(GSE98793))) %>% 
            t() %>%
            as.data.frame()
GSE98793$value <- as.numeric(GSE98793$value)


#GSE39653
gset <- getGEO("GSE39653", destdir=".", AnnotGPL = F, getGPL = F)  #设置后两个为T比较耗时，而且作用不大
exp<-exprs(gset[[1]])  #exp即为表达矩阵
pdata<- pData(gset[[1]]) #使用pData()函数提取临床信息
GPL <- gset[[1]]@annotation 
#pdata %>% filter(`treatment:ch1` == "none") -> pdata
GSE98793 <- exp3[which(rownames(exp3)=="ALOX5"),]



#临床信息中哪一列提供了分组信息需要自己去鉴别
#采用字符串分隔函数按空格进行分隔取第三位
group_list<- pdata$`subject group:ch1`
group_list<- ifelse(str_detect(pdata$`subject group:ch1`,"CNTL; healthy control"),"Control","Major Depressive Disorder")
#设置样本分组
#group_list <- factor(group_list,levels = c("Control","Major Depressive Disorder"))  #设置为因子变量
table(group_list)
GSE98793 <- rbind(group=group_list,value=GSE98793,dataset=rep("GSE98793",nrow(GSE98793))) %>% 
    t() %>%
    as.data.frame()
GSE98793$value <- as.numeric(GSE98793$value)

#箱型图
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)

if(T){mytheme <- theme(text = element_text(family = "serif"),
                       plot.title = element_text(size = 20,color="black",hjust = 0.5,family="serif"),
                       axis.title = element_text(size = 20,color ="black", family="serif"), 
                       axis.text = element_text(size= 20,color = "black", family="serif"),
                       #panel.grid.minor.y = element_blank(),
                       #panel.grid.minor.x = element_blank(),
                       axis.text.x = element_text(angle = 45,hjust = 1,family="serif"),
                       #panel.grid=element_blank(),
                       legend.position = "top",
                       axis.line = element_line(size = 0.8),
                       axis.ticks = element_line(size = 0.8),
                       legend.text = element_text(size= 20, family="serif"),
                       legend.title= element_text(size= 20, family="serif"),
                       axis.text.y = element_text(margin = margin(0,0,0,0.2,'cm')),
                       axis.ticks.length.y = unit(.25, "cm"),legend.direction = "vertical"
                       #axis.ticks = element_line(linewidth = 2,size = 2)
) }
pdf("外部验证.pdf", width=6.2, height=6,family = "Arial")
ggplot(GSE98793, aes(x = group, y = value))+
    geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
    scale_fill_manual(values = c("#BEBADA", "#FDB462"))+
    theme_classic() + mytheme + 
    stat_compare_means(aes(group =  group),
                       label = "p.signif",
                       method = "anova",
                       hide.ns = T,
                       label.y = 7.5,
                       label.x = 1.45,
                       size=10)+
    geom_signif(
        comparisons = list(c("Control", "Major Depressive Disorder")),  # 替换为你实际的组合
        y_position = 7.5,
        textsize = 0,
        size =  1# 调整横线的位置
        
    )+#"p.signif""p.format"
    scale_y_continuous(limits = c(0,9),expand = c(0,0),breaks = c(1:9))+
    facet_grid(. ~ dataset)+
    labs(x="")
dev.off()    
