#安装需要的R包
install.packages("ggplot2")
install.packages("tidyverse")
install.packages("ggrepel")
install.packages("corrplot")
install.pachages("reshape2")
#2.载入需要的R包

library(ggplot2)
library(tidyverse)
library(ggrepel)
library(corrplot)
library(reshape2)
#3.读取数据

cor <- gsva_inflam[,c(3,5,6)] %>%
    set_names(c("Inflammation", "Depression", "T_cells_CD8"))
cor[,2] <- case_when(cor[,2]=="Depression"~1,
                       cor[,2]=="Control"~0)
cor_f <- train %>%
    select(c("VDR","ALOX5","C5AR1","HMGCR")) %>%
    bind_cols(cor)
#目标基因集
genelist<- c("VDR","ALOX5","C5AR1","HMGCR")
#提取基因集的表达矩阵

#计算相关系数
comcor<-cor(cor_f)
#计算显著性差异
comp<-cor.mtest(comcor,conf.level=0.95)
pval<-comp$p
#获取目标基因相关性矩阵
goalcor<-select(as.data.frame(comcor),genelist)%>%rownames_to_column(var="celltype")
goalcor<-filter(goalcor,!(celltype %in% genelist))
##长宽数据转换
goalcor<-melt(goalcor,id.vars="celltype")
colnames(goalcor)<-c("celltype","Gene","correlation")
#获取目标基因集pvalue矩阵
pval<-select(as.data.frame(pval),genelist)%>%rownames_to_column(var="celltype")
pval<-filter(pval,!(celltype %in% genelist))
#长宽数据转换
pval<-melt(pval,id.vars="celltype")
colnames(pval)<-c("celltype","gene","pvalue")
#将pvalue和correlation两个文件合并
final<-left_join(goalcor,pval,by=c("celltype"="celltype","Gene"="gene"))
#5.开始绘图

#添加一列,来判断pvalue值范围
final$sign<-case_when(final$pvalue >0.05~">0.05",
                      final$pvalue <0.05 &final$pvalue>0.01 ~"<0.05",
                      final$pvalue <0.01 &final$pvalue>0.001 ~"<0.01",
                      final$pvalue <0.001 &final$pvalue>0.0001~"<0.001",
                      final$pvalue <0.0001 ~"<0.0001")
#添加一列来判断correlation的正负
final$core<-case_when(final$correlation >0 ~"positive",
                      final$correlation<0 ~"negtive")
#开始绘图
pdf("final_check.pdf", width=8, height=7,family = "serif")
ggplot(data=final,aes(x=Gene,y=celltype))+
    #筛选正相关的点
    geom_point(data=filter(final,core=="positive"),aes(x=Gene,y=celltype,size=abs(correlation),fill=sign),shape=21)+
    #筛选负相关的点
    geom_point(data=filter(final,core=="negtive"),
               aes(x=Gene,y=celltype,size=abs(correlation),color=sign),shape=16)+
    #自定义颜色
    scale_color_manual(values=c("#333366","#336699","#3399CC","#66CCFF","#CCCCCC"))+
    #自定义填充颜色
    scale_fill_manual(values=c("#FF6666","#FF9999","#FFCCCC","#FFCCCC","#CCCCCC"))+
    #去除x轴和y轴
    labs(x="",y="")+
    #修改主题
    theme_bw()+
    theme(axis.text.x=element_text(angle=45,hjust=1),
          axis.ticks.x=element_blank())+
    #修改legend,title表示标题，order表示标题的顺序，override.aes来设置大小
    guides(color=guide_legend(title="Negitive\ncorrelation\nFDR q-value",order=1,override.aes=list(size=4)),
           size=guide_legend(title="Spearman's p",order=2),
           fill=guide_legend(title="Positive\ncorrelation\nFDR q-value",order=3,override.aes=list(size=4)))+
    mytheme
#保存图片
dev.off()
if(T){mytheme <- theme(text = element_text(family = "serif"),
                     plot.title = element_text(size = 20,color="black",hjust = 0.5),
                     axis.title = element_text(size = 20,color ="black"), 
                     axis.text = element_text(size= 20,color = "black"),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1 ),
                     panel.grid=element_blank(),
                     legend.position = "right",
                     legend.text = element_text(size= 20,color = "black"),
                     legend.title= element_text(size= 20,color = "black"),
                     axis.line = element_line(size = 0.3),
                     axis.ticks = element_line(size = 0.5),
                     panel.border = element_rect(linewidth = 2.5)
) }
