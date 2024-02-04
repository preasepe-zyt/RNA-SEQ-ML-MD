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
gset <- getGEO("GSE76826", destdir=".", AnnotGPL = F, getGPL = F)  #设置后两个为T比较耗时，而且作用不大
exp<-exprs(gset[[1]])  #exp即为表达矩阵
pdata<- pData(gset[[1]]) #使用pData()函数提取临床信息
pdata[c(1:22),] -> pdata
exp[,c(1:22)] -> exp
GPL <- gset[[1]]@annotation 
#临床信息中哪一列提供了分组信息需要自己去鉴别
#采用字符串分隔函数按空格进行分隔取第三位
group_list<- str_split(pdata$title,'_',simplify = T)[1]
group_list <- pdata$`state:ch1`
group_list<- ifelse(str_detect(pdata$`state:ch1`,"Healthy"),"Control","Depression")
#设置样本分组
group_list <- factor(group_list,levels = c("Control","Depression"))  #设置为因子变量
table(group_list)

#还可以通过PCA来查看两组数据的分布差异:
pdf("PCA.pdf", width=6.2, height=5.5,family = "Arial")
library(FactoMineR)
library(factoextra)
dat<- as.data.frame(t(exp3))
dat.pca <- PCA(dat, graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = group_list,
             addEllipses = TRUE,
             palette = "jco", 
             legend.title = "Groups"
)
dev.off()

#check
expFile="ids_exprs.txt"      #数据读取
#读取表达矩阵数据
data=read.table(expFile, header=T, row.names = 1,sep="\t", check.names=F)
##

#data <- normalizeBetweenArrays(exp3, method="quantile")#批次校正
#data <- normalizeBetweenArrays(exp3, method="scale")
##看看数据分布
n.sample=ncol(data)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
pdf("批次效应.pdf", width=15, height=7,family = "Arial")
par(family="Arial",mfrow = c(2, 2),cex = 1)

boxplot(data, col = cols,main="expression value",las=2) 
#看一下
#还差点意思
library(limma)
##看看数据分布
data2 <- normalizeBetweenArrays(data, method="quantile")#批次校正
boxplot(data2, col = cols,main="expression value",las=2)##箱线图看一下 
dev.off()

#差异分析
design=model.matrix(~group_list)
fit=lmFit(data,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)  #提取所有基因的差异分析结果
colnames(deg)
##设置logFC值和P值来标记上下调基因
logFC=0.5
P.Value = 0.05
type1 = (deg$adj.P.Val < P.Value)&(deg$logFC < -logFC)
type2 = (deg$adj.P.Val< P.Value)&(deg$logFC > logFC)
deg$Group = ifelse(type1,"Down",ifelse(type2,"Up","Not-Sig"))
table(deg$Group)
write.table(deg,file = "diff.txt",sep = "\t",row.names=T,col.names = T)

##展示特定的gene在图中
library(ggplot2)
pdf("vol.pdf", width=8, height=8,family = "serif")
p <- ggplot(deg, aes(logFC, -log10(adj.P.Val),color=Group)) + 
  geom_point(alpha=0.8, size=2) +  # 设置点的透明度和大小
  theme_bw(base_size = 12) +  #设置一个主题背景
  xlab("Log2 (Fold change)") + # x轴名字
  ylab("-Log10 (adj.P.Val)") + # y轴名字
  scale_colour_manual(values = c("#BEBADA",'gray',"#FB8072"),
                      labels=c("Down: 603",
                               "Not-Sig: 28073",
                               "Up: 823")) + # 各自的颜色
  geom_hline(yintercept = -log10(0.05), lty = 4) + #定义p值和线形
  geom_vline(xintercept = c(-0.5, 0.5), lty = 4)+ #定义差异倍数和线形
  labs(title = "Control vs Major Depressive Disorder")#加上题目


#加基因标记
point.Pvalue=0.05
point.logFc=0.5
gene <- deg[deg$Group!="Not-Sig",]
gene <- rownames(gene)
for_label <- deg %>% filter(abs(logFC) >point.logFc & P.Value< point.Pvalue )
for_label$sig <- rownames(for_label)
for_label <- for_label[for_label$sig %in% final_genes$x,]
p+geom_point(for_label,mapping=aes(color=Group))+
    ggrepel::geom_label_repel(
    aes(label = sig),
    data = for_label[-5,],
    color="black",
    label.size=1,
    size=5,
    segment.size=0.5, nudge_x=3, direction="y", hjust=0
  )+geom_point(for_label,mapping=aes(color=Group))+
    ggrepel::geom_label_repel(
        aes(label = sig),
        data = for_label[5,],
        color="black",
        label.size=1,
        size=5,
        segment.size=0.5, nudge_x=-3, direction="y", hjust=0
    )+mytheme
    
  #geom_label_repel(aes(label=sign), fontface="bold", color="grey50", box.padding=unit(0.35, "lines"), point.padding=unit(0.5, "lines"), segment.colour = "grey50")
dev.off()
library(pheatmap)
library(ggsci)
#热图展示差异最大的前50个基因
deg=deg[order(as.numeric(as.vector(deg$logFC))),]
degGene=as.vector(rownames(deg))
degLength=length(degGene)
afGene=c()
if(degLength>(100)){
    afGene=degGene[c(1:25,(degLength-25+1):degLength)]
}else{
    afGene=degGene
}
afExp=exp3[afGene,]
#分组标签
Type=c(rep("Control",12),rep("Depression",10))
names(Type)=colnames(data)
Type=as.data.frame(Type)
#分组标签的注释颜色
anncolor=list(Type=c(Control=pal_npg()(1),Depression=pal_npg()(2)[2]))

pdf(file="GSE76826_heatmap.pdf",height=8,width=8)
pheatmap(afExp,                                                                      #热图数据
         annotation=Type,                                                            #分组
         color = c("#7FC97F","white", "#BEAED4"),     #热图颜色
         cluster_cols =F,                                                           #不添加列聚类树
         show_colnames = F,                                                         #展示列名
         scale="row", 
         fontsize = 10,
         #fontsize_row=6,
         #fontsize_col=8,
         annotation_colors=anncolor,
         border_color = "NA",
         labels_row = make_bold_names(afExp, rownames,rownames(afExp))
         )

dev.off()
#箱型图
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
t(data) %>% bind_cols(group_list) %>% melt() -> box
names(box) <- c("group","gene","value")
if(T){mytheme <- theme(
                     text = element_text(family = "serif"),
                     plot.title = element_text(size = 30,color="black",hjust = 0.5),
                     axis.title = element_text(size = 30,color ="black"), 
                     #axis.text = element_text(size= 20,color = "black"),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     axis.text.x = element_text(angle = 0, hjust = 1,size=30,color="black"),
                     axis.text.y = element_text(angle = 0, hjust = 1,size=30,color="black"),
                     panel.grid=element_blank(),
                     legend.position = "top",
                     legend.text = element_text(size= 20),
                     legend.title= element_text(size= 20),
                     axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)),
                     axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0))
    ) }

box_TME <- ggplot(box, aes(x = group, y = value))+ 
    labs(y="Cell composition",x= NULL,title = "TME Cell composition")+  
    geom_boxplot(aes(fill = group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
    scale_fill_manual(values = c("#1CB4B8", "#EB7369"))+
    theme_classic() + mytheme + 
    stat_compare_means(aes(group =  group),
                       label = "p.signif",
                       method = "t.test",
                       hide.ns = T)+
    ylim(0,-2)
