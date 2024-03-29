
library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
lapply(c('clusterProfiler','enrichplot','patchwork'), function(x) {library(x, character.only = T)})
library(org.Hs.eg.db)
library(patchwork)
library(ggridges)
library(extrafont)
library(tidyverse)
loadfonts(device = "pdf")
#导入字体进入pdf
#showtext包可给定字体文件，加载到 R环境中，生成新的字体家族名字，后期调用这个名字设定字体，并且支持中文写入pdf不乱码
library(showtext)
showtext_auto(enable=TRUE)
#输入文件
expFile="ids_exprs"         
sgene="VDR"       #输入进行单基因GSEA的基因名称

rt=ids_exprs <- read.delim("C:/Users/79403/Desktop/传森/ids_exprs.txt", row.names=1)
rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
#按基因表达分组
group <- ifelse(data[c(sgene),]>median(data[c(sgene),]), "High", "Low")   
group <- factor(group,levels = c("High","Low"))

#差异分析
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(High-Low,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
deg=topTable(fit2,adjust='fdr',number=nrow(data))
Diff=deg
#保存单基因分组的所有基因差异结果
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file=paste0("1.","DIFF_all.xls"),sep="\t",quote=F,col.names=F)


#展示差异最大的前30个基因
Diff=Diff[order(as.numeric(as.vector(Diff$logFC))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(60)){
  afGene=diffGene[c(1:30,(diffLength-30+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=data[afGene,]
#分组标签
Type1=as.data.frame(group)
Type1=Type1[order(Type1$group,decreasing = T),,drop=F]
Type=Type1[,1]
names(Type)=rownames(Type1)
Type=as.data.frame(Type)
#分组标签的注释颜色
anncolor=list(Type=c(High="red",Low="blue"  ))

pdf(file=paste0("1.", "DIFF_heatmap.pdf"),height=7,width=6,family = "Arial")
pheatmap(afExp[,rownames(Type1)],                                                                      #热图数据
         annotation=Type,                                                            #分组
         color = colorRampPalette(c("blue","white","red"))(50),     #热图颜色
         cluster_cols =F,                                                           #不添加列聚类树
         show_colnames = F,                                                         #展示列名
         scale="row", 
         fontsize = 10,
         fontsize_row=6,
         fontsize_col=8,
         annotation_colors=anncolor
)
dev.off()

#火山图差异标准设置
adjP=0.05
aflogFC=0.5
Significant=ifelse((Diff$P.Value<adjP & abs(Diff$logFC)>aflogFC), ifelse(Diff$logFC>aflogFC,"Up","Down"), "Not")
#开始绘制
p = ggplot(Diff, aes(logFC, -log10(P.Value)))+
  geom_point(aes(col=Significant),alpha=0.4, size=2)+               #size点的大小
  scale_color_manual(values=c(pal_npg()(2)[2], "#838B8B", pal_npg()(1)))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  geom_hline(aes(yintercept=-log10(adjP)), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=aflogFC), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=-aflogFC), colour="gray", linetype="twodash",size=1)+
    scale_colour_manual(values = c('steelblue','gray','#F2C661'))
#查看，不添加标记可以直接保存
p
#添加基因点标记，按照,可自行根据差异分析的结果进行标记
point.Pvalue=0.01
point.logFc=1.5
#继续绘制
Diff$symbol=rownames(Diff)
pdf(paste0("1.", "DIFF_vol.pdf"),width=7,height=6,family = "Arial")
p=p+theme_bw()
for_label <- Diff %>% 
  filter(abs(logFC) >point.logFc & P.Value< point.Pvalue )
p+geom_point(size = 1.5, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black",
    label.size =0.1
  )
dev.off()


logFC_t=0
deg$g=ifelse(deg$P.Value>0.05,'stable',
             ifelse( deg$logFC > logFC_t,'UP',
                     ifelse( deg$logFC < -logFC_t,'DOWN','stable') )
)
table(deg$g)

deg$symbol=rownames(deg)
df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
           toType = c( "ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=deg
DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
data_all_sort <- DEG %>% 
  arrange(desc(logFC))

geneList = data_all_sort$logFC #把foldchange按照从大到小提取出来
names(geneList) <- data_all_sort$ENTREZID #给上面提取的foldchange加上对应上ENTREZID
head(geneList)

#开始富集分析
kk2 <- gseKEGG(geneList     = geneList,
               organism     = "hsa",
               nPerm        = 10000,
               minGSSize    = 10,
               maxGSSize    = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none" )
class(kk2)
colnames(kk2@result)
kegg_result <- as.data.frame(kk2)
rownames(kk2@result)[head(order(kk2@result$enrichmentScore))]
af=as.data.frame(kk2@result)
write.table(af,file=paste0("2.","all_GSEA.xls"),sep="\t",quote=F,col.names=T)

#排序后分别取GSEA结果的前5个和后5个
num=5
pdf(paste0("2.","down_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[head(order(kk2@result$enrichmentScore),num)])
dev.off()
pdf(paste0("2.","up_GSEA.pdf"),width = 8,height = 8)
gseaplot2(kk2, geneSetID = rownames(kk2@result)[tail(order(kk2@result$enrichmentScore),num)])
dev.off()
#排序后取前5个和后5个一起展示
num=5
pdf(paste0("2.","all_GSEA.pdf"),width = 10,height = 20,family = "Arial")
gseaplot2(kk2, geneSetID = rownames(kk2@result)[c(head(order(kk2@result$enrichmentScore),num),tail(order(kk2@result$enrichmentScore),num))])
dev.off()

#单独展示,自行保存

#山脊图，展示10个，自行保存
library(stringr)
kk2@result$Description=gsub("HALLMARK_","",kk2@result$Description)
pdf(paste0("2.","ridgeplot.pdf"),width = 10,height = 15,family = "Arial")
ridgeplot(kk2,showCategory = 20)
dev.off()

