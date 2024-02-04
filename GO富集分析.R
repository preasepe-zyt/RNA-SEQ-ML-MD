final_genes <- read.csv("C:/Users/79403/Desktop/时瑞蝶/交集和kegg-go/final_genes.txt", row.names=1)
library(org.Hs.eg.db)  
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(enrichplot)
library(DOSE)
library(stringr)
library(extrafont)
loadfonts(device = "pdf")
#导入字体进入pdf
#showtext包可给定字体文件，加载到 R环境中，生成新的字体家族名字，后期调用这个名字设定字体，并且支持中文写入pdf不乱码
library(showtext)
showtext_auto(enable=TRUE)

pvalueFilter=0.05         
qvalueFilter=1  
showNum=8

rt <- final_genes  %>%
    as.data.frame()
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)

colorSel="qvalue"
if(qvalueFilter>0.05){
  colorSel="pvalue"
}


kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]

write.table(GO,file="GO.xls",sep="\t",quote=F,row.names = F)


if(nrow(GO)<30){
  showNum=nrow(GO)
}

pdf(file="GO_bubble.pdf",width = 20,height =22,family="serif")
bub=dotplot(kk,showCategory = 5, orderBy = "GeneRatio",split="ONTOLOGY", color = colorSel) + 
    facet_grid(ONTOLOGY~., scale='free')+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))+
    theme(text=element_text(face= "bold", family = "serif"),
          axis.text.y = element_text(size= 30,color = "black",face= "bold"),
          axis.text.x = element_text(size= 30,color = "black",face= "bold"),
          legend.text = element_text(size= 30,face= "bold"),
          legend.title= element_text(size= 30,face= "bold"),
          #axis.text.x.top = element_text(size= 30,face= "bold"),
          #axis.text.y.left = element_text(size= 30,face= "bold"),
          #strip.text.x = element_text(size = 30, face = "bold"),
          strip.text.y = element_text(size = 30, face = "bold"),#小标题大 
          legend.key.size = unit(2.5, "cm"),
          axis.title.x.bottom =  element_text(size = 30, face = "bold"),
          axis.line = element_line(size = 1.5),
          panel.border = element_rect(linewidth = 3))#图例大小
print(bub)
dev.off()

pdf(file="GO_cnet.pdf",width = 17,height = 10,family="Arial")
af=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(kk, showCategory = 5, categorySize="pvalue",colorEdge = TRUE,cex_label_category = 1.5,cex_label_gene = 1.5,circular = TRUE)+
    theme(text=element_text(family="Arial"),
          #axis.text.y = element_text(size= 10,color = "black",face= "bold"),
          #axis.text.x = element_text(size= 10,color = "black",face= "bold"),
          legend.text = element_text(size= 17,face= "bold"),
          legend.title= element_text(size= 17,face= "bold"),
          #axis.text.x.top = element_text(size= 30,face= "bold"),
          #axis.text.y.left = element_text(size= 30,face= "bold"),
          #strip.text.x = element_text(size = 30, face = "bold"),
          strip.text.y = element_text(size = 17, face = "bold"),#小标题大 
          legend.key.size = unit(2.3, "cm"))
          #axis.title.x.bottom =  element_text(size = 10, face = "bold"))#图例大小
dev.off()
