R.utils::setOption("clusterProfiler.download.method","auto")
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(pathview)
library(ggnewscale)
library(DOSE)
library(stringr)
library(R.utils)
library(extrafont)
loadfonts(device = "pdf")
#导入字体进入pdf
#showtext包可给定字体文件，加载到 R环境中，生成新的字体家族名字，后期调用这个名字设定字体，并且支持中文写入pdf不乱码
library(showtext)
showtext_auto(enable=TRUE)

pvalueFilter=0.05        
qvalueFilter=1        
showNum=10
keggId="hsa04659"

rt <- final_genes %>%
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
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]

write.table(KEGG,file="KEGG.xls",sep="\t",quote=F,row.names = F)

if(nrow(KEGG)<showNum){
	showNum=nrow(KEGG)
}

pdf(file="KEGG_bubble.pdf",width = 17,height = 15,family="Arial")
dotplot(kk, showCategory = showNum, orderBy = "GeneRatio",color = colorSel)+scale_y_discrete(labels=function(x) stringr::str_wrap(x, width=60))+
    theme(text=element_text(family="Arial"),
          axis.text.y = element_text(size= 25,color = "black",face= "bold"),
          axis.text.x = element_text(size= 25,color = "black",face= "bold"),
          legend.text = element_text(size= 25,face= "bold"),
          legend.title= element_text(size= 25,face= "bold"),
          #axis.text.x.top = element_text(size= 30,face= "bold"),
          #axis.text.y.left = element_text(size= 30,face= "bold"),
          #strip.text.x = element_text(size = 30, face = "bold"),
          strip.text.y = element_text(size = 25, face = "bold"),#小标题大 
          legend.key.size = unit(1.5, "cm"),
          axis.title.x.bottom =  element_text(size = 25, face = "bold"))#图例大小
dev.off()

source("cnetplot.R")
pdf(file="KEGG_cnet.pdf",width = 30,height = 20,family="serif")
#par(family = "serif")
#low='darkgreen', high='firebrick'
af=setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
p <- cnetplot(af, showCategory = 10,categorySize="pvalue", family = "serif", 
         circular = TRUE,
         color.params = list(edge = TRUE, category = "darkgreen", gene ="firebrick"),
         cex.params = list(category_node = 2, gene_node = 2, 
                           category_label = 3, gene_label = 3, family = "serif"
                           ),node_label = "all")
p$layers[[3]]$aes_params$family <- "serif"  
p+theme(
        text = element_text(face = "bold", family = "serif"),
        #axis.text.x.top = element_text(size = 30, face = "bold", family = "serif"),
        #axis.text.y.left = element_text(size = 30, face = "bold", family = "serif"),
        #strip.text.x = element_text(size = 30, face = "bold"),
        #strip.text.y = element_text(size = 25),  # 修改小标题大
        legend.text = element_text(size = 30, family = "serif"),
        legend.title = element_text(size = 40),
        legend.key.size = unit(2.5, "cm")
    )
dev.off()
