
library(tidyverse)
library(readxl)
library(xlsx)
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
library(ggbeeswarm)
library(cols4all)
library(cowplot)
library(extrafont)
loadfonts(device = "pdf")
#导入字体进入pdf
#showtext包可给定字体文件，加载到 R环境中，生成新的字体家族名字，后期调用这个名字设定字体，并且支持中文写入pdf不乱码
library(showtext)
showtext_auto(enable=TRUE)

#同源基因
library(biomaRt)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
drerio <- useMart("ensembl", dataset = "drerio_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
getLDS(attributes = "hgnc_symbol", # 要转换符号的属性，这里基因名（第3步是基因名）
       filters = "hgnc_symbol", #参数过滤
       mart = human, #需要转换的基因名的种属来源，也就是第2步的mouse
       values = gene_index_randomforest$Feature,
       attributesL = "zfin_id_symbol", #要同源转换的目标属性，这里还是转为基因名，也可加其他
       martL = drerio, #要同源转换的目标种属，也就是第2步的human
       uniqueRows = TRUE#要转换的基因集
) -> z2


rm(list=ls())
gene_count <- read_excel("gene_count.xls")
gene_count[gene_count$gene_name %in% z2$ZFIN.symbol,] %>%
    as.data.frame() %>%
    select(c("Ctl_1","Ctl_2","Ctl_3","M_1", "M_2","M_3","gene_name")) %>%
    set_names(c(rep("Ctl",3),rep("M",3),"gene_name")) %>%
    reshape2::melt() -> final
final$variable %>%  fct_inorder()  -> final$variable
ggplot(final, aes(x = gene_name, y = value,color = variable,fill = variable))+ 
    labs(y="",x= "",title = "")+  
    geom_point() +
    scale_color_manual(values = rev(mycol)) +
    theme_classic()+
    mytheme


if(T){mytheme <- theme(plot.title = element_text(size = 20,color="black",hjust = 0.5,family="Arial",face = "bold"),
                       axis.title = element_text(size = 20,color ="black", family="Arial",face = "bold"), 
                       axis.text = element_text(size= 20,color = "black", family="Arial",face = "bold"),
                       panel.grid.minor.y = element_blank(),
                       panel.grid.minor.x = element_blank(),
                       axis.text.x = element_text(angle = 45,hjust = 1,family="Arial",face = "bold" ),
                       panel.grid=element_blank(),
                       legend.position = "right",
                       legend.text = element_text(size= 20, family="Arial",face = "bold"),
                       legend.title= element_text(size= 20, family="Arial",face = "bold"),
                       axis.text.y = element_text(margin = margin(0,0,0,0.2,'cm')),
                       axis.ticks.length.y = unit(.25, "cm"),
                       axis.ticks = element_line(linewidth = 2,size = 2)
) }

c4a_gui()
mycol <- c4a('bold',8) #自定义配色挑选
gene_count[gene_count$gene_name %in% z2$ZFIN.symbol,] %>%
    as.data.frame() %>%
    select(c("Ctl_1","Ctl_2","Ctl_3","M_1", "M_2","M_3","gene_name")) %>%
    set_names(c(rep("Ctl",3),rep("M",3),"gene_name"))  -> final


gene_count[gene_count$gene_name %in% z2$ZFIN.symbol,] %>%
    as.data.frame() %>%
    select(c("Ctl_1","Ctl_2","Ctl_3","M_1", "M_2","M_3","gene_name")) %>%
    set_names(c(rep("Ctl",3),rep("M",3),"gene_name")) -> final
rownames(final) <- final$gene_name 
final <- final[,-7]
final <- log2(final+1) %>% set_names(c(rep("Ctl",3),rep("M",3)))
rownames(final) -> final$gene
final <- melt(final)
