library(ggvenn)
library(readxl)
library(extrafont)
loadfonts(device = "pdf")
#导入字体进入pdf
#showtext包可给定字体文件，加载到 R环境中，生成新的字体家族名字，后期调用这个名字设定字体，并且支持中文写入pdf不乱码
library(showtext)
showtext_auto(enable=TRUE)
library(tidyverse)
THBP1vsctl_deg_all <- read.delim("C:/Users/79403/Desktop/排序+取交集/THBP1vsctl_deg_all.xls")
THBP4vsctl_deg_all <- read.delim("C:/Users/79403/Desktop/排序+取交集/THBP4vsctl_deg_all.xls")
THBP1vsctl_deg_all %>% subset(gene_name != "-") -> THBP1vsctl_deg_all
THBP4vsctl_deg_all %>% subset(gene_name != "-") -> THBP4vsctl_deg_all
dat <- list( THBP1_vs_ctl = THBP1vsctl_deg_all$gene_name,  THBP4_vs_ctl = THBP4vsctl_deg_all$gene_name)
df1<-intersect(THBP1vsctl_deg_all$gene_name,THBP4vsctl_deg_all$gene_name) #取三个交集
pdf("共同变化基因.pdf",width = 8,height = 8,family = "Arial")
ggvenn(dat,show_percentage = T,
       stroke_color = "white",
       stroke_size = 0.5,
       fill_color = c("#CC78A6",'#0796E0'),
       set_name_color =c("#CC78A6",'#0796E0'), 
       set_name_size = 5,
       text_size=4)
dev.off()

#高浓度筛选
THBP4vsctl_deg_all %>% subset(gene_name != "-") %>% select("gene_name","log2FoldChange") -> all_dif
all_dif[all_dif$gene_name %in% df1,] %>% arrange(desc(abs(log2FoldChange))) -> ord
xlsx::write.xlsx(ord,"排序.xlsx")
xlsx::write.xlsx(df1,"共同变化基因.xlsx")