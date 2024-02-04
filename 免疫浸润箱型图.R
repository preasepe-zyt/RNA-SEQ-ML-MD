模块一. Cibersort算法的纯代码实现及结果绘图展示
remove(list = ls()) #一键清空
#加载包
library(ggplot2)
library(reshape2)
library(ggpubr)
library(dplyr)
'
1. Cibersort计算免疫细胞
加载cibersort的官方提供的源码，指定基准数据库文件 (LM22.txt，这是22种免疫细胞的marker基因，下载自Cibersort官网)。
'''
source('./assist/Cibersort.R')

# 设置分析依赖的基础表达文件
# 每类免疫细胞的标志性基因及其表达
# 基因名字为Gene symbol
LM22.file <- "./database/LM22.txt"
加载自己的数据用于分析计算免疫细胞

# 1. Cibersort

TCGA_exp.file <- "./Rawdata/TCGA_HNSC_mRNA_fpkm_paired_43vs43.txt"

TCGA_TME.cibersort <- CIBERSORT(LM22.file ,TCGA_exp.file, perm = 50, QN = F)  
# perm置换次数=1000
# QN如果是芯片设置为T，如果是测序就设置为F

write.csv(TCGA_TME.cibersort, "./Output/TCGA_CIBERSORT_cibersort_fromRcode.csv")

2. 分组信息
## 2. 分组信息

# TCGA的数据还可以从名字获取
#group_list <- ifelse(as.numeric(substring(rownames(TCGA_TME.cibersort),14,15)) < 10,
#                 "Tumor","Normal") %>% 
# factor(.,levels = c("Normal","Tumor"))

phenotype = read.csv("./Rawdata/TCGA_HNSC_paired_metadata.csv",header = T,row.names = 1)

group_list <- phenotype$group %>% 
    factor(.,levels = c("Nontumor","Tumor"))

'
3. 绘图
3.1 数据转换预处理，取前22列，忽略掉后面计算出的P-value,Correlation, RMSE单列信息。
'''
## 3. 绘图
# 3.1 数据粗处理
TME_data <- as.data.frame(cibersort[,2:ncol(LM22)])

TME_data$group <- group_list
TME_data$sample <- ID
TME_data <- TME_data %>% select(c("Macrophages_M0_CIBERSORT","Macrophages_M0_CIBERSORT","NK_cells_resting_CIBERSORT","T_cells_CD4_memory_resting_CIBERSORT","T_cells_CD4_memory_activated_CIBERSORT","T_cells_CD4_naive_CIBERSORT","T_cells_CD8_CIBERSORT","group"))
# 2.2 融合数据
TME_New = melt(TME_data)

## Using group, sample as id variables

colnames(TME_New)=c("Group","Celltype","Composition")  #设置行名
head(TME_New)

#3.2 按免疫细胞占比中位数排序绘图（可选）
# 3.3 按免疫细胞占比中位数排序绘图（可选）
plot_order = TME_New[TME_New$Group=="Healthy",] %>% 
    group_by(Celltype) %>% 
    summarise(m = median(Composition)) %>% 
    arrange(desc(m)) %>% 
    pull(Celltype)

## `summarise()` ungrouping output (override with `.groups` argument)

TME_New$Celltype = factor(TME_New$Celltype,levels = plot_order)
#3.3 绘制箱线图

# 3.3 出图
if(T){mytheme <- theme(text = element_text(family = "serif")
                       ,plot.title = element_text(size = 30,color="black",hjust = 0.5),
                       axis.title = element_text(size = 30,color ="black"), 
                       axis.text = element_text(size= 30,color = "black"),
                       #panel.grid.minor.y = element_blank(),
                       #panel.grid.minor.x = element_blank(),
                       axis.text.x = element_text(angle = 45,hjust = 1,color = "black" ),
                       #panel.grid=element_blank(),
                       legend.position = "top",
                       legend.text = element_text(size= 20),
                       legend.title= element_text(size= 20),
                       axis.text.y = element_text(margin = margin(0,0,0,0.2,'cm')),
                       axis.ticks.length.y = unit(0.25, "cm"),
                       axis.title.y = element_text(size = 30)
                       #axis.ticks = element_line(linewidth = 2,size = 2)
) }

gsub("Healthy","Control",TME_New$Group) -> TME_New$Group
gsub("_CIBERSORT","",TME_New$Celltype) -> TME_New$Celltype

pdf(file="immune_boxplot.pdf",height=8,width=10,family = "serif")
ggplot(TME_New, aes(x = Celltype, y = Composition))+ 
    labs(y="Cell composition",x= NULL,title = "Immune Cell composition")+  
    geom_boxplot(aes(fill = Group),position=position_dodge(0.5),width=0.5,outlier.alpha = 0)+ 
    scale_fill_manual(values = c("#3690C0", "#11A579"))+
    theme_classic() + 
    mytheme + 
    stat_compare_means(
        aes(group = Group),
        label = "p.signif",
        method = "wilcox.test",
        hide.ns = TRUE,
        size = 10,
        label.y = 0.15) +
    geom_signif(
        comparisons = list(c("Control", "Treatment")),  # 你要比较的组合
        y_position = 4,  # 调整标记线的纵坐标位置
        textsize = 5     )+ # 调整标记文本的大小
        scale_y_continuous(limits = c(0,0.3),expand = c(0,0),breaks = c(0,0.1,0.2,0.3))

dev.off()

write.table(TME_New,file = "immune_box.txt",sep = "\t",row.names=T,col.names = T)

