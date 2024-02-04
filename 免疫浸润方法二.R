
# install packages 这三个安装不成功的话，就安后面的bseqsc包也行
install.packages('e1071')
install.pacakges('parallel')
install.packages('preprocessCore')
library(e1071)
library(preprocessCore)
library(parallel)

install.packages('devtools')
library(devtools)
devtools::install_github('shenorrlab/bseqsc')
library(bseqsc)#这个包携带大量CIBERSORT的依赖，前三个安装不好可以安装他

################安装CIBERSORT包##########################################################
if(!require(CIBERSORT))devtools::install_github("Moonerss/CIBERSORT")
library(CIBERSORT)
# 包全部安装完成

# 画热图的包
install.packages("pheatmap")
install.packages("ComplexHeatmap")
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)


# 画热图的包
install.packages("pheatmap")
install.packages("ComplexHeatmap")
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(limma)
# 同时准备好LM22的TXT文件，注意自己以后的文件要和这个TXT的格式一样
# 加载CIBERSORT包成功后，系统内部会自带data(LM22)
data(LM22) 
data(mixed_expr)#TCGA的演示数据，正式情况下就用自己的数据

# 正式开始探索
# 看5*5的数据
LM22[1:5,1:5]
mixed_expr[1:5,1:5]
ids_exprs <- normalizeBetweenArrays(ids_exprs, method="quantile")
ids_exprs <- normalizeBetweenArrays(ids_exprs, method="scale")
# 分别定义signature矩阵LM22和我的数据（演示）矩阵mixed_expr
cibersort <- cibersort(sig_matrix = LM22, mixture_file = ids_exprs)
cibersort$ID -> ID
# 理解一下cibersort的结果
# 你可以理解为返回了一个列名为细胞类型、行名为样本名的细胞浸润程度（占比）的矩阵
# 此外result中还会多出三列：
# P-value: 用来展示去卷积的结果在所有细胞类群中是否具有差异
# Correlation:参考矩阵与输入矩阵的特征基因相关性
# RMSE: Root mean squared error，参考矩阵与输入矩阵的特征基因标准差
library(pheatmap)
library(ggsci)
#分组标签
Type=c(rep("Healthy",12),rep("Depression",10))
names(Type)= ID
Type=as.data.frame(Type)
anncolor=list(Type=c(Healthy=pal_npg()(1),Depression=pal_npg()(2)[2]))
#分组标签的注释颜色
columnscale <- cibersort[,2:ncol(LM22)]
as.vector(rownames(columnscale)) -> type_order
columnscale <- columnscale[type_order,]
columnscale <- columnscale[,apply(columnscale, 2, function(x){sum(x)>0})]#删除全是0的列
rownames(columnscale) <- ID

#字体
library(extrafont)
loadfonts(device = "pdf")
#导入字体进入pdf
#showtext包可给定字体文件，加载到 R环境中，生成新的字体家族名字，后期调用这个名字设定字体，并且支持中文写入pdf不乱码
library(showtext)
showtext_auto(enable=TRUE)

pdf(file="immune_heatmap.pdf",height=8,width=8,family = "Arial")
pheatmap(columnscale,
         scale = 'row',
         cluster_col=F,
         cluster_row=F,
         angle_col = "315",
         annotation_colors=anncolor,
         annotation_row=Type  
         )
dev.off()
# 堆积比例图
my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658','#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398','#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5',  '#CCC9E6', '#625D9E', '#68A180', '#3A6963','#968175'
)
cellnum <- cibersort[,2:ncol(LM22)]
cell.prop<- apply(cellnum, 1, function(x){x/sum(x)})
colnames(cell.prop) <- ID
data4plot <- data.frame()
for (i in 1:ncol(cell.prop)) {
    data4plot <- rbind(
        data4plot,
        cbind(cell.prop[,i],rownames(cell.prop),
              rep(colnames(cell.prop)[i],nrow(cell.prop)
              )
        )
    )
}
colnames(data4plot)<-c('proportion','celltype','sample')
data4plot$proportion <- as.numeric(data4plot$proportion)
gsub("_CIBERSORT","",data4plot$celltype) -> data4plot$celltype
pdf(file="堆积图2.pdf",height=12,width=10,family = "serif")
ggplot(data4plot,aes(sample,proportion,fill=celltype))+
    geom_bar(stat="identity",position="fill")+
    scale_fill_manual(values=my36colors)+#自定义fill的颜色
    ggtitle("Immune Cell Proportion")+
    labs(x="",y="Cell Proportion",fill="Cell Type")+
    theme_bw()+mytheme

dev.off()

write.table(data4plot,file = "data4plot.txt",sep = "\t",row.names=T,col.names = T)

if(T){mytheme <- theme(text = element_text(family = "serif"),
                       plot.title = element_text(size = 40,color="black",hjust = 0.5),
                       axis.title = element_text(size = 20,color ="black"), 
                       axis.text = element_text(size= 20,color = "black"),
                       #panel.grid.minor.y = element_blank(),
                       #panel.grid.minor.x = element_blank(),
                       axis.text.x = element_text(angle = 45,hjust = 1 ),
                       #panel.grid=element_blank(),
                       legend.position = "top",
                       legend.text = element_text(size= 20),
                       legend.title= element_text(size= 40),
                       axis.text.y = element_text(margin = margin(0,0,0,0.2,'cm')),
                       axis.ticks.length.y = unit(0.25, "cm"),
                       #axis.ticks = element_line(linewidth = 2,size = 2),
                       legend.direction = "vertical",
                       legend.key.size = unit(0.5, "cm"),
                       axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0),size=40)
) }
