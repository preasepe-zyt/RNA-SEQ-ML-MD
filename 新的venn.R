install.packages("ggvenn")
library(ggvenn)
library(tidyverse)
library(extrafont)
library(readxl)
loadfonts(device = "pdf")
#导入字体进入pdf
#showtext包可给定字体文件，加载到 R环境中，生成新的字体家族名字，后期调用这个名字设定字体，并且支持中文写入pdf不乱码
library(showtext)
showtext_auto(enable=TRUE)
turquoise <- read.csv("C:/Users/79403/Desktop/时瑞蝶/交集和kegg-go/turquoise.txt", row.names=1)
diff <- read.delim("C:/Users/79403/Desktop/时瑞蝶/交集和kegg-go/diff.txt", row.names=1)
targets_quchong <- read_excel("targets_quchong.xlsx")
anti_inflammatory <- read_excel("抗炎gene.xls")
genes_diff <- rownames(diff[diff$Group!="Not-Sig",])
targets_quchong$x -> d_targets
dat <- list( Arisaematis_Rhizoma = d_targets,  Depression = genes_diff)

pdf("天南星和差异.pdf",width = 8,height = 8,family = "serif")
p <- ggvenn(dat,show_percentage = T,
         stroke_color = "white",
         stroke_size = 0.5,
         fill_color = c("#00BFC8","#B7AE50"),
         #set_name_color = c("#00BFC8","#B7AE50"), 
         set_name_size = 8,
         text_size=8)
p[["layers"]][[1]][["constructor"]][["family"]] <- "simhei.ttf"
p[["layers"]][[2]][["constructor"]][["family"]] <- "simhei.ttf"
p[["layers"]][[3]][["constructor"]][["family"]] <- "simhei.ttf"
p[["layers"]][[4]][["constructor"]][["family"]] <- "simhei.ttf"
p$theme$text$family <- "simhei.ttf"
p
dev.off()

tian_deg <- intersect(d_targets,genes_diff)
dat2 <- list( Inflammatory = anti_inflammatory$Symbol,  Turquoise = colnames(turquoise), Arisaematis_Rhizoma = d_targets,DEGs = genes_diff)
pdf("venn.pdf",width = 8,height = 6,family = "Arial")
ggvenn(dat2,
       show_percentage = T,
       stroke_color = "white",
       stroke_size = 1,
       fill_color = c("#BF5B17","#CC78A6","darkgreen",'#0796E0'),
       set_name_color =c("#BF5B17","#CC78A6","darkgreen",'#0796E0'), 
       set_name_size = 7,
       text_size=7)
dev.off()
reduce(dat2,intersect) -> final_genes
write.csv(fianl_genes,"final_genes.txt")
dat3 <- list( DEGs = genes_diff, anti_inflammatory = 抗炎gene$Symbol)
pdf("DEGS和anti-inflammatory.pdf",width = 8,height = 8,family = "Arial")
ggvenn(dat3,show_percentage = T,
       stroke_color = "white",
       stroke_size = 0.5,
       fill_color = c("#BC80BD","darkgreen"),
       set_name_color =c("#BC80BD","darkgreen"), 
       set_name_size = 5,
       text_size=4)
dev.off()



library(tidyverse)
library(ggvenn)

# Assuming 'dat' is a list and you want to convert it to a data frame
# If 'dat' is already a data frame, you can skip this step
dat <- list_to_data_frame(dat)

# Assuming your data frame has columns named "Arisaematis_Rhizoma" and "Depression"
# You can adjust the column names based on your actual data frame structure
pdf("天南星和差异.pdf",width = 8,height = 8,family = "serif")
ggplot(dat) +
    geom_venn(aes(A = Arisaematis_Rhizoma, B = Depression), 
              stroke_size = 0.5,
              stroke_color = "white",
              fill_color = c("#00BFC8","#B7AE50"),
              set_name_color = c("#00BFC8","#B7AE50"), 
              set_name_size = 8,
              text_size=8) +
    coord_fixed() +
    theme_void()
# 修改文本字体和大小
library(ggplot2)
library(grid)
library(gridExtra)
grid.newpage()
grid.draw(venn)

for (i in 1:4) {
    grid.text(
        label = gsub("arial", "Times Bold", grid.get(sprintf("label-%d-1", i))$label),
        gp = gpar(fontsize = 14, fontface = "bold"),
        vp = grid.get(sprintf("label-%d-1", i))
    )
}
dev.off()

'
3.图形美化
填充
fill_color：填充颜色
fill_alpha：填充透明度
边框
stroke_color：边框颜色
stroke_alpha：边框透明度
stroke_size：边框粗细
stroke_linetype：边框线的类型
集合名
set_name_color：集合名颜色
set_name_size：集合名字号
集合内文本
text_color：文本颜色
text_size：文本字号
百分比
show_percentage：TRUE or FALSE

'''