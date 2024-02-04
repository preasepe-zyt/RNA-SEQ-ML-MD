library(tidyverse)
library(reshape2)
# 创建一个 2x2 的列联表 
data <- matrix(c(19, 210, 1407, 0), nrow = 2) 
colnames(data) <- c("In Ingredient", "Not In Ingredient") 
rownames(data) <- c("In Disease", "Not In Disease")

# 执行费舍尔精确检验 
fisher.test(data)  

# 或者指定表格和检验方向的替代语法 
fisher.test(x = data, alternative = "two.sided")
#plot
mosaicplot(data,
           main = "Mosaic plot",
           color = TRUE
)
#检验期望次数
chisq.test(data)$expected
#卡方检验
test <- chisq.test(data)
test

# load packages
library(ggstatsplot)
library(ggplot2)
# plot
#pre deal
herb <- c(rep("In Intersection",19),rep("Not In Intersection",210))
disease <- c(rep("In Intersection",19),rep("Not In Intersection",1407))
group <- c(rep("Arisaematis Rhizoma",229),rep("Major Depressive Disorder",1426))
dat <- data.frame( Conditions=c(herb,disease),group)
#卡方的
ggbarstats(
    data = dat,
    x =  Conditions,
    y = group
) +
    labs(caption = NULL)+ # remove caption
   mytheme
library(vcd)

mosaic(~ group + Conditions,
       direction = c("v", "h"),
       data = dat,
       shade = TRUE
)
#fisher的
pdf("fisher.pdf", width=8, height=9,family = "serif")
ggbarstats(
    dat, 
    Conditions, 
    group,
    results.subtitle = FALSE,
    sample.size.label.args = list(size = 10),
    subtitle = paste0(
        "Fisher's exact test", ", p-value = ",
        ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
    )
)+mytheme
dev.off()

if(T){mytheme <- theme(text = element_text(family = "serif"),
                     plot.subtitle = element_text(size = 20,hjust = 0.5,color="black",face = "bold"),
                     axis.title = element_text(size = 20,color="black",face = "bold"), 
                     axis.text = element_text(size= 20,color="black",face = "bold"),
                     panel.grid.minor.y = element_blank(),
                     panel.grid.minor.x = element_blank(),
                     axis.text.x = element_text(angle = 45, hjust = 1,size= 20,color="black",face = "bold"),
                     panel.grid=element_blank(),
                     legend.position = "top",
                     legend.text = element_text(size= 20,color="black",face = "bold"),
                     legend.title= element_text(size= 20,color="black",face = "bold"),
                     title = element_text(size = 20,hjust = -2,color="black",face = "bold")
    ) }
