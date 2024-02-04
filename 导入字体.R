library(extrafont)
font_import(pattern = "simsun.ttf")  
#loadfonts(device = "all")
font_import(pattern = "simsun.ttf")  # 这里以 Times New Roman 为例，你可以根据需要替换为新罗马字体的名称
loadfonts()

windowsFonts(Arial = windowsFont("arial"))
#导入字体进入pdf
#showtext包可给定字体文件，加载到 R环境中，生成新的字体家族名字，后期调用这个名字设定字体，并且支持中文写入pdf不乱码

library(showtext)
## 载入黑体
font_add("simhei", regular = "simhei.ttf")
font_add("simsun", regular = "simsun.ttc")
showtext_auto(enable = TRUE, record = TRUE)
font_families()

plot(1:3, ann = F, axes = F, type = "n", frame.plot = T)
text(2,2.5, "Times New Roman", 
     family = "serif", font = 2)
text(2,2, "Arial", 
     family = "sans", font = 3, cex=5)
text(2,1.5, "Courier New", 
     family = "mono", font = 4)

font_paths("C:\\Windows\\Fonts")

## 开始使用
showtext_begin()
text(2,2.5, "宋体", 
     family = "ST", font = 2)
text(2,2, "黑体", 
     family = "heiti", font = 3)
