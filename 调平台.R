library(devtools)
#安装注释探针的R包
install_github("jmzeng1314/AnnoProbe")
#加载AnnoProbe这个包
library(AnnoProbe)
#选择要注释的探针类型
gpl='GPL570'

#得到探针对应的基因名字
ids=idmap(GPL,type = 'pipe')
#展示前10条结果
head(probe2gene)
