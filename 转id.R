library(GEOquery)
'
#转id
GPL <- gset[[1]]@annotation  # 平台信息——提取芯片平台编号
GPL
http://www,bio-info-trainee.com/1399.html#里面找对应的平台
BiocManager::install("illuminaHumanv3") 
BiocManager::install("hgu133plus2.db") 
'''
BiocManager::install("hgu133plus2.db") 
library(hgu133plus2.db) #GPL570
#大致查看R包中的信息，寻找我们需要的SYMBOL
ls("package:hgu133plus2.db")
#用toTable（）函数提取
ids <- toTable(hgu133plus2SYMBOL) # 提取
head(ids)  # 查看提取内容 

#install.packages("dplyr")
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")
exp<-as.data.frame(exp)
exp$probe_id=rownames(exp)  # 将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
exp2<- merge(exp,ids,by.x="probe_id", by.y="probe_id")  # 合并数据
exp3<- exp2[!duplicated(exp2$symbol),]  # 按照symbol列去重
# 数据框probe_exp的行名变成symbol
rownames(exp3)=exp3$symbol
exp3=exp3[,c(-1,-ncol(exp2))]
#输出文件
write.table(exp3,file = "ids_exprs.txt",sep = "\t",row.names=T,col.names = T)
write.csv(exp2,file = "ids_exprs.csv") 



#直接读gpl平台
GPL=getGEO(filename = "GPL1211_family.soft.gz") 
gpl=GPL@dataTable@table
colnames(gpl) 
ids=gpl[,c(1,4)]
"
我只要ID和symbol，这两列在1，11列，所以下面的数字为1，11
每次你自己可以先看看gpl这个矩阵，那一列是ID和symbol,这个不固定
ids=gpl[,c(1,11)]
输出文件
"""
write.table(ids,file = "ids.txt",sep = "\t",row.names=F,col.names = T)
write.csv(ids,file = "ids.csv")
"
通过上述方法我们得到了ID对应的symbol,
接下来只需要将symbol列对应进矩阵就好了 ,下面的可以直接运
如果打不开这个包就下载一
"""
#install.packages("dplyr")
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")
exp=as.data.frame(exp)
exp$probe_id=rownames(exp)  # 将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id")  # 合并数据
# 按照symbol列去重
exp3=exp2[!duplicated(exp2$symbol),]
# 数据框probe_exp的行名变成symbol
rownames(exp3)=exp3$symbol
exp3=exp3[,c(-1,-ncol(exp2))]
#输出文件
write.table(exp3,file = "ids_exprs.txt",sep = "\t",row.names=T,col.names = T)
write.csv(exp3,file = "ids_exprs.csv") 
'
方法三、通过 bitr（）函数进行转化
有些GSE的ID实可恨，既没有对应的R包，官网的平台数据也没有对应的symbol
（例如GSE42872平台信息没有对应的symbol），寻寻觅觅，
只能放出终极大招——通过bitr（）函数进行转化
强制将ID转成我们需要的symbol格式
！注意，这个方法不一定可以将全部的ID转化成gene symbol
'''
rm(list = ls())
下载数据
library(GEOquery)
eSet <- getGEO("GSE42872",
               destdir = '.',
               getGPL = F)
exp <- exprs(eSet[[1]])  # 表达矩阵
首先打开这两个R包
library(clusterProfiler)
library(org.Hs.eg.db)
如果没有这两个R包，就运下面的代码进行下载
if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
查看数据类型（就是这个R包提供哪些个ID类型可供转化）
keytypes(org.Hs.eg.db)
下面这个是R包org.Hs.eg.db拥有的ID类型，可供选择，对应原来的ids里面的类型
ID的格式，你挑一个出来和下面的是对应的

#[1] "ACCNUM" "ALIAS" "ENSEMBL" "ENSEMBLPROT" "ENSEMBLPROT"

#"ENTREZID"

#[7] "ENZYME" "EVIDENCE" "EVIDENCEALL" "GENENAME" "GO" "GOALL"

#[13] "IPI" "MAP" "OMIM" "ONTOLOGY" "ONTOLOGYALL" "PATH"

#[19] "PFAM" "PMID" "PROSITE" "REFSEQ" "SYMBOL" "UCSCKG"
确保数据格式为数据框

ids1=rownames(exp)
ids1=as.data.frame(ids1)
colnames(ids1)="ID"
注意这一步一定会有红字出现的，说有百分几几fail，只要百分之10就是对的，50%就谢天谢地了
就可以继续下去
ids <- bitr(ids1$ID,  # 你的数据框
            fromType = "PMID",   # 你的ID的数据类型
            toType = c("SYMBOL","ENSEMBL"),  # 转化的数据类型
            OrgDb = org.Hs.eg.db)  # org.Hs.eg.db——人类
#最后得出的ids就是结果
输出文件
write.table(ids,file = "ids.txt",siep = "\t",row.names=F,col.names = T)
write.csv(ids,file = "ids.csv")
通过上述方法我们得到了ID对应的symbol,

接下来只需要将symbol列对应进矩阵就好了 ,下面的可以直接运

如果打不开这个包就下载一下

#install.packages("dplyr")
library(dplyr)
colnames(ids) = c("probe_id" ,"symbol")
exp=as.data.frame(exp)  # 转化为数据框格式
exp$probe_id=rownames(exp)  # 将行名变为列名为probe_id的一列
# exp是原来的表达矩阵
exp2= merge(exp,ids,by.x="probe_id", by.y="probe_id")  # 合并数据
# 按照symbol列去重
exp2=exp2[!duplicated(exp2$symbol),]
# 数据框probe_exp的行名变成symbol
rownames(exp2)=exp2$symbol
exp2=exp2[,c(-1,-ncol(exp2))]
#输出文件
write.table(exp2,file = "ids_exprs.txt",sep = "\t",row.names=T,col.names = T)
write.csv(exp2,file = "ids_exprs.csv")
#最后提醒一下,方法一二不分上下，方法三是在一二没有办法下才用的