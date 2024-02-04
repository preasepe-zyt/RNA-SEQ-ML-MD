
#1.1 读取数据
library(WGCNA)
library(xlsx)
library(tidyverse)
enableWGCNAThreads()#查看WGCNA可以调用多少线程工作
## Allowing parallel execution with up to 7 working processes.
options(stringsAsFactors = FALSE);
#Read in the female liver data set
datExpr0 = as.data.frame(t(ids_exprs));

#1.2 检查是否包含缺失值
gsg = goodSamplesGenes(datExpr0, verbose = 3);

#反之需要删除有问题的样品与基因
if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
    # Remove the offending genes and samples from the data:
    datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
datExpr0[1:5,1:5]

#1.3 聚类 观察是否有离群值
sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
pdf(file="聚类和极端组排除.pdf",height=10,width=20,family = "Arial")
par(cex = 0.6);
par(mar = c(0,4,2,0),mfrow=c(2,1))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#1.4 去除聚类高度高于15的样本
# Plot a line to show the cut

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)+abline(h = 160, col = "red");
dev.off()
clust = cutreeStatic(sampleTree, cutHeight = 160, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
datExpr[1:5,1:5]

#1.5 加载注释信息
traitData <- pdata2
dim(traitData)#有三十多列注释
names(traitData)
traitData$...1 -> ID
traitData[,-1] -> traitData
traitData$samples <- ID
allTraits <- traitData
# remove columns that hold information we do not need.

# Form a data frame analogous to expression data that will hold the clinical traits.
femaleSamples = rownames(datExpr);#取出样本名
traitRows = match(femaleSamples, allTraits$samples);
datTraits = allTraits[traitRows,];
rownames(datTraits) = datTraits$samples
case_when(datTraits$`gender:ch1`=="Female"~ "1",
          datTraits$`gender:ch1`=="Male"~"2") -> datTraits$`gender:ch1`
case_when(datTraits$`state:ch1`=="Healthy"~"1",
          datTraits$`state:ch1`=="Depression"~"2") -> datTraits$`state:ch1`
datTraits[,-5] -> datTraits
apply(datTraits, 2, as.numeric)-> datTraits
as.data.frame(datTraits) -> datTraits
names(datTraits) <- c("Gender","Age","Depression","T_cells_CD8")
collectGarbage();

sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry

traitColors = numbers2colors(datTraits, signed = FALSE);
#给注释中的数字变量赋予颜色，感觉与我的hcluster热图无异
# Plot the sample dendrogram and the colors underneath.
pdf(file = "sampleClustering_related_info.pdf", width = 5, height = 5,family = "Arial");
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
#构建network、识别module

#2.1 自动化分块式

#2.1.1 确定power（β）


#Choosing the soft-thresholding power: analysis of network topology
#Constructing a Depressioned gene network entails the choice of the soft thresholding power 
# to which co-expression similarity is raised to calculate adjacency
powers = c(c(1:10), seq(from = 12, to=20, by=2))#给定一些值用于测试
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

pdf(file = "soft_threshhlod.pdf", width = 6, height = 6,family = "serif");
# Plot the results:
#sizeGrWindow(9, 5)
# Scale-free topology fit index as a function of the soft-thresholding power
par(family = "serif",cex.axis = 2,cex.lab = 1.5, font.lab = 2,cex.main = 1.5)  
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"),cex=10);
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=2,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.98,col="red")
text(sft$fitIndices[,1][4], -sign(sft$fitIndices[,3][4])*sft$fitIndices[,2][4],col="blue",cex=2, labels=4)
dev.off()

pdf(file = "Mean_connectivity_threshhlod.pdf", width = 6, height = 6,family = "serif");
par(family = "serif",cex.axis = 2,cex.lab = 1.5, font.lab = 2,cex.main = 1.5)  
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=2,col="red")
dev.off()
#要选择lowest power for which the scale-free topology index curve  attens out upon reaching a high value
sft$powerEstimate
## [1] 6
#也就是说这里要选择6
#2.1.2 一步法网络构建
net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       TOMType = "unsigned", 
                       minModuleSize = 30,#最小模块的基因数
                       reassignThreshold = 0, 
                       mergeCutHeight = 0.25,#threshold for merging of modules
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "femaleMouseTOM",
                       verbose = 3,
                       maxBlockSize=ncol(datExpr))
#对应内存需求如下
#16GB workstation should handle up to 20000 probes; 
#32GB workstation should handle perhaps 30000. 
#4GB standard  8000-10000 probes,
#这一步参数众多，需要的话可以?一下
#返回值中net$colors 包含了模块的分配情况
#net$MEs 包含模块的特征基因.
table(net$colors)#看看这些基因被分为了哪些模块

#注意，0是不属于任何模块的基因，也就是后面画出灰色的部分，其余模块颜色按模块基因从大到小排列
net$dendrograms[[1]]#hierarchical clustering dendrogram存放位置

mergedColors = labels2colors(net$colors)#把数字返回成颜色
# Plot the dendrogram and the module colors underneath
pdf(file="聚类树0.1.pdf",height=5,width=12,family = "serif")
par(family = "serif",cex.axis = 2,cex.lab = 3, font.lab = 2,cex.main = 3)  
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    cex.colorLabels = 1.2)
dev.off()
#2.1.3 修改、存储
#recutBlockwiseTrees可以用于重新计算，比从头计算节省时间
# 存储文件用于下游分析
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,net,
     file = "FemaleLiver-02-networkConstruction-auto.RData")


#2.2 步进法 有助于自定义各类参数

#2.2.1 还是一样要确定β

powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# Plot the results:

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))+
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red")+abline(h=0.90,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))+
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#2.2.2 adjacency
adjacency = adjacency(datExpr, power = sft$powerEstimate)
adjacency[1:5,1:5]

#2.2.3 Topological Overlap Matrix (TOM)
#去噪、防止假关联的出现将adjacency转换为Topological Overlap Matrix（反应similarity）, 并用1-TOM反应dissimilarity:
TOM = TOMsimilarity(adjacency);

dissTOM = 1-TOM
TOM[1:5,1:5]


# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average")#利用dissTOM聚类更快
# Plot the resulting clustering tree (dendrogram)

plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)#这种聚类树每个leaf代表一个基因

#2.2.4 色块聚类树
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 30)#每个模块最小的基因数量

table(dynamicMods)

dynamicColors = labels2colors(dynamicMods)# Convert numeric lables into colors
table(dynamicColors)

# Plot the dendrogram and colors underneath
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

#2.2.5 合并相似表达的模块
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEs[1:5,1:5]


MEDiss = 1-cor(MEs)# Calculate dissimilarity of module eigengenes,注意计算的是不详细性
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")+abline(h=0.25, col = "red")

#height cut 是 0.25, 也就是说对应的相关性是 0.75以上的模块用于merge
#
####开始merge
# Call an automatic merging function
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)

#merge为返回的新的WGCNA对象
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)#查看合并前后的模块聚类图

#2.2.6 保存核心文件用于下游分析
moduleColors = mergedColors# Rename to moduleColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1;
MEs = mergedMEs;
# Save module colors and labels for use in subsequent parts
save(MEs, moduleLabels, moduleColors, geneTree, file = "FemaleLiver-02-networkConstruction-stepByStep.RData")

#2.2.7 在机器配置不允许的情况下如何处理大型矩阵？
#### 2.3 
#还是要先确定power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))+
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         labels=powers,cex=cex1,col="red")+abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))+text(sft$fitIndices[,1], sft$fitIndices[,5], 
                                             labels=powers, cex=cex1,col="red")

#测试数据中仅包含3600个探针，而实际情况中可能会有超过五万个探针的矩阵出现
#本方法并不是真正的处理这么大的矩阵，而是采取两级聚类的方法：
#先用一个crude clustering method将基因划分至一些block中，而后对每个block中的基因进行完整的network分析
#最后再用2.2中谈到的merge功能将相似的模块拟合起来，
#优点：运行速度大大加快并节省内存，缺点:由于第一次聚类方法比较简单，离群的基因可能也会被分配到block中
bwnet = blockwiseModules(datExpr, maxBlockSize = 2000,
                         power = 3, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "femaleMouseTOM-blockwise",
                         verbose = 3)

table(net$colors)
bwLabels = matchLabels(bwnet$colors, moduleLabels)#把这次的label和2.2中的merge后的label合并 
bwModuleColors = labels2colors(bwLabels)# Convert labels to colors for plotting
table(bwLabels)

#最大的处理能力是2000，那么这些基因在这里会把3600个基因划分为2个block
# Plot the dendrogram and the module colors underneath for block 1
plotDendroAndColors(bwnet$dendrograms[[1]], bwModuleColors[bwnet$blockGenes[[1]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 1",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

# Plot the dendrogram and the module colors underneath for block 2
plotDendroAndColors(bwnet$dendrograms[[2]], bwModuleColors[bwnet$blockGenes[[2]]],
                    "Module colors", main = "Gene dendrogram and module colors in block 2",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#把这两个画在一张聚类树上
plotDendroAndColors(geneTree,
                    cbind(moduleColors, bwModuleColors),
                    c("Single block", "2 blocks"),
                    main = "Single block gene dendrogram and module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

#把两次的特征基因取出来
singleBlockMEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
blockwiseMEs = moduleEigengenes(datExpr, bwModuleColors)$eigengenes
single2blockwise = match(names(singleBlockMEs), names(blockwiseMEs))
signif(diag(cor(blockwiseMEs[, single2blockwise], singleBlockMEs)), 3)


#探寻各模块与表型间的相关性

# Define numbers of genes and samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs0[1:5,1:5]

MEs = orderMEs(MEs0)
MEs[1:5,1:5]

moduleTraitCor = cor(MEs, datTraits, use = "p")
#他这里是加入了临床数据，是数值型的可以做相关性分析，但是更多的时候加入的是分组变量
moduleTraitCor[1:5,1:5]

moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
moduleTraitPvalue[1:5,1:5]


textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
textMatrix = paste(signif(moduleTraitCor, 2))
textMatrix[1:5,1:5]#绘制热图时也可以参考这一展示方式
pdf(file="临床相关.pdf",height=8,width=7,family = "serif")
par(mar = c(8, 10, 3, 4),family = "serif")
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(100),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"),
               cex.lab.x = 1,
               cex.lab.y = 1)#这个图
dev.off()
## Warning in greenWhiteRed(50): WGCNA::greenWhiteRed: this palette is not suitable for people
## with green-red color blindness (the most common kind of color blindness).
## Consider using the function blueWhiteRed instead.

#计算感兴趣的特征，比如这里的Depression与各个模块基因的关联性
Depression = as.data.frame(datTraits$Depression)
names(Depression) = "Depression"
Depression[1:5,]
## [1] 38.0 33.5 33.9 44.3 32.9
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
geneModuleMembership[1:5,1:5]

names(MMPvalue) = paste("p.MM", modNames, sep="")
MMPvalue[1:5,1:5]

geneTraitSignificance = as.data.frame(cor(datExpr, Depression, use = "p"))
geneTraitSignificance[1:5,]
## [1] -0.06937286 -0.07525085  0.24342785  0.17121953 -0.21325726
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(Depression), sep="")
names(GSPvalue) = paste("p.GS.", names(Depression), sep="")


#通过Gene Significance(GS)和Module Membership(MM)测量，我们可以识别出体重与体重significance高的基因以及模块特征基因之间的关系
#下面的例子是在棕色模块中绘制了significance基因与模块成员关系的散点图
module = "turquoise"
column = match(module, modNames);
moduleGenes = moduleColors==module;

pdf(file="模块和AD的相关.pdf",height=6,width=6,family = "Arial")

verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Depression",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()

#模块表达图
pdf(file="模块和AD的相关.pdf",height=6,width=6,family = "Arial")
barplot(MEs[,21],col="turquoise",cex.axis=2,cex.names=3,
        ylab="Eigengene Expression",xlab="Samples(Control-Depression)", main ="turquoise Module")
dev.off()

##### 分析结果输出 #####
names(datExpr)[moduleColors=="turquoise"]#会返回包含在棕色模块中的探针名
datExpr[,moduleColors=="turquoise"] %>%
    t() -> turquoise
datExpr[,tan] -> tan
write.xlsx(turquoise,"turquoise.xlsx")
annot = read.csv(file = "./testdata/GeneAnnotation.csv");
dim(annot)

names(annot)

probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
## [1] 0
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$gene_symbol[probes2annot],
                       LocusLinkID = annot$LocusLinkID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
geneInfo0[1:5,]

# Order modules by their significance for Depression
modOrder = order(-abs(cor(MEs, Depression, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
    oldNames = names(geneInfo0)
    geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]]);
    names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.Depression));
geneInfo = geneInfo0[geneOrder, ]

write.csv(geneInfo, file = "geneInfo.csv")

04

#主要是把模块中的基因提取出来做富集分析

# Read in the probe annotation
annot = read.csv(file = "./testdata/GeneAnnotation.csv");
# Match probes in the data set to the probe IDs in the annotation file
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# Get the corresponding Locuis Link IDs
allLLIDs = annot$LocusLinkID[probes2annot];
# $ Choose interesting modules
intModules = c("brown", "red", "salmon")
for (module in intModules)
{
    # Select module probes
    modGenes = (moduleColors==module)
    # Get their entrez ID codes
    modLLIDs = allLLIDs[modGenes];
    # Write them into a file
    fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
    write.table(as.data.frame(modLLIDs), file = fileName,
                row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)
后面就是富集分析的内容，大家可以参考下面这两篇： 
基因集富集分析
富集分析App

本文参考

‍https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/
    
    文件下载链接

https://pan.baidu.com/s/1ZyjcJ0qkRNwZtMA07s9Tig?pwd=enoh
t(tan) -> tan
