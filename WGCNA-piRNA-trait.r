####### WGCNA #############
library(WGCNA)
library(reshape2)
library(stringr)
library(flashClust)
library(GO.db)
options(stringsAsFactors = FALSE)
# 打开多线程
enableWGCNAThreads()
setwd('D:/项目/精浆cfRNA/分析结果/Exp/WGCNA/Result/piRNA/')

exprMat <- "D:/项目/精浆cfRNA/分析结果/Exp/Filter/piRNA_TPM.txt"
type = "unsigned"
#type = 'signed'
#corType = "pearson"
corType = 'bicor'
corFnc = ifelse(corType=="pearson", cor, bicor)
maxPOutliers = ifelse(corType=="pearson",1,0.05)
#关联样品性状的二元变量时设置
robustY = ifelse(corType=="pearson",T,F)
##导入数据
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T, 
                       quote="", comment="", check.names=F)
dim(dataExpr)
dataExpr <- dataExpr[,1:423]
#dataExpr <- log2(dataExpr[,1:423]+1)
head(dataExpr)[,422:423]

###筛选中位绝对偏差前80%的基因，至少MAD大于0.01
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                                max(quantile(m.mad, probs=seq(0, 1, 0.2))[2],0.01)),]
###转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))


### 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)
if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
head(dataExpr)[,1:8]
# 随机选取4000列
#selected_cols <- sample(ncol(dataExpr), 4000)
# 生成新的矩阵
#dataExpr <- dataExpr[, selected_cols]
## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average")
#pdf('QC-outline.pdf')
pdf(file="QC-outliner.pdf", onefile=F, paper="special", 
    width=20, height=14, bg="white", pointsize=6)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
#绘制线显示剪切
#abline(h = 40000, col = "red");
dev.off()
################### No need to filter sample
#确定线下的群集
#clust = cutreeStatic(sampleTree, cutHeight = 40000, minSize = 10)
#table(clust)
# clust1包含我们要保留的样本
#keepSamples = (clust==1)
#dataExpr = dataExpr[keepSamples, ]
#nGenes = ncol(dataExpr)
#nSamples = nrow(dataExpr)
#sampleTree = hclust(dist(dataExpr), method = "average")
#pdf(file="QC-filtered.pdf", onefile=F, paper="special", 
#    width=20, height=14, bg="white", pointsize=6)
#plot(sampleTree, main = "Sample clustering after QC", sub="", xlab="")
#dev.off()
#读入表型数据
trait <- read.table("D:/项目/精浆cfRNA/分析结果/Exp/WGCNA/20230416_zq_精浆WGCNA/inputdata/clinical.txt",header = T)
#intersect with Exp matrix
Samples = rownames(dataExpr);
traitRows = match(Samples, trait$sample);
datTraits = trait[traitRows, -1];
rownames(datTraits) = trait[traitRows, 1];
collectGarbage()
#plot cluster of samples with correlation with traits
sampleTree2 = hclust(dist(dataExpr), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度
pdf(file="sample-trait.pdf", onefile=F, paper="special", 
    width=30, height=14, bg="white", pointsize=8)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
####### soft thread hold
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers, 
                        networkType=type, verbose=5)
pdf('Power.pdf',width = 10,height = 6)
par(mfrow = c(1,2))
cex1 = 0.9
#横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准R-square=0.85
abline(h=0.85,col="red")
# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
     cex=cex1, col="red")
dev.off()
power = sft$powerEstimate
power
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                 ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                        ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                               ifelse(type == "unsigned", 6, 12))       
                 )
  )
}
cor <- WGCNA::cor
type='unsigned'
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = 'unsigned', minModuleSize = 20,networkType = type,
                       reassignThreshold = 0, mergeCutHeight = 0.3,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = paste0(exprMat, ".tom"),
                       verbose = 3)
cor <- stats::cor
table(net$colors)
## 灰色的为**未分类**到模块的基因。
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
dynamicColors <- labels2colors(net$unmergedColors)
pdf('Gene_dendrogram_module.pdf')
plotDendroAndColors(net$dendrograms[[1]], cbind(dynamicColors,moduleColors),
                    c("Dynamic Tree Cut", "Module colors"),
                    dendroLabels = FALSE, hang = 0.03,)
dev.off()

###合并模块
MEList = moduleEigengenes(dataExpr, colors = moduleColors)
Mes = MEList$eigengenes
MEDiss = 1 - cor(Mes)
METree = hclust(as.dist(MEDiss), method = "average")
pdf('module_cluster.pdf')
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
dev.off()
MEDissThres = 0.25
abline(h = MEDissThres, col = "red")
merge = mergeCloseModules(dataExpr, moduleColors, cutHeight = MEDissThres, verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf('Gene_dendrogram_merged_module.pdf')
plotDendroAndColors(net$dendrograms[[1]], cbind(moduleColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), 
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)
dev.off()
####save module and genes
gene_module <- data.frame(ID=colnames(dataExpr), module=moduleColors)
gene_module = gene_module[order(gene_module$module),]
write.table(gene_module,file=paste0("TPM",".gene_module.xls"),
            sep="\t",quote=F,row.names=F)

### module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
#save result
MEs_colt = as.data.frame(t(MEs_col))
colnames(MEs_colt) = rownames(dataExpr)
write.table(MEs_colt,file=paste0("FPKM.log2",".module_eipgengene.xls"),
            sep="\t",quote=F)
# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
pdf(file="plotEigengeneNetworks.pdf", onefile=F, paper="special", 
  bg="white", pointsize=6)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
### select hub genes
hubs = chooseTopHubInEachModule(dataExpr, colorh=moduleColors, power=power, type=type)
write.table(hubs,file=paste0("FPKM.log2",".module_hubgengene.xls"),
            sep="\t",quote=F)

#############################################################
con <- nearestNeighborConnectivity(dataExpr, nNeighbors=50, power=power,
                                   type=type, corFnc = 'cor')
load(net$TOMFiles[1], verbose=T)
TOM <- as.matrix(TOM)
dissTOM = 1-TOM
plotTOM = dissTOM^7
diag(plotTOM) = NA


###TOM PLOT 行列同时做层级聚类
TOMplot(plotTOM, net$dendrograms, moduleColors, 
        main = "Network heatmap plot, selected genes")
probes = colnames(dataExpr)
dimnames(TOM) <- list(probes, probes)


###基因网络
cyt = exportNetworkToCytoscape(TOM,
                               edgeFile = paste(exprMat, ".edges.txt", sep=""),
                               nodeFile = paste(exprMat, ".nodes.txt", sep=""),
                               weighted = TRUE, threshold = 0,
                               nodeNames = probes, nodeAttr = moduleColors)
#trait <- "C:/Users/ZD123/Desktop/WGCNA/clinical.txt"
trait <- "D:/项目/精浆cfRNA/分析结果/Exp/WGCNA/20230416_zq_精浆WGCNA/inputdata/clinical.txt"
#trait2 <- "D:/项目/精浆cfRNA/分析结果/Exp/WGCNA/20230416_zq_精浆WGCNA/inputdata/group.txt"

# 读入表型数据
trait <-"D:/项目/精浆cfRNA/分析结果/Exp/WGCNA/20230416_zq_精浆WGCNA/inputdata/clinical.txt"

if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}
#表型和模块相关性聚类
MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
pdf(file="plotEigengeneNetworks_trait.pdf", onefile=F, paper="special", 
    bg="white", pointsize=6)
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
                      marDendro = c(10,10,2,4),
                      marHeatmap = c(10,10,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)
dev.off()
### 模块与表型数据关联
#corType='pearson'
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}

# signif表示保留几位小数
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 2), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
pdf(file="module_trait_relationship.pdf", onefile=F, paper="special", 
    pointsize=10)
#pdf(file="module_group_relationship.pdf", onefile=F, paper="special", 
 #   pointsize=10)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()
## 模块内基因与表型数据关联
### 计算模块与基因的相关性矩阵
if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue = geneModuleMembershipA$p
}
# 计算性状与基因的相关性矩阵
## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。
if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}


geneTraitCorMelt = as.data.frame(geneTraitCor)
write.table(geneTraitCorMelt,file=paste0("FPKM",".gene_trait_correlation.xls"),
            sep="\t",quote=F)
geneTraitCorMelt$ID = rownames(geneTraitCor)
geneTraitCorMelt = melt(geneTraitCorMelt)
colnames(geneTraitCorMelt) <- c("Gene","Trait","PersonCorrelationValue")
geneTraitPMelt = as.data.frame(geneTraitP)
write.table(geneTraitPMelt,file=paste0("FPKM",".gene_trait_correlationPvalue.xls"),
            sep="\t",quote=F)
geneTraitPMelt$ID = rownames(geneTraitP)
geneTraitPMelt = melt(geneTraitPMelt)
colnames(geneTraitPMelt) <- c("Gene","Trait","Pvalue")
#geneTraitCorP = cbind(geneTraitCorMelt, Pvalue=geneTraitPMelt$Pvalue)
geneTraitCorP = merge(geneTraitCorMelt, geneTraitPMelt, by=c("Gene","Trait"))
write.table(geneTraitCorP,
            file=paste0("FPKM",".gene_trait_correlationPvalueMelt.xls"),
            sep="\t",quote=F,row.names=F)
geneTraitCorColor <- numbers2colors(geneTraitCor)
pdf(file = "Cluster_dendrogram_total.pdf")
plotDendroAndColors(net$dendrograms[[1]],
                    cbind(dynamicColors,moduleColors,geneTraitCorColor),
                    c("Dynamic Tree Cut", "Module colors", colnames(geneTraitCor)),
                    dendroLabels = FALSE, hang = 0.5,
                    addGuide = TRUE, guideHang = 0.05,marAll = c(1, 7, 3, 1))
dev.off()
# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
module = "green"
pheno = "survival_rate"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
# 获取模块内的基因
moduleGenes = moduleColors == module
sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


###筛选brown模块中hub基因
geneModuleMembership<-as.data.frame(geneModuleMembership)
hub<- abs(geneModuleMembership$MEgreen)>0.8 & abs(geneTraitCor$survival_rate)>0.6
table(hub)
hub<-as.data.frame(dimnames(data.frame(dataExpr))[[2]][hub]) 
write.csv(hub, "hubgene_GSMM_green_survival_rate.csv")
##HubGenes <- chooseTopHubInEachModule(dataExpr,moduleColors)
##write.csv (HubGenes,file = "TopHubGenes_of_each_module.csv",quote=F)
##HubGenes



###去除grey模块
restGenes= (moduleColors != "grey")
diss1=1-TOMsimilarityFromExpr(dataExpr[,restGenes], power = power)
colnames(diss1) =rownames(diss1) =SubGeneNames[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )
pdf('Gene_dendrogram_module_nogrey.pdf')
plotDendroAndColors(hier1, moduleColors[restGenes], "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, 
                    addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
dev.off()

save.image(file='WGCNA_trait_FPKM.RData')

###模块在样本中的表达量
module_colors= setdiff(unique(moduleColors), "grey")
for (color in module_colors){
  module=SubGeneNames[which(moduleColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,
              quote=FALSE)
}
module.order <- unlist(tapply(1:ncol(dataExpr),as.factor(moduleColors),I))
m<-t(t(dataExpr[,module.order])/apply(dataExpr[,module.order],2,max))
heatmap(t(m),zlim=c(0,1),col=gray.colors(100),Rowv=NA,Colv=NA,labRow=NA,scale="none",
        RowSideColors=moduleColors[module.order])
