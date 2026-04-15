library('ggplot2')
library('pheatmap')
library('colorspace')



#########   Figure S1 : Quality control of all samples ############
#PCA of samples
smRNA_plots <- function(RNA_type)
{
  #read exp matrix
  tpm_path <- paste(RNA_type,'_TPM.txt',sep='')
  meta_path <- "Sample_info.xlsx"
    # 📊 读取表达矩阵（TPM）
  tpm_data <- read_table(tpm_path, show_col_types = FALSE)
  gene_ids <- tpm_data[[1]]
  expr_data <- as.data.frame(tpm_data[, -1])
  rownames(expr_data) <- gene_ids
  colnames(expr_data) <- colnames(tpm_data)[-1]
    # 转置表达矩阵
  expr_t <- as.data.frame(t(expr_data))
  expr_t$sample <- rownames(expr_t)
   # 📋 读取样本注释
  meta <- read_excel(meta_path)
  meta <- meta %>% filter(sample %in% expr_t$sample)
  merged <- left_join(expr_t, meta, by = "sample")
  exp <- merged[, 1:(ncol(merged)-2)]
##PCA
  pca <- prcomp(log2(exp+1))
  sum=as.data.frame(summary(pca)$importance)
  pc1=round(sum$PC1[2]*100,2)
  pc2=round(sum$PC2[2]*100,2)
  PCA = data.frame(group=merged$Group,sam=merged$sample,PC1=pca$x[,1],PC2=pca$x[,2])
#  write.csv(PCA,paste(RNA_type,"noCTL__PCA.txt",sep="_"),row.names=T)
  mytheme <- theme_bw() + 
    theme(axis.text.x=element_text(size=rel(1.3),colour="black",face="bold"), 
          axis.text.y=element_text(size=rel(1.3),colour="black",face="bold"),
          axis.title = element_text(size=rel(1.5),colour="black",face="bold"),
          legend.text = element_text(size=rel(1.5),colour="black"),
          plot.title = element_text(hjust = 0.5,size=rel(1.5),face="bold"),
          panel.grid = element_blank(),
          strip.text = element_text(size = rel(1.5)))
  p=ggplot(PCA,aes(x=PC1,y=PC2,color=factor(group,level=c('NOR','AZS','OLI','AZO','CTL')))) + mytheme +  
    geom_hline(yintercept=0,linetype="dashed",color="grey",size=1)+
    geom_vline(xintercept=0,linetype="dashed",color="grey",size=1)+
    geom_point(size=2) + 
    scale_color_brewer(palette = "Set1")+
    xlab(paste("PC1 (", pc1,"% var explained)",sep = "")) + 
    ylab(paste("PC2 (", pc2,"% var explained)",sep = "")) + 
    ggtitle(RNA_type) +
    labs(colour="",shape="")
#PCA plot
  pdf(paste(RNA_type,"_PCA.pdf",sep="_"),width = 6,height = 5)
  print(p)
  dev.off()
}
smRNA_plots(RNA_type = 'tRNA')
smRNA_plots(RNA_type = 'miRNA')
smRNA_plots(RNA_type = 'piRNA')
smRNA_plots(RNA_type = 'mRNA')
smRNA_plots(RNA_type = 'lncRNA')

#desity of clean_reads
pdf('FigS1B.Density_clean_reads.pdf',pointsize = 20,width = 6,height = 4)
QC_stat <- data.frame(read.table('TableS1_Statistics of sequencing data and quality control.txt',header = T,row.names = 1))
ggplot(QC_stat, aes(x = log10(clean_reads))) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 50, 
                 fill = "lightblue", 
                 alpha = 0.5) +
  geom_density(color = "red", linewidth = 1)+
  theme_classic()
dev.off()

#############FigS2. Profile of miRNA in healthy individuals
#violin of marker miRNAs 
tpm_path <- 'D:/项目/精浆cfRNA/分析结果/Exp/Filter/miRNA_TPM.txt'
tpm_path <- 'TableS3. Origin score of CIBERSOTx using scRNA data of testis.xlsx'
meta_path <- "Sample_info.xlsx"
meta_path <- "TableS2. Clinical information of samples_debug.xlsx"
# 🧾 异常样本排除列表
outlier_samples <- c("CF-S094","CF-S121","CF-S127","CF-S144")
# 📊 读取表达矩阵（TPM）
#tpm_data <- read_table(tpm_path, show_col_types = FALSE)
tpm_data <- read_xlsx(tpm_path)
gene_ids <- tpm_data[[1]]
expr_data <- as.data.frame(tpm_data[, -1])
rownames(expr_data) <- gene_ids
colnames(expr_data) <- colnames(tpm_data)[-1]
# 删除异常样本
expr_data <- expr_data[, !colnames(expr_data) %in% outlier_samples]
# 转置表达矩阵
expr_t <- as.data.frame(t(expr_data))
expr_t$sample <- rownames(expr_t)
# 📋 读取样本注释
meta <- read_excel(meta_path)
meta <- meta %>% filter(sample %in% expr_t$sample)
merged <- left_join(expr_t, meta, by = "sample")

#expr_data$sample <- rownames(expr_data)
#meta <- meta %>% filter(sample %in% expr_data$sample)
#merged <- left_join(expr_data, meta, by = "sample")
#merged_order <- merged[order(merged$number),]
#write.table(merged_order,file='Origin_score.txt',quote = F,sep="\t")

Normal <- merged %>%
  filter(Group %in% c("NOR"))
exp_N <- as.data.frame(log2(Normal[,1:2212]+1))
marker<-data.frame(RPM=c(exp_N$`hsa-miR-1246`,exp_N$`hsa-miR-449a`,exp_N$`hsa-miR-210-3p`,exp_N$`hsa-miR-135a-5p`,exp_N$`hsa-miR-141-5p`,
                         exp_N$`hsa-miR-34c-5p`,exp_N$`hsa-miR-6075`),
                   miRNA=c(rep("hsa-miR-1246",101),rep("hsa-miR-449a",101),rep("hsa-miR-210-3p",101),rep("hsa-miR-135a-5p",101),rep("hsa-miR-141-5p",101),
                           rep("hsa-miR-34c-5p",101),rep("hsa-miR-6075",101)))
pdf('miRNA_marker_violin.pdf',height = 4)
p <- ggplot(marker, aes(x=miRNA, y=RPM,fill=miRNA)) + 
  geom_violin(scale = 'width')
p+theme_classic() + theme(axis.text.x = element_blank()) +
  labs(y='log2(TPM+1)')
dev.off()

## FigS2B. seqlogo of seed region(top300)
library(ggseqlogo)
dat<-read.table("miRNA_seed",header = T)
seqs_miRNA <- as.character(dat$family)
dat<-read.table("miRNA_N_true_seed.txt",header = T)
seqs_miRNA_N <- as.character(dat$Normal)

pdf('miRNA_seqlogo.pdf')
ggseqlogo(seqs_miRNA, method="prob")+
  theme_classic()+
  ggtitle('miRNAs')
ggseqlogo(seqs_miRNA_N, method="prob")+
  theme_classic()+
  ggtitle('miRNAs')
dev.off()

#########   Figure S3 :expression of merged tRNA ###########
tRNA <- data.frame(read.table('tRNA_merge_TPM.txt',header = T,check.names = FALSE),check.names = FALSE)
tRNA_exp <- tRNA[, 1:(ncol(tRNA)-3)]
par(mar=c(8,4,2,2))
pheatmap(t(log2(tRNA_exp+1)),fontsize_col = 8,
         file='tRNA_heatmap_merge_v2.pdf',height = 16,width = 10,
         show_rownames = T,show_colnames = F,scale='column'
         ,cluster_cols = F
         ,treeheight_row = 25
         ,gaps_col = c(101,229,322)
)

##############FigS4. Flowchart generated by Biorender

##########   FigS5. Molecular signatures of altered piRNAs
len_5base <- data.frame(read.table('R_DE_len_5base.txt',header = T,row.names = 1))
pdf('piRNA_len_5base_R_UP.pdf',width = 8,height = 5,pointsize = 10)
barplot(t(data.matrix(len_5base))/1000/126,ylab='Average number of piRNA(X1000)',col=c('red','blue','green','orange'))
legend("topleft",legend =c('G','U','C','A'),cex=0.8,pch=15,col=c('orange','green','blue','red'))
dev.off()
len_5base <- data.frame(read.table('R_down_len_5base.txt',header = T,row.names = 1))
pdf('piRNA_len_5base_R_DOWN.pdf',width = 8,height = 5,pointsize = 10)
barplot(t(data.matrix(len_5base))/1000/126,ylab='Average number of piRNA(X1000)',col=c('red','blue','green','orange'))
legend("topleft",legend =c('G','U','C','A'),cex=0.8,pch=15,col=c('orange','green','blue','red'))
dev.off()


##########    Figure S6-S8. WGCNA of miRNA and piRNA 
#### details in WGCNA code

#heatmap of important module, FigS6B & FigS8B
module_heatmap_ave <- function(data,RNA_type,module_name)
{
  exp <- data.frame(t(data))
  groups = rep(c("N","R","S","W"),times=c(101,128,93,97))
  exp$group <- groups
  #pheatmap(log2(exp+1),scale = 'column',fontsize =8,fontsize_col = 2,show_rownames = T,show_colnames = T,cluster_cols = F)
  #draw heatmap using mean value for each group
  mir_ave <- aggregate(exp[,1],by=list(type=exp$group),mean)[,2]
  for (i in seq(2,dim(data)[1]))
  {
    mir_ave = cbind(mir_ave,aggregate(exp[,i],by=list(type=exp$group),mean)[,2])
  }
  colnames(mir_ave) <- rownames(data)
  pheatmap(t(log2(mir_ave+1)),file=paste(RNA_type,module_name,'_module_ave_heatmap.pdf',sep="_"),height = 4,width = 4,scale = 'row',fontsize = 6,show_rownames = T,show_colnames = T,cluster_cols = F)
}
module <- data.frame(read.table('D:/项目/精浆cfRNA/分析结果/Exp/WGCNA/Result/miRNA/green_module_tpm.txt',header = T,row.names = 1))
module_heatmap_ave(module,RNA_type='miRNA',module_name = 'green')
module <- data.frame(read.table('D:/项目/精浆cfRNA/分析结果/Exp/WGCNA/Result/piRNA/paleturquoise_module_tpm.txt',header = T,row.names = 1))
module_heatmap_ave(module,RNA_type='piRNA',module_name = 'paleturquoise')


##########   Figure S9. boxplot of origin score of all cell types
origin <- data.frame(read_xlsx('TableS3. Origin score of CIBERSOTx using scRNA data of testis.xlsx'))
origin_score <- origin %>%
  mutate(across(2:13, as.numeric))
origin_score$group <- factor(origin$Group, 
                   levels = c('NOR','AZS','OLI','AZO'))
pdf('FigS8_origin_score_celltype_bxp_CR_fpkm.pdf',width=16,height = 12)
par(mfrow=c(3,4))
for (i in seq(2,13))
{
  scores <- data.frame(score = origin_score[,i],group = origin_score$group)
  boxplot(score ~ group,data=scores,pch='.',names=c('NOR','AZS','OLI','AZO'),main=colnames(origin_score)[i],ylab='Score')
}
dev.off()
#wilcox test of origin scores by cell type between groups
for (i in seq(2,13))
{
  #p <- wilcox.test(as.numeric(origin[1:101,i]),as.numeric(origin[102:229,i]))$p.value
  #p <- wilcox.test(as.numeric(origin[1:101,i]),as.numeric(origin[230:322,i]))$p.value
  p <- wilcox.test(as.numeric(origin[1:101,i]),as.numeric(origin[323:419,i]))$p.value
  print(colnames(origin)[i]);print(p)
}

############FigS10 Venn diagram generated from online tool

#############Figure S11 barchart of shap value of features
setwd('D:/项目/精浆cfRNA/分析结果/建模/模型文章结果/num_article')
library(ggplot2)
library(dplyr)
value <- data.frame(read.table('num_feature_importance_byclass_top20.txt',header = T))
# 基础正负条形图
pdf('FigS10A.Num_feature_top20_barchart.pdf',width=6,height = 3)
ggplot(value, aes(x=reorder(Feature,Mean,decreasing=T), y = Mean,fill = Mean < 0)) +
  geom_col() +
  scale_fill_manual(
    values = c("FALSE" = "#E74C3C", "TRUE" = "#3498DB") )+
  labs(title = "",x='',y='Average SHAP value') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),legend.position="none")
dev.off()

setwd('D:/项目/精浆cfRNA/分析结果/建模/模型文章结果/pr_article/')
library(ggplot2)
library(dplyr)
value <- data.frame(read.table('PR_feature_importance_byclass_top20.txt',header = T))
# 基础正负条形图
pdf('FigS10C.PR_feature_top20_barchart.pdf',width=6,height = 3)
ggplot(value, aes(x=reorder(Feature,Mean,decreasing=T), y = Mean,fill = Mean < 0)) +
  geom_col() +
  scale_fill_manual(
    values = c("FALSE" = "#E74C3C", "TRUE" = "#3498DB") )+
  labs(title = "",x='',y='Average SHAP value') +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),legend.position="none")
dev.off()

#### Violin of certain feature
Genes_TPM <- data.frame(read.table("Genes_TPM.txt",header = T,row.names = 1))
groups = rep(c("NOR","AZS","OLI","AZO"),times=c(101,128,93,97))

miRNA_TPM <- data.frame(read.table('miRNA_TPM.txt',header = T,row.names = 1))
exp_gene <- as.data.frame(t(log2(Genes_TPM[,1:423]+1)))
exp_gene$group <- as.factor(groups)
pdf('WT1_violin.pdf',height = 4)
p <- ggplot(exp_gene, aes(x=group, y=WT1,fill=group)) + 
  geom_violin(scale = 'width')
p+theme_classic() + theme(axis.text.x = element_blank()) +
  labs(y='log2(TPM+1)') +scale_x_discrete(limits=c("NOR", "AST", "OLI","AZO"))
dev.off()
pdf('KDM2B_violin.pdf',height = 4)
p <- ggplot(exp_gene, aes(x=group, y=KDM2B,fill=group)) + 
  geom_violin(scale = 'width')
p+theme_classic() + theme(axis.text.x = element_blank()) +
  labs(y='log2(TPM+1)') +scale_x_discrete(limits=c("NOR", "AST", "OLI","AZO"))
dev.off()
pdf('SPATA42_violin.pdf',height = 4)
p <- ggplot(exp_gene, aes(x=group, y=SPATA42,fill=group)) + 
  geom_violin(scale = 'width')
p+theme_classic() + theme(axis.text.x = element_blank()) +
  labs(y='log2(TPM+1)') +scale_x_discrete(limits=c("NOR", "AST", "OLI","AZO"))
dev.off()
pdf('SPANXA2_violin.pdf',height = 4)
p <- ggplot(exp_gene, aes(x=group, y=SPANXA2,fill=group)) + 
  geom_violin(scale = 'width')
p+theme_classic() + theme(axis.text.x = element_blank()) +
  labs(y='log2(TPM+1)') +scale_x_discrete(limits=c("NOR", "AZS", "OLI","AZO"))
dev.off()
pdf('KIAA1210_violin.pdf',height = 4)
p <- ggplot(exp_gene, aes(x=group, y=KIAA1210,fill=group)) + 
  geom_violin(scale = 'width')
p+theme_classic() + theme(axis.text.x = element_blank()) +
  labs(y='log2(TPM+1)') +scale_x_discrete(limits=c("NOR", "AZS", "OLI","AZO"))
dev.off()
pdf('KAT6A_violin.pdf',height = 4)
p <- ggplot(exp_gene, aes(x=group, y=KAT6A,fill=group)) + 
  geom_violin(scale = 'width')
p+theme_classic() + theme(axis.text.x = element_blank()) +
  labs(y='log2(TPM+1)') +scale_x_discrete(limits=c("NOR", "AZS", "OLI","AZO"))
dev.off()

exp_gene <- as.data.frame(cbind(t(log2(Genes_TPM[,1:423]+1)),t(log2(miRNA_TPM[,1:423]+1))))
exp_gene$group <- as.factor(groups)

ggplot(exp_gene, aes(x=`hsa-miR-125a-3p`, y=YBX2,color=group)) + 
  geom_point()
ggplot(exp_gene, aes(x=`hsa-miR-125a-3p`, y=`CDKN2B-AS1`,color=group)) + 
  geom_point()
p <- ggplot(exp_gene, aes(x=group, y=`hsa-miR-125a-3p`,fill=group)) + 
  geom_violin(scale = 'width')
p+theme_classic() + theme(axis.text.x = element_blank()) +
  labs(y='log2(TPM+1)') +scale_x_discrete(limits=c("NOR", "AZS", "OLI","AZO"))
p <- ggplot(exp_gene, aes(x=group, y=XIST,fill=group)) + 
  geom_violin(scale = 'width')
p+theme_classic() + theme(axis.text.x = element_blank()) +
  labs(y='log2(TPM+1)') +scale_x_discrete(limits=c("NOR", "AZS", "OLI","AZO"))

