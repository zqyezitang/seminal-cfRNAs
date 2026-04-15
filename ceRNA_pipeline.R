### ---------------
###
### Create: qyl
### Date: 2023-09-11 16:09:00
### Email: 524583919@qq.com
### Function：
### - Constructing ceRNA network according to mRNA and lncRNA expression matrix
### - 
###
### Update Log: 2023-09-11
###
### ---------------

#### ----- Function ------ ####
ceRNA <- function(mRNA_exp,miRNA_exp,
                  group,case="case",ctrl="ctrl",
                  grp_nm = "GO_KEGG",dir_nm = "ceRNA"){
  
  # 需要修改源代码中部分代码
  
  ### 1.library ####
  suppressMessages({
    library(tidyverse)
    library(GDCRNATools) 
    library(clusterProfiler) # enrichGo/enrichKEGG/bitr
    library(org.Hs.eg.db) # 人类gene ID
  })
  output_dir <- paste0("./outputdata/",dir_nm,"/",grp_nm)
  photo_dir <- paste0("./photo/",dir_nm,"/",grp_nm)
  dir.create(output_dir,recursive = T)
  dir.create(photo_dir,recursive = T)
  
  ### 2.normalization ####
  rnaCounts2 <- mRNA_exp %>% 
    round() %>% 
    na.omit() %>% 
    filter(rowSums(.)>0) %>% 
    mutate(symbol=rownames(.),.before=1)
  rnaCounts2_ENSG <- rnaCounts2 %>% 
    filter(grepl("^ENSG.*",symbol))
  rnaCounts2_sym <- rnaCounts2 %>% 
    filter(!grepl("^ENSG.*",symbol))
  ent_3 <- bitr(rnaCounts2_sym$symbol, 
                fromType = 'SYMBOL', 
                toType = 'ENSEMBL',
                OrgDb = org.Hs.eg.db)
  rnaCounts <- merge(ent_3,rnaCounts2_sym,by.y="symbol",by.x="SYMBOL") %>% 
    dplyr::select(-1) %>% 
    dplyr::rename(symbol=ENSEMBL) %>% 
    rbind(rnaCounts2_ENSG) %>% 
    group_by(symbol) %>% 
    summarise(across(everything(),mean)) %>% 
    column_to_rownames("symbol") %>% 
    round() %>% 
    na.omit() %>% 
    filter(rowSums(.)>0)
  
  rnaExpr <- gdcVoomNormalization(counts = rnaCounts, filter = FALSE)
  mirExpr <- gdcVoomNormalization(counts = miRNA_exp, filter = FALSE)
  
  ### 3.DEGs ####
  DEGAll <- gdcDEAnalysis(
    counts     = rnaCounts, 
    group      = group, # 必须：组名
    comparison = paste0(case,"-",ctrl),# 肿瘤比正常
    method     = 'DESeq2', # 建议DESeq2
    filter     = F # 只用cpm>1
  )
  
  deLNC <- gdcDEReport(
    deg = DEGAll, 
    gene.type = 'long_non_coding' 
  )
  dePC <- gdcDEReport(deg = DEGAll, gene.type = 'protein_coding')
  
  ### 4.ceRNA ####
  for (i in c('starBase','spongeScan','miRcode')) {
    print(i)
    tryCatch({
      ceOutput <- gdcCEAnalysis(
        lnc         = rownames(deLNC), # 行名必须是ENSG格式，基因太少会报错
        pc          = rownames(dePC), 
        deMIR       = NULL, # 可以给定差异的miRNA基因向量
        lnc.targets = i, # 'spongeScan', 'starBase', and 'miRcode'
        pc.targets  = i, # 同上
        rna.expr    = rnaExpr, # 经过voom转换的表达矩阵，用于判断各个样本的三者间表达趋势
        mir.expr    = mirExpr
      ) 
      
      ### 5.output ####
      ceOutput2 <- ceOutput[ceOutput$hyperPValue<0.01 & 
                              ceOutput$corPValue<0.01 & ceOutput$regSim != 0,]
      edges <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'edges')
      nodes <- gdcExportNetwork(ceNetwork = ceOutput2, net = 'nodes')
      table(nodes$type)
      write.table(edges, 
                  file=paste0(output_dir,"/",i,"_edges.txt"), 
                  sep='\t',row.names = T,col.names = NA, quote=F)
      write.table(nodes, 
                  file=paste0(output_dir,"/",i,"_nodes.txt"), 
                  sep='\t', row.names = T,col.names = NA,quote=F)
    },
    error=function(e) print(e))
  }
}

#### ----- running analysis ------ ####

### 1.ceRNA ####
#source(file = "E:/Bioinformatic_project/02_R/code_R/A_Script_Function/pipeline/ceRNA_pipeline.R")
mlncRNA_count <- read.table(file = "./inputdata/Exp_count/mlncRNA_count.txt",
                            sep = "\t",row.names = 1,header = T,check.names = F)
mi_count <- read.table(file = "./inputdata/Exp_count/miRNA_count.txt",
                       sep = "\t",row.names = 1,header = T,check.names = F)
mlncRNA_count <- mlncRNA_count[,1:(ncol(mlncRNA_count)-6)]
mi_count <- mi_count[,1:(ncol(mi_count)-6)]
(colnames(mlncRNA_count)==colnames(mi_count)) %>% 
  table()

grp_df <- read.table(file = "./inputdata/TableS2. Clinical information of samples_debug.txt", # test or debug
                     sep = "\t",header = T,row.names = 1)
grp_df <- grp_df %>% 
  mutate(Group = factor(Group,levels = c("NOR", "AZS", "OLI", "AZO"))) %>% 
  arrange(Group) %>% 
  mutate(Group = as.character(Group))
intersect_saps <- intersect(rownames(grp_df),colnames(mlncRNA_count))
mlncRNA_count <- mlncRNA_count[,intersect_saps,drop=F]
mi_count <- mi_count[,intersect_saps,drop=F]
grp_df <- grp_df[intersect_saps,,drop=F]
grp <- grp_df$Group

for(j_num in 2:length(unique(grp))){
  # j_num <- 4
  tryCatch({
    j_grp <- unique(grp)[j_num]
    index_N <- grep("NOR",grp)
    index_j <- grep(j_grp,grp)
    j_mlnc <- mlncRNA_count[,c(index_j,index_N)]
    j_mi <- mi_count[,c(index_j,index_N)]
    j_grps <- grp[c(index_j,index_N)]
    ceRNA(mRNA_exp = j_mlnc,miRNA_exp = j_mi,group = j_grps,
          case = j_grp,ctrl = "NOR",grp_nm = j_grp,dir_nm = "ceRNA")
  },
  error=function(e) print(e))
}
