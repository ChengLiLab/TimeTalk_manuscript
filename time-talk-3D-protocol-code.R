###-----log-----
###Project: tele-embryo
### test ZGA datasets
### author wlt
### 20200830
### 20200831
### 20200831:v0.0.2
### 20200931
### 20200902-04
### 20200907
### 20200908:v0.0.3
### 20200909:v0.0.4
### 20200911:v0.0.5 don't lost your lane and main purpose
### 20200913:v0.0.5 add mature tissue sample(fibroblasts)
### 20200914:v0.0.6 don't be overflow
### 20200916:v0.0.7 add GO plot 
### 20200918:v0.0.8 organize your results
### 20200919-22:v0.0.9 generate cellphonedb input
### 20200923:v0.0.1 pheatmap
### 20200924-0929: v0.1.0 add reprogramming
### 20200930:v0.1.1 reformat your code
### 20201001-20201003: v.0.1.2
### 20201006:plot,add protein coding gene
### 20201007-20201012: massive,chaos;
### 20201013:v0.1.3, step by step, calm down;
### 20201014-1023:v0.1.4, step by step, from QC metric;
### 20201026-1102:v0.1.5, update LR genelist by ConnectDB 2.0
### 20201103-1108:v0.1.6, filter eLR
### 20201109-1113:v0.1.7, add zga trends
### 20201113:v.0.1.8, pretty figure
### 20201125:v.0.1.9 add atac-seq anlysis; coordination;
### 20201129-30:v.0.1.893, coordination
### 20201201:v0.1.91: the loser now will later to win 
### 20201210:v.0.1: rename to tele-embryo, and formate the code
### 20201214:v.0.1: move forward, forward
### 20201215-20201225: v0.1.1 the road to convince the audience;
### 20210108: v0.1.2 investigate eLR function
### 20210111: v0.1.2.2 update correlation failed
### 20210114: v0.1.3 update for correlation and genes, eg Bmp4
### 20210201: v0.1.4 use nichenetr!
###-----keep you on R------
for(ii in 1:1000000){
  cat(paste("round",ii),sep = "\n")
  Sys.sleep(ii/10)
}

####----0.Load package and function-------

####----0.1 Load package-------
library(tidyverse)
#library(SingleCellExperiment)
library(ggplot2)
###library(LRBase.Mmu.eg.db)
library(cowplot)
#library(ggpubr)
library(ggbeeswarm)
#library(ggridges)
#library(ggalluvial)
library(ggpubr)
library(rstatix)
library(pheatmap)
library(VennDiagram)
library(Seurat)
library(UpSetR)
library(clusterProfiler)
library(RDAVIDWebService)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(factoextra)
library(energy)
library(RcisTarget)
library(igraph)
library(ggnetwork)
##library(ESEA)
library(RColorBrewer)
library(msigdbr)
library(clusterProfiler)
library(GOSemSim)
library(nichenetr)
library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)
source("code/myUtils.R")


####----0.2 define color------
rdbu <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu"))) 
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))
col.spectral <- colorRampPalette(brewer.pal(11,'Spectral')[-6])
test.color.3 <- colorRampPalette(c("#f86e11","#e9970a","#71a701","#62b474","#71c3ac","#9fc4ca"))
rdwhbu <- colorRampPalette(c("navy", "white", "brown3"))
pheatmap(volcano,color = rev(col.spectral(100)),
         border_color = NA,
         cluster_cols = F,
         cluster_rows = F)

RColorBrewer::display.brewer.all()
?pheatmap

###-----1. prepare data------
####GSE38495
###-----1.0 prepare Oocyte------
tmp.1 <- read_delim(file = "data/GSM967487_MII_oocyte_CASTEi_1_RPKM.txt.gz",delim = "\t",skip = 1) %>%
  dplyr::rename(Oocyte1 = `RPKM/FPKM`,Gene_symbol=`Gene symbol`) %>%
  dplyr::select(Gene_symbol,Oocyte1)

tmp.2 <- read_delim(file = "data/GSM967488_MII_oocyte_CASTEi_2_RPKM.txt.gz",delim = "\t",skip = 1) %>%
  dplyr::rename(Oocyte2 = `RPKM/FPKM`,Gene_symbol=`Gene symbol`) %>%
  dplyr::select(Gene_symbol,Oocyte2)

tmp.3 <- read_delim(file = "data/GSM967489_MII_oocyte_CASTEi_3_RPKM.txt.gz",delim = "\t",skip = 1) %>%
  dplyr::rename(Oocyte3 = `RPKM/FPKM`,Gene_symbol=`Gene symbol`) %>%
  dplyr::select(Gene_symbol,Oocyte3)


####check the symbol coordination
sum(tmp.1$Gene_symbol==tmp.2$Gene_symbol)
sum(tmp.2$Gene_symbol==tmp.3$Gene_symbol)

Oocyte.df <- cbind(tmp.1,Oocyte2=tmp.2$Oocyte2,Oocyte3=tmp.3$Oocyte3) 
Oocyte.df <- Oocyte.df[!duplicated(Oocyte.df$Gene_symbol),]

#### fix gene name, eg. Bmp4|Gm15217	
Oocyte.df$Gene_symbol <- gsub(pattern = "\\|.*",replacement = "",x = Oocyte.df$Gene_symbol)



saveRDS(Oocyte.df,file = myFileName(prefix = "res/R/Oocyte_GSE38495_raw_rpkm_",suffix = ".rds"))


###-----1.1 prepare Deng et.al 2013---------

###-----1.1.1 meta data-----
dir <- "data/test/GSE45719_RAW/"
samplefiles <- list.files(dir)

# tmp.idx <- which(unlist(lapply(strsplit(samplefiles,split = "_"),length))==3)
# samplefiles[tmp.idx]
# 
# tmp.idx <- which(unlist(lapply(strsplit(samplefiles,split = "_"),length))==5)
# samplefiles[tmp.idx]

sampleinfor <- gsub(pattern = "_expression.txt","",x = samplefiles)
sampleinfor <- gsub(pattern = "GSM[0-9]{7}_",replacement = "",sampleinfor)
GSM <- regmatches(x = samplefiles,m = regexpr(pattern = "GSM[0-9]{7}",text = samplefiles))
GSM

sample.df <- data.frame(samplefiles=samplefiles,
                        sampleinfor=sampleinfor,
                        GSM=GSM,
                        embryo_id=gsub(pattern = "-[0-9]*",replacement = "",sampleinfor),
                        stringsAsFactors = F)
sample.df$Stage <- gsub(pattern = "_[0-9].*",replacement = "",sample.df$embryo_id)
length(unique(sample.df$sampleinfor))

sample.df$embryo_id <- gsub(pattern = "fibroblast_.*_",replacement = "",sample.df$embryo_id)
sample.df$embryo_id[grep("BxC|CxB",sample.df$embryo_id,perl = T)] <- paste0("fibroblast","_",sample.df$embryo_id[grep("BxC|CxB",sample.df$embryo_id,perl = T)])

table(sample.df$Stage)




###-----1.1.2 get gene expression-----
# ii <- 1
# tmp.df <- read.delim(paste0(dir,samplefiles[ii]),stringsAsFactors = F)
# tmp.df <- tmp.df[,c(1,3)]
# colnames(tmp.df) <- c("Gene_symbol",sampleinfor[ii])
# tmp.df <- tmp.df[!duplicated(tmp.df$Gene_symbol),]
# data.exp <- tmp.df
# for (ii in 2:length(sampleinfor)) {
#   ### show messages during run loops
#   cat("merge sample",ii,"\n")
#   tmp.df <- read.delim(paste0(dir,samplefiles[ii]),stringsAsFactors = F)
#   tmp.df <- tmp.df[,c(1,3)]
#   colnames(tmp.df) <- c("Gene_symbol",sampleinfor[ii])
#   tmp.df <- tmp.df[!duplicated(tmp.df$Gene_symbol),]
#   data.exp <- merge(data.exp,tmp.df,by = "Gene_symbol")
# }
# 
# ###length(unique(data.exp$Gene_symbol))
# rownames(data.exp) <- data.exp$Gene_symbol
# data.exp <- data.exp[,-1]
# gene_symbols <- rownames(data.exp)
# saveRDS(data.exp,file = "res/R/deng_science_2013_rpkm.Rds")


####-------1.2 prepare early gastrulation data-----------
#### don't run this again
# ####-------1.2.1 load data--------
# 
# ### expression data
# early_gastrulation <- read_delim(file = "data/GSE109071_rpkm.txt.gz",delim = "\t")
# early_gastrulation.1 <- read_delim(file = "data/GSE109071_rpkmBlastocyst.txt.gz",delim = "\t")
# #### This step will introduce NA, so should remove NA
# early_gastrulation <- full_join(early_gastrulation,early_gastrulation.1,by="Gene") %>%
#   column_to_rownames("Gene")
# early_gastrulation <- myRemoveNA(early_gastrulation)
# 
# data.exp <- myNormalize(data.exp = early_gastrulation)
# ### remove ERCC
# data.exp <- data.exp[grep(pattern = "^ERCC",rownames(data.exp),invert = T),]
# 
# 
# ###meta data
# test <- GEOquery::getGEO("GSE109071")
# early_gastrulation.meta <- test$GSE109071_series_matrix.txt.gz@phenoData@data
# early_gastrulation.meta <- early_gastrulation.meta %>%
#   mutate(cell_id=gsub(pattern  = "Mouse embryo cell ",replacement = "",title)) %>%
#   mutate(Stage=gsub(pattern  = "embryonic day ",replacement = "E",`age:ch1`)) %>%
#   dplyr::select(cell_id,Stage)
# rownames(early_gastrulation.meta) <- early_gastrulation.meta$cell_id
# 
# 
# 
# ###---1.2.2 run seurat work flow------
# seu <- myCreateSeurat(data = data.exp,meta.data = early_gastrulation.meta)
# levels(seu)
# ### "blastocyst" "EB"    
# seu <- subset(seu,idents ="blastocyst",invert=T)
# saveRDS()
# #### remove ERCC gene
# 
# nrow(seu)
# ncol(seu)
# ###1724 is right
# 
# seu <- ScaleData(seu)
# seu <- FindVariableFeatures(seu,selection.method = "vst",nfeatures = 3000)
# seu <- RunPCA(seu,features = VariableFeatures(seu))
# tmp.mat <- as.matrix(seu@assays$RNA@data)
# tmp.genes <- VariableFeatures(seu)
# seu.dist <- 1- cor(tmp.mat[tmp.genes,],method = "pearson")
# seu <- RunTSNE(seu,
#                distance.matrix = seu.dist,
#                dims = 1:20,
#                seed.use = 42,
#                perplixity=100)
# 
# seu <- FindNeighbors(seu,dims = c(1:20))
# seu <- FindClusters(seu,resolution = 0.1)
# 
# TSNEPlot(seu,label=T)
# 
# 
# 
# 
# 
# 
# #####------1.2.3 Set cluster-----------
# p1 <- TSNEPlot(seu,label=T)
# 
# p2 <- FeaturePlot(seu,features = c("Pou5f1", "Bmp4", "Amn"))
# 
# plot_grid(p1,p2)
# ?SetIdent
# ###Pou5f1+, EPI
# cell.use <- WhichCells(seu,idents = c("0","4","5","6"))
# seu <- SetIdent(seu,cells = cell.use,value = "EPI")
# ###BMP4+,ExE
# cell.use <- WhichCells(seu,idents = c("2"))
# seu <- SetIdent(seu,cells = cell.use,value = "ExE")
# 
# cell.use <- WhichCells(seu,idents = c("1","3","7"))
# seu <- SetIdent(seu,cells = cell.use,value = "VE")
# 
# TSNEPlot(seu,label=T)
# 
# 
# seu$cell.type <- Idents(seu)
# 
# p1 <- TSNEPlot(seu,group="Stage")
# p1
# p2 <- FeaturePlot(seu,features = c("Pou5f1", "Bmp4", "Amn"))
# p3 <- TSNEPlot(seu,label=T)
# 
# plot_grid(p3,p2)
# ggsave(filename = myFileName(prefix = "res/fig/early_gastru_tSNEPlot",
#                              suffix = ".jpg"),
#        width = 12,height = 6, dpi = 350)
# 
# 
# #####----1.2.4 save Reuslts-----
# seu$cell_type <- Idents(seu)
# saveRDS(seu,file = myFileName(prefix = "res/R/early_gastrulation",suffix = ".rds"))

#####---1.3 load bulk RNA-seq data---------

xiewei.data.exp <- read.delim(file = "data/test/GSE66582_stage_FPKM.txt",stringsAsFactors = F)
xiewei.data.exp <- xiewei.data.exp %>%
  column_to_rownames("gene")
boxplot(xiewei.data.exp)
colnames(xiewei.data.exp)[4:6] <- c("2cell","4cell","8cell")
colnames(xiewei.data.exp)
saveRDS(object = xiewei.data.exp,file = "res/R/xiewei.data.exp_20201012.rds")



#### the above code needn't change

#### 20210114

####-----1.4 merge to get beforeEPI data-------

#### load before exp data
data.exp <- readRDS("res/R/deng_science_2013_rpkm.Rds")
tmp.select <- c("zy1","zy2","zy3","zy4",
                "early2cell_0r-1","early2cell_0r-2",
                "early2cell_1-1","early2cell_1-2", 
                "early2cell_2-1","early2cell_2-2",
                "early2cell_3-1","early2cell_3-2",
                "mid2cell_0r-1","mid2cell_0r-2",  
                "mid2cell_3-1","mid2cell_3-2",
                "mid2cell_4-1","mid2cell_4-2",   
                "mid2cell_5-1","mid2cell_5-2",
                "mid2cell_6-1","mid2cell_6-2",   
                "mid2cell_7-1","mid2cell_7-2",
                "late2cell_5-1","late2cell_5-2",
                "late2cell_6-1","late2cell_6-2",  
                "late2cell_7-1","late2cell_7-2",
                "late2cell_8-1","late2cell_8-2", 
                "late2cell_9-1","late2cell_9-2",
                "4cell_1-1","4cell_1-2",
                "4cell_1-4","4cell_2-1",
                "4cell_2-2","4cell_2-3",
                "4cell_2-4","4cell_3-1",
                "4cell_3-3","4cell_3-4",
                "4cell_4-1","4cell_4-2",
                "4cell_4-3","4cell_4-4",
                "8cell_1-1","8cell_1-2",            
                "8cell_1-4","8cell_1-5",            
                "8cell_1-6","8cell_1-7",
                "8cell_1-8","8cell_2-1",            
                "8cell_2-2","8cell_2-3",            
                "8cell_2-4","8cell_2-6",            
                "8cell_2-7","8cell_2-8",            
                "8cell_5-1","8cell_5-2",            
                "8cell_5-3","8cell_5-4",            
                "8cell_5-6","8cell_5-7",            
                "8cell_5-8","8cell_8-1",            
                "8cell_8-2","8cell_8-3",            
                "8cell_8-4","8cell_8-6",            
                "8cell_8-7","8cell_8-8")

tmp.select <- c(tmp.select,grep("split",grep("16cell",sampleinfor,perl = T,value = T),value = T,invert = T))
tmp.select <- c(tmp.select,grep("earlyblast",sampleinfor,perl=T,value = T))
tmp.select <- c(tmp.select,grep("midblast",sampleinfor,perl=T,value = T))

tmp.select <- c(tmp.select,grep("lateblast",sampleinfor,perl=T,value = T))
tmp.select <- c(tmp.select,grep("fibroblast",sampleinfor,perl = T,value = T))

tmp.sample.df <- data.frame(cell_id=tmp.select,
                            cell_type=tmp.select,
                            embryo_id=gsub(pattern = "-.*",
                                           replacement = "",
                                           gsub(pattern = ".*_",
                                                replacement = "",
                                                tmp.select)),
                            embryo_cell_id = gsub(pattern = ".*-",
                                                  replacement = "",
                                                  tmp.select),
                            Stage=gsub(pattern = "_.*",
                                       replacement = "",
                                       tmp.select),
                            stringsAsFactors = F)
tmp.sample.df$Stage[1:4] <- "zygote"


data.exp <- data.exp[,tmp.select]
beforeEPI.exp <- data.exp %>%
  rownames_to_column("Gene_symbol")


### fix gene id
### 20210114
beforeEPI.exp$Gene_symbol <- gsub(pattern = "\\|.*",replacement = "",beforeEPI.exp$Gene_symbol)


####load Oocyte
Oocyte.exp <- readRDS(file = "res/R/Oocyte_GSE38495_raw_rpkm__20210114.rds")
Oocyte.meta <- data.frame(cell_id=colnames(Oocyte.exp)[-1],
                        cell_type=colnames(Oocyte.exp)[-1],
                        embryo_id=colnames(Oocyte.exp)[-1],
                        embryo_cell_id=colnames(Oocyte.exp)[-1],
                        Stage=rep("Oocyte",3),stringsAsFactors = F)


####merge data
beforeEPI.exp <- inner_join(Oocyte.exp,beforeEPI.exp,by="Gene_symbol")
beforeEPI.exp <- beforeEPI.exp[!duplicated(beforeEPI.exp$Gene_symbol),]

rownames(beforeEPI.exp) <- NULL
beforeEPI.exp <- beforeEPI.exp %>%
  column_to_rownames("Gene_symbol")


beforeEPI.meta <- rbind(Oocyte.meta,tmp.sample.df)



####save to rds

saveRDS(beforeEPI.exp,file = "res/R/beforeEPI_exp_20210114.rds")
saveRDS(beforeEPI.meta,file = "res/R/beforeEPI_meta_20210114.rds")



####------2.QC----------
####-------2.1 load data-----------

####-----2.1.1 load deng et.al data and normalized, get mean expression-----
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)
Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                 "4cell","8cell","16cell","earlyblast","midblast",  
                 "lateblast","fibroblast")
### test if normalized
# head(beforeEPI.meta)
# boxplot(beforeEPI.exp[,110:100])
beforeEPI.exp.mean <- beforeEPI.exp %>%
  rownames_to_column("Gene") %>%
  gather(key = "cell_id",value = "gene.exp",-Gene) %>%
  left_join(beforeEPI.meta,"cell_id") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp)) %>%
  ungroup() %>%
  mutate(Stage=factor(Stage,levels = Stage.level)) %>%
  spread(key="Stage",value = "gene.exp.mean") %>%
  column_to_rownames("Gene")

boxplot(beforeEPI.exp.mean)
saveRDS(beforeEPI.exp.mean,file = "res/R/beforeEPI_exp_mean_20210114.rds")





####------2.1.2 load xiewei data--------
##### RNA-seq xiewei data
xiewei.data.exp <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
xiewei.data.exp <- myNormalize(myRemoveNA(xiewei.data.exp))


####-------2.2 scatter plot--------

####-------2.2.1 prepare for plot-----

beforeEPI.exp.mean <- readRDS(file = "res/R/beforeEPI_exp_mean_20210114.rds")
gene_symbols <- intersect(rownames(xiewei.data.exp),rownames(beforeEPI.exp.mean))

####rename for plot
colnames(beforeEPI.exp.mean)[1] <- "MII_oocyte"
colnames(xiewei.data.exp)[4] <- "late2cell"
colnames(xiewei.data.exp)[3] <- "early2cell"

print(colnames(xiewei.data.exp),quote=T)



####-------revised version-----------
Stage <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell")
data.plot.1 <- beforeEPI.exp.mean[gene_symbols,]
data.plot.2 <- xiewei.data.exp[gene_symbols,]


plot.list <- list()
for(ii in 1:length(Stage)){
  for(jj in 1:length(Stage)){
    df.plot <- data.frame(x=data.plot.1[,Stage[jj]],
                          y=data.plot.2[,Stage[ii]],
                          stringsAsFactors = F)
    res.cor <- cor(df.plot$x,df.plot$y)
    res.p <- cor.test(df.plot$x,df.plot$y)
    p <- ggplot(df.plot,aes(x,y))+
      geom_point(alpha=0.3)+
      xlab(paste0("pseudo-bulk ",Stage[jj]))+
      ylab(paste0("bulk ",Stage[ii]))+
      ggtitle(paste0("R=",round(res.cor,3)))+
      scale_x_continuous(limits = c(0,25))+
      scale_y_continuous(limits = c(0,25))+
      coord_fixed()+
      theme_cowplot()+
      theme(plot.title = element_text(hjust = 0.5))
    plot.list <- c(plot.list,list(p))
  }
}

####Figure S2
pp  <- plot_grid(plotlist = plot.list,nrow = length(Stage))
myggsave(pp,prefix = "res/fig/figS2_QC_beforeEPI_cor_revised",
         suffix = ".jpg",
         width = 12,height = 12,dpi=350)




#####-----pheatmap-----
Stage <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell")

res.mat <- matrix(0,nrow = length(Stage),ncol = length(Stage))
for(ii in 1:length(Stage)){
  for(jj in 1:length(Stage)){
    df.plot <- data.frame(x=data.plot.1[,Stage[jj]],
                          y=data.plot.2[,Stage[ii]],
                          stringsAsFactors = F)
    res.cor <- cor(df.plot$x,df.plot$y)
    res.p <- cor.test(df.plot$x,df.plot$y)
    res.mat[ii,jj] <- res.cor
  }
}
rownames(res.mat) <- paste0("bulk ",Stage)
colnames(res.mat) <- paste0("pseudo-bulk ",Stage)

tmp.stage  <- c("MII oocyte","zygote","early 2-cell","late 2-cell","4-cell","8-cell")

jpeg(filename = myFileName(prefix = "res/fig/Fig1_QC_pheatmap_",suffix = ".jpg"),
     width = 8,height = 8,res = 350,units = "in")
pheatmap(res.mat,fontsize = 16,
         show_rownames = T,
         display_numbers = T,
         number_color = "white",
         fontsize_number = 16,
         labels_row = paste0("bulk ",tmp.stage),
         labels_col = paste0("pseudo-bulk ",tmp.stage),
         show_colnames = T,angle_col = 315,
         color = rdbu(100))
dev.off()

####------2.3 pca---------

####------2.3.1 replicate fig1A deng.et.al-----------

median_center <- function(x){
  res <- x - apply(x, 2, median)
  return(res)
}

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage!="fibroblast") %>%
  pull(cell_id)

data.plot <- data.plot[,tmp.select]
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)

#fviz_eig

eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(beforeEPI.meta,by = "cell_id")

data.plot$Stage[1:3] <- "MII oocyte" 


Stage <- c("MII oocyte","zygote","early2cell","mid2cell","late2cell",
           "4cell","8cell","16cell","earlyblast","midblast",  
           "lateblast")

data.plot$Stage <- factor(data.plot$Stage,levels = Stage)

###edit embryo id
data.plot$embryo_id[1:3] <- 1:3
data.plot$embryo_id[4:7] <- 1:4

data.plot <- data.plot %>%
  mutate_cond(Stage == "early2cell",
              embryo_id = plyr::mapvalues(embryo_id,from = c("0r",1,2,3),to = 1:4)) %>%
  mutate_cond(Stage == "mid2cell",
              embryo_id = plyr::mapvalues(embryo_id,from = c("0r",3:7),to = 1:6)) %>%
  mutate_cond(Stage == "late2cell",
              embryo_id = plyr::mapvalues(embryo_id,from = 5:9,to = 1:5)) %>%
  mutate_cond(Stage == "8cell",
              embryo_id = plyr::mapvalues(embryo_id,from = c(1,2,5,8),to = c(1,2,3,4))) %>%
  mutate_cond(Stage == "16cell",
              embryo_id = plyr::mapvalues(embryo_id,from = c(1,4,5,6),to = c(1,2,3,4))) %>%
  mutate_cond(Stage == "earlyblast",
              embryo_id = plyr::mapvalues(embryo_id,from = 2:4,to = c(1,2,3)))



ggplot(data.plot,aes(PC1,PC2,color=Stage))+
  geom_point(size=3)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(24,25,16,17,18,19))+
  theme_cowplot()
  
ggsave(filename = myFileName(prefix = "res/fig/QC_replicate_deng_fig1A_v1",suffix = ".jpg"),
       width = 8,height = 8,dpi = 300)

### more precise
ggplot(data.plot,aes(PC1,PC2,color=Stage,shape=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  theme_cowplot()

ggsave(filename = myFileName(prefix = "res/fig/QC_replicate_deng_fig1A_v2",suffix = ".jpg"),
       width = 8,height = 8,dpi = 300)


######-------------2.3.2  add xiewei's data-----------------

######-------------2.3.2.1 boxplot-----------

######pseudo-bulk 

Stage <- c("MII_oocyte","zygote","early2cell","mid2cell","late2cell",
           "4cell","8cell","16cell","earlyblast","midblast",  
           "lateblast")

data.plot.1 <- beforeEPI.exp.mean[,Stage] %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = "gene.exp",-Gene) 


data.plot.1$Stage <- factor(data.plot.1$Stage,Stage)

tmp.stage <-  c("MII oocyte","zygote","early 2-cell","mid 2-cell","late 2-cell",
                "4-cell","8-cell","16-cell","early blastocyst","mid blastocyst",  
                "late blastocyst")
p <- myBoxplot_advanced(data = data.plot.1,
                   x = "Stage",
                   y = "gene.exp",
                   fill = "Stage",color = NULL,
                   mycolor = test.color.2(length(Stage)),
                   ylimts = c(0,20),ybreaks = seq(0,20,5))+
  scale_x_discrete(labels=tmp.stage)+
  ylab("gene.exp.mean")+
  theme(legend.position = "none")
p
myggsave(p = p,prefix = "res/fig/figS1_QC_pseudobulk_boxplot",suffix = ".jpg",width = 8,height = 8,dpi = 350)


tmp.stage <- c("MII_oocyte","zygote","early2cell","late2cell",
               "4cell","8cell")
data.plot.2 <- xiewei.data.exp[,tmp.stage] %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = gene.exp,-Gene) %>%
  mutate(group = "bulk")


data.plot.2$Stage <- factor(data.plot.2$Stage,levels = tmp.stage)

tmp.stage <- c("MII oocyte","zygote","early 2-cell","late 2-cell",
               "4-cell","8-cell")
p <- myBoxplot_advanced(data = data.plot.2,
                        x = "Stage",
                        y = "gene.exp",
                        fill = "Stage",
                        color = NULL,
                        mycolor = test.color.2(length(Stage)),
                        ylimts = c(0,20),ybreaks = seq(0,20,5))+
  scale_x_discrete(labels=tmp.stage)+
  ylab("log2(RPKM+1)")+
  theme(legend.position = "none")
p
myggsave(p = p,prefix = "res/fig/figS1_QC_bulk_boxplot",suffix = ".jpg",width = 8,height = 8,dpi = 350)



#####-----2.3.2.2 PCA---------

gene.use <- intersect(rownames(xiewei.data.exp),rownames(beforeEPI.exp.mean))
data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
     filter(Stage!="fibroblast") %>%
     pull(cell_id)
data.plot <- data.plot[gene.use,tmp.select]

tmp.select <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell")
data.plot <- cbind(xiewei.data.exp[gene.use,tmp.select],data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####add meta data for xieweidata
tmp.select <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell")
ttttt <- data.frame(cell_id=tmp.select,
                    cell_type=tmp.select,
                    embryo_id="x",
                    embryo_cell_id = tmp.select,
                    Stage=tmp.select,
                    stringsAsFactors = F)

beforeEPI.meta$embryo_id[1:3] <- 1:3
beforeEPI.meta$embryo_id[4:7] <- 1:4

# beforeEPI.meta <- beforeEPI.meta %>%
#   mutate_cond(Stage == "early2cell",
#               embryo_id = plyr::mapvalues(embryo_id,from = c("0r",1,2,3),to = 1:4)) %>%
#   mutate_cond(Stage == "mid2cell",
#               embryo_id = plyr::mapvalues(embryo_id,from = c("0r",3:7),to = 1:6)) %>%
#   mutate_cond(Stage == "late2cell",
#               embryo_id = plyr::mapvalues(embryo_id,from = 5:9,to = 1:5)) %>%
#   mutate_cond(Stage == "8cell",
#               embryo_id = plyr::mapvalues(embryo_id,from = c(1,2,5,8),to = c(1,2,3,4))) %>%
#   mutate_cond(Stage == "16cell",
#               embryo_id = plyr::mapvalues(embryo_id,from = c(1,4,5,6),to = c(1,2,3,4))) %>%
#   mutate_cond(Stage == "earlyblast",
#               embryo_id = plyr::mapvalues(embryo_id,from = 2:4,to = c(1,2,3)))

ttttt$data.type <- "RNA-seq"
beforeEPI.meta$data.type <- "scRNA-seq"
tmp.meta <- rbind(ttttt,beforeEPI.meta)


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")

data.plot$Stage[c(1,7:9)] <- "MII_oocyte" 

Stage <- c("MII_oocyte","zygote","early2cell","mid2cell","late2cell",
           "4cell","8cell","16cell","earlyblast","midblast",  
           "lateblast")

Stage.1 <- c("MII oocyte","zygote",
           "early 2-cell","mid 2-cell","late 2-cell",
           "4-cell","8-cell","16-cell",
           "early blastocyst","mid blastocyst","late blastocyst")
data.plot$Stage <- plyr::mapvalues(data.plot$Stage,from = Stage,to = Stage.1)

unique(data.plot$Stage)

data.plot$Stage <- factor(data.plot$Stage,levels = Stage.1)


p <- ggplot(data.plot,aes(PC1,PC2,color=Stage,shape=data.type))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))
p
myggsave(p,prefix = "res/fig/Fig1_QC_pca",suffix = ".jpg",width = 8,height = 8,dpi = 350)


####------2.3.3.1 scRNA-seq data merge--------
####beforeEPI
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)
####EPI
seu <- readRDS("res/R/early_gastrulation_20210115.rds")
EPI.exp <- seu@assays$RNA@data
EPI.meta <- seu[[]]
####merge metadata
head(beforeEPI.meta)
head(EPI.meta)
beforeEPI.meta$cell_type <- beforeEPI.meta$Stage
tmp.1 <- beforeEPI.meta %>%
  dplyr::select(cell_id,cell_type,Stage)
tmp.2 <- EPI.meta %>%
  dplyr::select(cell_id,cell_type,Stage)
early.scRNAseq.meta <- rbind(tmp.1,tmp.2)
####merge expression
tmp.1 <- beforeEPI.exp %>%
  rownames_to_column("Gene")
tmp.2 <- EPI.exp %>%
  as.data.frame() %>%
  rownames_to_column("Gene")

###inner_join
early.scRNAseq.exp <- inner_join(tmp.1,tmp.2,by="Gene")
####test na
which(is.na(early.scRNAseq.exp$Gene))
#####test duplicate
which(duplicated(early.scRNAseq.exp$Gene))
early.scRNAseq.exp <- early.scRNAseq.exp %>%
  column_to_rownames("Gene")
early.scRNAseq.exp <- myRemoveNA(early.scRNAseq.exp)

early.scRNAseq.exp["Igf2",1000:1003]

# ####save
# saveRDS(early.scRNAseq.exp,file = "res/R/early.scRNAseq.exp_20210114.rds")
# saveRDS(early.scRNAseq.meta,file = "res/R/early.scRNAseq.meta_20210114.rds")

####save
saveRDS(early.scRNAseq.exp,file = "res/R/early.scRNAseq.exp_20211218.rds")
saveRDS(early.scRNAseq.meta,file = "res/R/early.scRNAseq.meta_20211218.rds")


#####------2.4 add EPI QC--------
seu <- readRDS("res/R/early_gastrulation_20210115.rds")
seu@assays$RNA@counts <- 2^(seu@assays$RNA@data)-1

levels(seu) <- rev(c("EPI","ExE","VE"))
TSNEPlot(seu,label=T)

p1 <- TSNEPlot(seu,group="Stage")



mytmp.plot <- function(tmp.gene){
  p2 <- FeaturePlot(seu,features = tmp.gene,pt.size=6)+
    theme(axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank())+
    scale_color_gradientn(colors = rdwhbu(100))
  p2
  p2.2 <- VlnPlot(seu,
                  features = tmp.gene,
                  slot = "counts")+
    xlab(NULL)+
    ylab("RPKM")+
    ggtitle(NULL)+scale_y_continuous(expand = c(0,0))+
    scale_fill_manual(values = rev(c("#e90b8e","#3b53a5","#6ebe45")))+
    NoLegend()+
    coord_flip()
  p2 <- p2 / p2.2
  return(p2)
}

p2 <- mytmp.plot(tmp.gene = "Pou5f1")
p2
p3 <- mytmp.plot(tmp.gene = "Bmp4")
p3
p4 <- mytmp.plot(tmp.gene = "Amn")
p4

pp <- plot_grid(p2,p3,p4,ncol = 3)
pp


p <- plot_grid(p1+xlab("t-SNE Dimension 1")+ylab("t-SNE Dimension 2")+ggtitle(NULL)+
                 theme(legend.position = "top",
                       axis.line = element_line(size=1),
                       axis.ticks = element_blank(),
                       axis.text  = element_blank()),
               pp,rel_widths = c(0.5,1))
p
myggsave(p , prefix = "res/fig/early_gatrulation_Fig1C",suffix = ".jpg",width = 16,height = 6)

####The code below will be stable in the future version

table(seu$cell_type)

#####-------2.5 all dimension reduction------------

#####-------2.5.1 load data-------------
early.scRNAseq.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20211218.rds")
early.scRNAseq.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20211218.rds")
rownames(early.scRNAseq.meta) <- early.scRNAseq.meta$cell_id

tmp.levels.1 <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                "4cell","8cell","16cell","earlyblast","midblast","lateblast",
                "E5.25","E5.5","E6.25","E6.5","fibroblast")
tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst",
                  "E5.25","E5.5","E6.25","E6.5","fibroblast")
early.scRNAseq.meta$Stage <- plyr::mapvalues(early.scRNAseq.meta$Stage,
                                             from = tmp.levels.1,
                                             to = tmp.levels.2)
early.scRNAseq.meta$Stage <- factor(early.scRNAseq.meta$Stage,
                                    levels = tmp.levels.2)

tmp.levels.1 <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                  "4cell","8cell","16cell","earlyblast","midblast","lateblast",
                  "EPI","ExE","VE","fibroblast")
tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst",
                  "Epiblast(EPI),","Ectoderm(ExE)",
                  "Visceral Endoderm (VE)","fibroblast")
early.scRNAseq.meta$cell_type <- plyr::mapvalues(early.scRNAseq.meta$cell_type,
                                             from = tmp.levels.1,
                                             to = tmp.levels.2)
early.scRNAseq.meta$cell_type <- factor(early.scRNAseq.meta$cell_type,
                                    levels = tmp.levels.2)


#####-------2.5.2 create seurat object---------
seu <- myCreateSeurat(data = early.scRNAseq.exp,
                      meta.data = early.scRNAseq.meta)

seu <- ScaleData(seu)
seu <- FindVariableFeatures(seu,selection.method = "vst",nfeatures = 3000)
seu <- RunPCA(seu,features = VariableFeatures(seu))
tmp.mat <- as.matrix(seu@assays$RNA@data)
tmp.genes <- VariableFeatures(seu)
seu.dist <- 1- cor(tmp.mat[tmp.genes,],method = "pearson")
seu <- RunTSNE(seu,
               distance.matrix = seu.dist,
               dims = 1:20,
               seed.use = 42,
               perplixity=100)

seu <- FindNeighbors(seu,dims = c(1:20))
seu <- FindClusters(seu,resolution = 0.1)
saveRDS(seu,file = myFileName(prefix = "res/R/early.scRNAseq.seurat",suffix = ".rds"))
####---------2.5.3 tsne plot-------------

# tmp.colors <- c("#FCB03B","#0F6937","#67BD45",
#   "#169192","#B09DBA","#954496",
#   "#4E1550","#8EAF3D","#F05153",
#   "#5358A5","#79A7AA","#62BDA6","#4098B6","#4173B3","#5E4FA2","#000000")
p <- TSNEPlot(seu,label=F,group="Stage",pt.size=1)+
  ggtitle(NULL)+
  scale_color_manual(values = col.spectral(length(tmp.levels.2)))
p
ggsave(filename = myFileName(prefix = "res/fig/Fig1_tsne",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)




p1 <- FeaturePlot(seu,features = c("Sox2","Pou5f1","Nanog","Gata6","Zscan4c","Bmp4")) & 
  scale_color_gradientn(colors = rdwhbu(100)) 
p1 | p
ggsave(filename = myFileName(prefix = "res/fig/FigS2_tsne_feature",suffix = ".jpg"),
       width = 12,height = 9,dpi = 350)


FeaturePlot(seu,features = c("Nodal"))

beforeEPI.exp["Nodal",]
# #####------8.5.3 we need to order the cells---------
# 
# beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20201019.rds")
# beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20201019.rds")
# 
# beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
# beforeEPI.exp <- myNormalize(beforeEPI.exp)
# 
# Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
#                  "4cell","8cell","16cell","earlyblast","midblast","lateblast")
# 
# 
# data.plot <- beforeEPI.exp
# tmp.select <- beforeEPI.meta %>%
#   filter(Stage!="fibroblast") %>%
#   pull(cell_id)
# 
# data.plot <- data.plot[,tmp.select]
# data.plot <- median_center(data.plot)
# pca.res <- prcomp(t(data.plot),center = T,scale. = F)
# rd.pca <- pca.res$x[,1:2] 
# 
# 
# ct <- beforeEPI.meta %>%
#   filter(Stage!="fibroblast")
# ct.label <- ct$Stage
# names(ct.label) <- ct$cell_id
# ct.label <- plyr::mapvalues(ct.label,from = Stage.level,to = 1:length(Stage.level))
# 
# 
# beforeEPI <- SingleCellExperiment(assays=List(counts=beforeEPI.exp[,tmp.select],
#                                               norm=beforeEPI.exp[,tmp.select]))
# reducedDims(beforeEPI) <- SimpleList(PCA = rd.pca)
# colData(beforeEPI)$Stage <- ct.label
# 
# 
# beforeEPI <- slingshot(beforeEPI, clusterLabels = 'Stage', reducedDim = 'PCA')
# summary(beforeEPI$slingPseudotime_1)
# 
# sim <- beforeEPI
# colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
# plotcol <- colors[cut(sim$slingPseudotime_1, breaks=100)]
# plot(reducedDims(sim)$PCA, col = plotcol, pch=16, asp = 1)
# lines(SlingshotDataSet(sim), lwd=2, col='black')





####--------3. Ligand, Receptor expression-----------
#####-------3.1 show Ligand, Receptor gene expression landscape--------
#####--------3.1.1 prepare for Ligand, Receptor database-----------
LRpairs.df <- read.delim(file = "database/Ligand-Receptor-Pairs/Mouse/Mouse-2020-Shao-LR-pairs.txt",stringsAsFactors = F)
LRpairs <- LRpairs.df$lr_pair
Lgenelist <- LRpairs.df$ligand_gene_symbol
Rgenelist <- LRpairs.df$receptor_gene_symbol 

data.merge.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20210114.rds")
data.merge.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20210114.rds")

gene_symbols <- rownames(data.merge.exp)
colnames(data.merge.meta)[1] <- "Cell"
l.remove <- setdiff(Lgenelist,gene_symbols)
r.remove <- setdiff(Rgenelist,gene_symbols)
index.remove <- c(which(Lgenelist %in% l.remove),which(Rgenelist %in% r.remove))
LRpairs <- LRpairs[-index.remove]
Lgenelist <- Lgenelist[-index.remove]
Rgenelist <- Rgenelist[-index.remove]
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
saveRDS(LRpairs,file = myFileName(prefix = "res/R/cytotalkdb_LRpairs",suffix = ".rds"))

#####------3.1.2 prepare data for plot-------
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)

data.merge.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20210114.rds")
data.merge.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20210114.rds")
colnames(data.merge.meta)[1] <- "Cell"

data.merge.exp.L <- data.merge.exp[Lgene,]
data.merge.exp.R <- data.merge.exp[Rgene,]

data.plot.L <- data.merge.exp.L %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>% 
  gather(key = Cell,gene.exp,-Gene) %>%
  full_join(y=data.merge.meta,by = "Cell")
data.plot.L$group <- "Lgene"

#### this was used as test
data.plot.L %>%
  rowid_to_column() %>%
  filter(is.na(gene.exp)) %>%
  head()

data.plot.R <- data.merge.exp.R %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>% 
  gather(key = Cell,gene.exp,-Gene) %>%
  full_join(y=data.merge.meta,by = "Cell")
data.plot.R$group <- "Rgene"

data.plot <-  data.merge.exp %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>% 
  gather(key = Cell,gene.exp,-Gene) %>%
  full_join(y=data.merge.meta,by = "Cell")
data.plot$group <- "Background"
nrow(data.plot)
data.plot <- rbind(data.plot.L,data.plot.R,data.plot)
nrow(data.plot)
cat(paste0("\"",unique(data.plot$Stage),"\""),sep = ",")

### prepare for data.plot
data.plot$Stage <- factor(data.plot$Stage,
                          levels = c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                                     "4cell","8cell","16cell","earlyblast","midblast","lateblast",
                                     "E5.25","E5.5","E6.25","E6.5","fibroblast"))

saveRDS(data.plot,file = myFileName(prefix = "res/R/data_merge_plot",suffix = ".rds"))



#####---------3.1.3 plot combine geom_boxplot----------------


#####---------3.1.3.1 load data and checkout----------
####load data
data.plot <- readRDS(file = "res/R/data_merge_plot_20210115.rds")
###checkout the group assignment
data.plot %>%
  group_by(Cell,group) %>%
  summarise(n=n()) %>%
  head()
###check normalization
summary(data.plot$gene.exp)
###check na
data.plot %>%
  rowid_to_column() %>%
  filter(is.na(gene.exp)) %>%
  head()
###
levels(data.plot$Stage)
head(data.plot)
data.plot$group <- factor(data.plot$group,levels = c("Lgene","Rgene","Background"))


#####-------3.1.3.2 raw gene.exp boxplot-------

tmp.data.plot <- data.plot
mypalette <- c("#ffc000","#00b050","#f2f2f2")
head(tmp.data.plot)
scales::show_col(mypalette)

x <- "Stage"
y <- "gene.exp"

p1 <- ggplot(data = tmp.data.plot ,aes(Stage,gene.exp))+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp=gene.exp),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp=gene.exp),
               fill = mypalette[1], 
               color=mypalette[1],
               outlier.shape = NA,
               width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp=gene.exp), 
               geom = "point", 
               fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp=gene.exp), 
               geom = "line", 
               color="black",
               aes(group=""),
               fun = mean)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp=-gene.exp),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp=-gene.exp),
               fill = mypalette[2], 
               color=mypalette[2],
               outlier.shape = NA,
               width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp=-gene.exp), 
               geom = "point", 
               fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp=-gene.exp), 
               geom = "line", 
               color="black",
               aes(group=""),
               fun = mean)+
  geom_hline(yintercept = 0,linetype="dashed")+
  scale_y_continuous(limits = c(-6,6),
                     breaks = seq(-6,6,2),
                     labels = as.character(c(-seq(-6,-2,by = 2),0,seq(2,6,by = 2))))+
  ylab("log2(RPKM+1)")+
  theme_cowplot(font_size  = 18)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))
myggsave(p = p1,prefix = "res/fig/Fig2_landscape_gene_exp",suffix = ".jpg",width = 8,height = 6,dpi=350)


####--------------3.2 gene mean expression changes------------------

####-------------3.2.1 calc mean expression-------------
####mean expression
####load data
data.plot <- readRDS(file = "res/R/data_merge_plot_20210115.rds")

#####Background
data.plot.mean.Background <- data.plot %>%
  filter(group=="Background") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp),
            gene.exp.sd=sd(gene.exp),
            gene.exp.cv=ifelse(mean(gene.exp)!=0,sd(gene.exp)/mean(gene.exp),0)) %>%
  ungroup()
#####Ligand
data.plot.mean.Ligand <- data.plot %>%
  filter(group=="Lgene") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp),
            gene.exp.sd=sd(gene.exp),
            gene.exp.cv=ifelse(mean(gene.exp)!=0,sd(gene.exp)/mean(gene.exp),0)) %>%
  ungroup()
  
#####Receptor
data.plot.mean.Receptor <- data.plot %>%
  filter(group=="Rgene") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp),
            gene.exp.sd=sd(gene.exp),
            gene.exp.cv=ifelse(mean(gene.exp)!=0,sd(gene.exp)/mean(gene.exp),0)) %>%
  ungroup()
  


data.plot.mean.Background$group <- "Background"
data.plot.mean.Ligand$group <- "Lgene"
data.plot.mean.Receptor$group <- "Rgene"

data.plot.mean <- rbind(data.plot.mean.Background,
                        data.plot.mean.Ligand,
                        data.plot.mean.Receptor)

data.plot.mean$group <- factor(data.plot.mean$group,
                               levels = c("Lgene","Rgene","Background"))


#####------3.2.2 plot-------


#####-----3.2.2.1 mean--------
tmp.data.plot <- data.plot.mean
levels(tmp.data.plot$group)
mypalette <- c("#ffc000","#00b050","#f2f2f2")

p2 <- ggplot(data = tmp.data.plot,aes(Stage,gene.exp.mean))+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp.mean=gene.exp.mean),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.mean=gene.exp.mean),
               fill = mypalette[1], color=mypalette[1],outlier.shape = NA,width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.mean=gene.exp.mean), geom = "point", fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.mean=gene.exp.mean), geom = "line", color="black",aes(group=""),fun = mean)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp.mean=-gene.exp.mean),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.mean=-gene.exp.mean),
               fill = mypalette[2], color=mypalette[2],outlier.shape = NA,width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.mean=-gene.exp.mean), geom = "point", fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.mean=-gene.exp.mean), geom = "line", color="black",aes(group=""),fun = mean)+
  geom_hline(yintercept = 0,linetype="dashed")+
  ylab("gene.exp.mean")+
  scale_y_continuous(limits = c(-6,6),
                     breaks = seq(-6,6,2),
                     labels = as.character(c(-seq(-6,-2,by = 2),0,seq(2,6,by = 2))))+
  theme_cowplot(font_size  = 18)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))


p2 <- p2 + ylab("gene.exp.mean")
myggsave(p = p2,prefix = "res/fig/Fig2_landscape_gene_mean_exp",suffix = ".jpg",width = 8,height = 6,dpi=350)




p3 <- ggplot(data = tmp.data.plot,aes(Stage,gene.exp.sd))+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp.sd=gene.exp.sd),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.sd=gene.exp.sd),
               fill = mypalette[1], color=mypalette[1],outlier.shape = NA,width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.sd=gene.exp.sd), geom = "point", fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.sd=gene.exp.sd), geom = "line", color="black",aes(group=""),fun = mean)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp.sd=-gene.exp.sd),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.sd=-gene.exp.sd),
               fill = mypalette[2], color=mypalette[2],outlier.shape = NA,width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.sd=-gene.exp.sd), geom = "point", fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.sd=-gene.exp.sd), geom = "line", color="black",aes(group=""),fun = mean)+
  geom_hline(yintercept = 0,linetype="dashed")+
  scale_y_continuous(limits = c(-5,5),
                     breaks = seq(-5,5,1),
                     labels = as.character(c(-seq(-5,-1,by = 1),0,seq(1,5,by = 1))))+
  ylab("gene.exp.sd")+
  theme_cowplot(font_size  = 18)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))

myggsave(p = p3,prefix = "res/fig/FigS3_landscape_gene_sd",suffix = ".jpg",width = 8,height = 6,dpi=350)




p4 <- ggplot(data = tmp.data.plot,aes(Stage,gene.exp.cv))+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp.cv=gene.exp.cv),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.cv=gene.exp.cv),
               fill = mypalette[1], color=mypalette[1],outlier.shape = NA,width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.cv=gene.exp.cv), geom = "point", fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.cv=gene.exp.cv), geom = "line", color="black",aes(group=""),fun = mean)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Background") %>%
                 mutate(gene.exp.cv=-gene.exp.cv),
               fill = mypalette[3], 
               color=mypalette[3],
               outlier.shape = NA,
               coef=0,
               width = 0.8)+
  geom_boxplot(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.cv=-gene.exp.cv),
               fill = mypalette[2], color=mypalette[2],outlier.shape = NA,width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.cv=-gene.exp.cv), geom = "point", fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(gene.exp.cv=-gene.exp.cv), geom = "line", color="black",aes(group=""),fun = mean)+
  geom_hline(yintercept = 0,linetype="dashed")+
  scale_y_continuous(limits = c(-30,30),
                     breaks = seq(-30,30,10),
                     labels = as.character(c(-seq(-30,-10,by = 10),0,seq(10,30,by = 10))))+
  ylab("gene.exp.cv")+
  theme_cowplot(font_size  = 18)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))
#p3 
myggsave(p = p4,prefix = "res/fig/FigS3_gene_cv",suffix = ".jpg",width = 8,height = 6,dpi=300)

#####---------3.2.2.2 rank---------------

LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)

data.plot.rank.Background <- data.plot.mean.Background %>%
  group_by(Stage) %>%
  mutate(rank=row_number(gene.exp.mean)) %>%
  ungroup() %>%
  dplyr::select(-group,-gene.exp.mean,-gene.exp.sd,-gene.exp.cv) %>%
  spread(key = "Stage",value = "rank") %>%
  column_to_rownames("Gene")

data.plot.rank.L <- data.plot.rank.Background[Lgene,] %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = "rank",-Gene) %>%
  mutate(group="Lgene")

data.plot.rank.R <- data.plot.rank.Background[Rgene,] %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = "rank",-Gene) %>%
  mutate(group="Rgene")

data.plot.rank <- data.plot.rank.Background %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = "rank",-Gene) %>%
  mutate(group="Background")

data.plot.rank <- rbind(data.plot.rank,data.plot.rank.L,data.plot.rank.R)

data.plot.rank$Stage <- factor(data.plot.rank$Stage,
                               levels = c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                                          "4cell","8cell","16cell","earlyblast","midblast","lateblast",
                                          "E5.25","E5.5","E6.25","E6.5","fibroblast"))

data.plot.rank$group <- factor(data.plot.rank$group,
                               levels = c("Lgene","Rgene","Background"))



tmp.data.plot <- data.plot.rank
head(tmp.data.plot)
mypalette <- c("#ffc000","#00b050","#f2f2f2")



tmp.data.plot.bar <- subset(tmp.data.plot,group == "Background") %>%
  mutate(rank=rank) %>%
  group_by(Stage) %>%
  summarise(rank=n()) %>%
  ungroup()


p5 <- ggplot(data = tmp.data.plot,aes(Stage,rank))+
  # geom_bar(data = tmp.data.plot.bar,
  #          mapping = aes(x = Stage,y = rank),fill = mypalette[3], color=mypalette[3],
  #          stat = "identity",width = 0.8)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(rank=rank), geom = "point",shape = 15, fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(rank=rank), geom = "line", color=mypalette[1],aes(group=""),fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(rank=rank), geom = "point",shape = 15, fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
                 mutate(rank=rank), geom = "line", color=mypalette[2],aes(group=""),fun = mean)+
  #geom_hline(yintercept = 0,linetype="dashed")+
  scale_y_continuous(breaks = seq(6500,8500,500))+
  ylab("rank")+
  theme_cowplot(font_size  = 18)+
  theme(legend.position = "right",
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))
p5
myggsave(p = p5,prefix = "res/fig/Fig2_gene_exp.mean_rank",suffix = ".jpg",width = 8,height = 6,dpi=350)

tmp.stage <- c("MII oocyte","zygote",
               "early 2-cell","mid 2-cell","late 2-cell",
               "4-cell","8-cell","16-cell",
               "early blastocyst","mid blastocyst","late blastocyst",
               "E5.25","E5.5","E6.25","E6.5","fibroblast")
plot.list <- plot_grid(p2+scale_x_discrete(labels=tmp.stage),
                       p5+scale_x_discrete(labels=tmp.stage),
                       ncol = 2,align = "hv")
myggsave(p = plot.list,prefix = "res/fig/Fig2_gene_landscape_combine",suffix = ".jpg",width = 16,height = 6,dpi=350)


#####--------3.3 heatmap of ligand, receptor ----------

#####--------3.3.1 pheatmap---------
####load data
data.plot <- readRDS(file = "res/R/data_merge_plot_20210115.rds")


LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)

#####Ligand
data.plot.mean.Ligand <- data.plot %>%
  filter(group=="Lgene") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp),
            gene.exp.sd=sd(gene.exp),
            gene.exp.cv=ifelse(mean(gene.exp)!=0,sd(gene.exp)/mean(gene.exp),0)) %>%
  ungroup()

#####Receptor
data.plot.mean.Receptor <- data.plot %>%
  filter(group=="Rgene") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp),
            gene.exp.sd=sd(gene.exp),
            gene.exp.cv=ifelse(mean(gene.exp)!=0,sd(gene.exp)/mean(gene.exp),0)) %>%
  ungroup()

data.plot.mean.Ligand.mat <- data.plot.mean.Ligand %>%
  dplyr::select(Gene,Stage,gene.exp.mean) %>%
  spread(key = "Stage",value = "gene.exp.mean") %>%
  column_to_rownames("Gene")

data.plot.mean.Receptor.mat <- data.plot.mean.Receptor %>%
  dplyr::select(Gene,Stage,gene.exp.mean) %>%
  spread(key = "Stage",value = "gene.exp.mean") %>%
  column_to_rownames("Gene")
rdbu <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))



num_clusters <- 6 
tmp.data.plot <- pheatmap:::scale_rows(data.plot.mean.Ligand.mat)
tmp.data.plot <- myRemoveNA(tmp.data.plot)
range(tmp.data.plot)
### filter 3
scale.max <- 3
scale.min <- -3
tmp.data.plot[tmp.data.plot > scale.max] = scale.max
tmp.data.plot[tmp.data.plot < scale.min] = scale.min

ph <- pheatmap(tmp.data.plot,
               color = rdbu(100),
               cluster_cols = F,cutree_rows = num_clusters ,
               scale = "none",
               show_rownames = F)

annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row,num_clusters)))

pheatmap::pheatmap(tmp.data.plot,
                   scale = "none", 
                   fontsize = 16,
                   labels_col = tmp.stage,
                   angle_col = 45,
                   color = rdbu(100),
                   cluster_rows = T,
                   cluster_cols = F,
                   show_rownames = F,
                   show_colnames = T,
                   cutree_rows = num_clusters,
                   main = paste0(length(Lgene)," Ligand gene"),
                   annotation_row = annotation_row)

annotation_row <- data.frame(Cluster = cutree(ph$tree_row,num_clusters))
annotation_row$Cluster <- plyr::mapvalues(annotation_row$Cluster,from = 1:6,to = c(5,6,2,4,3,1))
annotation_row$Cluster <- factor(annotation_row$Cluster)



tmp.stage <- c("MII oocyte","zygote",
               "early 2-cell","mid 2-cell","late 2-cell",
               "4-cell","8-cell","16-cell",
               "early blastocyst","mid blastocyst","late blastocyst",
               "E5.25","E5.5","E6.25","E6.5","fibroblast")
jpeg(filename = myFileName(prefix = "res/fig/Fig2_gene.mean.exp_L_heatmap",suffix = ".jpg"),
     width = 8,height = 8,res=350,units = "in")
pheatmap::pheatmap(tmp.data.plot,
                   scale = "none", 
                   fontsize = 16,
                   labels_col = tmp.stage,
                   angle_col = 45,
                   color = rdbu(100),
                   cluster_rows = T,
                   cluster_cols = F,
                   show_rownames = F,
                   show_colnames = T,
                   cutree_rows = num_clusters,
                   main = paste0(length(Lgene)," Ligand gene"),
                   annotation_row = annotation_row)
dev.off()



num_clusters <- 6
tmp.data.plot <- pheatmap:::scale_rows(data.plot.mean.Receptor.mat)
tmp.data.plot <- myRemoveNA(tmp.data.plot)
range(tmp.data.plot)
### filter 3
scale.max <- 3
scale.min <- -3
tmp.data.plot[tmp.data.plot > scale.max] = scale.max
tmp.data.plot[tmp.data.plot < scale.min] = scale.min

max(tmp.data.plot)
ph <- pheatmap(tmp.data.plot,color = rdbu(100), 
               cluster_cols = F,cutree_rows = num_clusters ,
               scale = "none",
               show_rownames = F,silent = T)

annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                     num_clusters)))
pheatmap::pheatmap(tmp.data.plot,
                   scale = "none", 
                   fontsize = 16,
                   labels_col = tmp.stage,
                   angle_col = 45, 
                   color = rdbu(100),
                   cluster_rows = T,
                   cluster_cols = F,
                   show_rownames = F,
                   show_colnames = T,
                   cutree_rows = num_clusters,
                   main = paste0(length(Rgene)," Receptor gene"),
                   annotation_row = annotation_row)

annotation_row <- data.frame(Cluster = cutree(ph$tree_row, 
                                              num_clusters))
annotation_row$Cluster <- plyr::mapvalues(annotation_row$Cluster,
                                          from = 1:6 ,
                                          to = c(4,1,2,6,5,3))
annotation_row$Cluster <- factor(annotation_row$Cluster)



jpeg(filename = myFileName(prefix = "res/fig/Fig2_gene.mean.exp_R_heatmap",suffix = ".jpg"),
     width = 8,height = 8,res=350,units = "in")
pheatmap::pheatmap(tmp.data.plot,
                   scale = "none", 
                   fontsize = 16,
                   labels_col = tmp.stage,
                   angle_col = 45, 
                   color = rdbu(100),
                   cluster_rows = T,
                   cluster_cols = F,
                   show_rownames = F,
                   show_colnames = T,
                   cutree_rows = num_clusters,
                   main = paste0(length(Rgene)," Receptor gene"),
                   annotation_row = annotation_row)
dev.off()


####----3.3.2 GO----


####using mysidb C5 for enrichment
Mm_msigdb <- msigdbr(species = "Mus musculus")
Mm_C5 <-  Mm_msigdb %>%
  dplyr::filter(gs_cat == "C5" & gs_subcat == "GO:BP") %>%
  dplyr::select(gs_name,gs_exact_source,entrez_gene) 
Mm_C5_GO <- Mm_C5[,1:2] %>%
  unique()


####Lgene
num_clusters <- 6
tmp.data.plot <- pheatmap:::scale_rows(data.plot.mean.Ligand.mat)
tmp.data.plot <- myRemoveNA(tmp.data.plot)
scale.max <- 3
scale.min <- -3
tmp.data.plot[tmp.data.plot > scale.max] = scale.max
tmp.data.plot[tmp.data.plot < scale.min] = scale.min

ph <- pheatmap(tmp.data.plot,
               color = rdbu(100),
               cluster_cols = F,cutree_rows = num_clusters ,
               scale = "none",
               show_rownames = F)

annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                     num_clusters)))
### perform GO
for (ii in 1:num_clusters) {
  test.gene <- annotation_row %>%
    rownames_to_column("Gene") %>%
    filter(Cluster==ii) %>%
    pull(Gene)
  
  gene.list <- test.gene
  eg = bitr(gene.list, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)
  em <- enricher(eg$ENTREZID, TERM2GENE=Mm_C5[,c(1,3)],pvalueCutoff = -1,qvalueCutoff = -1)
  em <- setReadable(em, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  tmp.res.df <- em@result %>%
    dplyr::filter(p.adjust < 0.05 ) %>%
    arrange(p.adjust) %>%
    rownames_to_column("rownames.id")
  em.res <- merge(tmp.res.df,Mm_C5_GO,by.x="ID",by.y="gs_name")
  em.res <- em.res %>%
    arrange(p.adjust)
  #em.res$Description
  write.table(em.res,file = myFileName(prefix = paste0("res/txt/","Ligand_Cluster_",ii,"_msigdb_GO"),suffix = ".txt"),quote = F,sep = "\t",row.names = F)
  p <- myenrichr_GOplot(df = em.res,
                        fill.color = "#C6DBEFFF",
                        show_number = 20,
                        title = paste0("Ligand gene ","cluster ",ii," ",length(gene.list)," gene"),
                        font.size = 18,
                        p.adjust.cut = 0.05,
                        plot.ylab = NULL,
                        term.pos.adjust = 0)+
    theme(axis.line.y = element_blank())
  ggsave(filename = myFileName(prefix = paste0("res/fig/","Ligand_Cluster_",ii,"_msigdb_GO"),suffix = ".jpg"),width = 8,height = 8)
  
}

####RGene
num_clusters <- 6
tmp.data.plot <- pheatmap:::scale_rows(data.plot.mean.Receptor.mat)
tmp.data.plot <- myRemoveNA(tmp.data.plot)
scale.max <- 3
scale.min <- -3
tmp.data.plot[tmp.data.plot > scale.max] = scale.max
tmp.data.plot[tmp.data.plot < scale.min] = scale.min

ph <- pheatmap(tmp.data.plot,
               color = rdbu(100),
               cluster_cols = F,cutree_rows = num_clusters ,
               scale = "none",
               show_rownames = F)

annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                     num_clusters)))
test.gene <- annotation_row %>%
  rownames_to_column("Gene") %>%
  filter(Cluster==5) %>%
  pull(Gene)


### perform GO
for (ii in 1:num_clusters) {
  test.gene <- annotation_row %>%
    rownames_to_column("Gene") %>%
    filter(Cluster==ii) %>%
    pull(Gene)
  
  gene.list <- test.gene
  eg = bitr(gene.list, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)
  em <- enricher(eg$ENTREZID, TERM2GENE=Mm_C5[,c(1,3)],pvalueCutoff = -1,qvalueCutoff = -1)
  em <- setReadable(em, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
  tmp.res.df <- em@result %>%
    dplyr::filter(p.adjust < 0.05 ) %>%
    arrange(p.adjust) %>%
    rownames_to_column("rownames.id")
  em.res <- merge(tmp.res.df,Mm_C5_GO,by.x="ID",by.y="gs_name")
  em.res <- em.res %>%
    arrange(p.adjust)
  #em.res$Description
  write.table(em.res,file = myFileName(prefix = paste0("res/txt/","Receptor_Cluster_",ii,"_msigdb_GO"),suffix = ".txt"),quote = F,sep = "\t",row.names = F)
  p <- myenrichr_GOplot(df = em.res,
                        fill.color = "#C6DBEFFF",
                        show_number = 20,
                        title = paste0("Receptor gene ","cluster ",ii," ",length(gene.list)," gene"),
                        font.size = 18,
                        p.adjust.cut = 0.05,
                        plot.ylab = NULL,
                        term.pos.adjust = 0)+
    theme(axis.line.y = element_blank())
  ggsave(filename = myFileName(prefix = paste0("res/fig/","Receptor_Cluster_",ii,"_msigdb_GO"),suffix = ".jpg"),width = 8,height = 8)
  
}



###--------4. eLR anlysis---------

#####---------4.1 Get eLR, filter in before implantation------------

#####---------4.1.1 load exp data----------
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")

beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)

Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
           "4cell","8cell","16cell","earlyblast","midblast",  
           "lateblast","fibroblast")

head(beforeEPI.meta)

beforeEPI.exp.mean <- beforeEPI.exp %>%
  rownames_to_column("Gene") %>%
  gather(key = "cell_id",value = "gene.exp",-Gene) %>%
  left_join(beforeEPI.meta,"cell_id") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp)) %>%
  ungroup() %>%
  mutate(Stage=factor(Stage,levels = Stage.level)) %>%
  spread(key="Stage",value = "gene.exp.mean") %>%
  column_to_rownames("Gene")

boxplot(beforeEPI.exp.mean)
  

#####---------4.1.2 load and filter LR pairs---------------
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210114.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)

####test weather there exist Lgene.list == Rgene.list
sum(Lgene.list==Rgene.list)

####none, that's good!

#####-------------------4.1.3 calc cell stage correlation----------------------------------

LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210114.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
           "4cell","8cell","16cell","earlyblast","midblast","lateblast")

tmp.cell_id <- beforeEPI.meta %>%
  filter(Stage %in% Stage.level) %>%
  pull(cell_id)

tmp.L.exp <- beforeEPI.exp[Lgene.list,tmp.cell_id] 
tmp.R.exp <- beforeEPI.exp[Rgene.list,tmp.cell_id]

tmp.res <- sapply(1:nrow(tmp.L.exp), function(i) cor(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,])))
tmp.res.rate.L <- sapply(1:nrow(tmp.L.exp), function(i) sum(as.numeric(tmp.L.exp[i,])>0)/length(as.numeric(tmp.L.exp[i,])>0))
tmp.res.rate.R <- sapply(1:nrow(tmp.R.exp), function(i) sum(as.numeric(tmp.R.exp[i,])>0)/length(as.numeric(tmp.R.exp[i,])>0))
tmp.res.cor.1 <- data.frame(Lgene=Lgene.list,
                            Rgene=Rgene.list,
                            LRpairs=paste0(Lgene.list,"_",Rgene.list),
                            PCC=tmp.res,
                            detection.rate.L=tmp.res.rate.L,
                            detection.rate.R=tmp.res.rate.R,
                            stringsAsFactors = F) %>%
  dplyr::filter(detection.rate.L > 0.05 & detection.rate.R > 0.05) %>%
  arrange(-PCC)

cutoff <- 0.2
tmp.res.LR.pos <- tmp.res.cor.1 %>%
  filter(PCC > cutoff)
tmp.res.LR.neg <- tmp.res.cor.1 %>%
  filter(PCC < -cutoff)

tmp.res.LR <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff)

saveRDS(tmp.res.cor.1,file = "res/R/putative_eLR_correlation_20210129.Rds")

#####-----------4.1.4 visualization-----------

hist(tmp.res.cor.1$PCC)

cor.test(rnorm(100),rnorm(100))

test.color <- colorRampPalette(RColorBrewer::brewer.pal(name = "OrRd",n=9))
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))
Stage.level <-  c("Oocyte","zygote","early2cell",
                  "mid2cell","late2cell","4cell",
                  "8cell","16cell","earlyblast",
                  "midblast","lateblast")
paper.color <- c("#FCB03B","#0F6937","#67BD45",
                 "#169192","#B09DBA","#954496",
                 "#4E1550","#8EAF3D","#F05153",
                 "#5358A5","#79A7AA")
scales::show_col(paper.color)
myCorplot <- function(mat=beforeEPI.exp[,tmp.cell_id],
                      mat.meta=beforeEPI.meta,
                      tmp.Lgene=tmp.res.cor.1$Lgene[1],
                      tmp.Rgene=tmp.res.cor.1$Rgene[1],
                      PCC=tmp.res.cor.1$PCC[1],
                      Stage.level=Stage.level,
                      color=test.color(length(Stage.level))){
  
  mat.L <- mat[tmp.Lgene,]
  mat.R <- mat[tmp.Rgene,]
  
  res.df <- data.frame(t(mat.L),
                       t(mat.R),
                       stringsAsFactors = F) %>%
    rownames_to_column("cell_id") %>%
    left_join(mat.meta,by = "cell_id") %>%
    mutate(Stage=factor(Stage,levels = Stage.level))
  
  p<- ggplot(data = res.df,
             aes_string(x=tmp.Lgene,y=tmp.Rgene,color="Stage"))+
    geom_point(size=8)+
    scale_color_manual(values = color)+
    ggtitle(paste0(tmp.Lgene,"_",tmp.Rgene,",cor:",PCC))+
    theme_cowplot(font_size = 15)+
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text( size = 15 ) ,
          axis.text = element_text(size = 15),
          axis.line.x = element_line( size = 1 ),
          axis.line.y = element_line( size = 1 ))
  return(p)
}
tmp.idx <- nrow(tmp.res.LR.pos)
p <- myCorplot(mat=beforeEPI.exp[,tmp.cell_id],
          mat.meta=beforeEPI.meta,
          tmp.Lgene=tmp.res.LR.pos$Lgene[1],
          tmp.Rgene=tmp.res.LR.pos$Rgene[1],
          PCC=tmp.res.LR.pos$PCC[1],
          Stage.level=Stage.level,
          color=test.color.2(length(Stage.level)))
p
myggsave(p = p,prefix = paste0("res/fig/eLR_",tmp.res.LR.pos$LRpairs[1]),suffix = ".jpg",width = 8,height = 8,dpi=350)

tmp.idx <- nrow(tmp.res.LR.neg)
p <- myCorplot(mat=beforeEPI.exp[,tmp.cell_id],
          mat.meta=beforeEPI.meta,
          tmp.Lgene=tmp.res.LR.neg$Lgene[tmp.idx],
          tmp.Rgene=tmp.res.LR.neg$Rgene[tmp.idx],
          PCC=tmp.res.LR.neg$PCC[tmp.idx],
          Stage.level=Stage.level,
          color=test.color.2(length(Stage.level)))
p
myggsave(p = p,prefix = paste0("res/fig/eLR_",tmp.res.LR.neg$LRpairs[tmp.idx]),suffix = ".jpg",width = 8,height = 8,dpi=350)


####--------4.2 Stage specific LR ----------

####--------4.2.1 eLR function----------

####--------4.2.1.1 load data and calculation strength----------

####load LRpairs 1744
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20201214.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)
####load mean mat and calculate strength
data.mean.plot.mat <- readRDS(file = "res/R/data.plot.mean.mat_20201214.rds")

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")
tmp.mat.L <- data.mean.plot.mat[Lgene.list,] 
tmp.mat.R <- data.mean.plot.mat[Rgene.list,]
tmp.res <- (tmp.mat.L * tmp.mat.R) %>%
  rownames_to_column("LRid") %>%
  mutate(LRid = LRpairs) %>%
  column_to_rownames("LRid")
saveRDS(tmp.res,file = myFileName(prefix = "res/R/early.scRNAseq_cytotalk_LR_IS",suffix = '.rds'))

####-------4.2.1.2 pheatmap--------------

####load eLR
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
tmp.res.LR.pos <- tmp.res.cor.1 %>%
  filter(PCC > cutoff)
tmp.res.LR.neg <- tmp.res.cor.1 %>%
  filter(PCC < -cutoff)
tmp.res.LR <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff)


####load mat
IS.mat <- readRDS(file = "res/R/early.scRNAseq_cytotalk_LR_IS_20201215.rds")

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")
####LR.pos
tmp.mat.plot <- IS.mat[tmp.res.LR.pos$LRpairs,Stage.levels]
tmp.mat.plot <- log10(tmp.mat.plot+1)
tmp.mat.plot <- myRemoveNA(tmp.mat.plot)
num_clusters <- 3
ph <- pheatmap(tmp.mat.plot,
               color = rdbu(100), 
               cluster_cols = F,
               cutree_rows = num_clusters ,
               scale = "none",
               show_rownames = F,
               silent = T)
annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                     num_clusters)))


p_eLR <- annotation_row %>%
  rownames_to_column("LRpairs")
saveRDS(object = p_eLR,file = myFileName(prefix = "res/R/p_eLR",suffix = ".rds"))

get_plot_dims <- function(heat_map){
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}

tmp.stage <- c("MII oocyte","zygote",
                            "early 2-cell","mid 2-cell","late 2-cell",
                            "4-cell","8-cell","16-cell",
                            "early blastocyst","mid blastocyst","late blastocyst",
                            "E5.25","E5.5","E6.25","E6.5","fibroblast")
ph <- pheatmap::pheatmap(tmp.mat.plot,
                         #scale = "row", 
                         color = rdbu(100),
                         border_color = NA,
                         cluster_rows = T,
                         cellwidth = 16,
                         cellheight = 4,
                         cluster_cols = F,
                         show_rownames = F,
                         show_colnames = T,
                         cutree_rows = num_clusters,
                         fontsize = 15,
                         main = paste0("LR pos ",cutoff," ",nrow(tmp.mat.plot)),
                         silent =T,
                         annotation_row = annotation_row,
                         labels_col = tmp.stage,
                         angle_col = 45)

plot.dims <- get_plot_dims(ph)

jpeg(file = myFileName(prefix = paste0("res/fig/eLR_pos_",cutoff,"_",nrow(tmp.mat.plot)),suffix = ".jpg"),
     width = plot.dims$width+0.5,height = plot.dims$height+0.5,res = 350,units = "in")
ph
dev.off()





####LR.neg
tmp.mat.plot <- IS.mat[tmp.res.LR.neg$LRpairs,Stage.levels]
tmp.mat.plot <- log10(tmp.mat.plot+1)
tmp.mat.plot <- myRemoveNA(tmp.mat.plot)

num_clusters <- 3
ph <- pheatmap(tmp.mat.plot,
               color = rdbu(100), 
               cluster_cols = F,
               cutree_rows = num_clusters ,
               #scale = "row",
               show_rownames = F,silent = T)

tmp.cluster <- cutree(ph$tree_row, 
                      num_clusters)
tmp.cluster <- plyr::mapvalues(tmp.cluster,from = c(1,2,3),to = c(3,2,1))
annotation_row <- data.frame(Cluster = factor(tmp.cluster))


n_eLR <- annotation_row %>% 
  rownames_to_column("LRpairs")
saveRDS(object = n_eLR,file = myFileName(prefix = "res/R/n_eLR",suffix = ".rds"))


ph <- pheatmap::pheatmap(tmp.mat.plot,
                   #scale = "row", 
                   color = rdbu(100),
                   cluster_rows = T,
                   cluster_cols = F,
                   show_rownames = F,
                   show_colnames = T,
                   border_color = NA,
                   cellwidth = 16,
                   cellheight = 4, 
                   fontsize = 15,
                   cutree_rows = num_clusters,
                   main = paste0("LR neg ",cutoff," ",nrow(tmp.mat.plot)),
                   annotation_row = annotation_row,
                   silent = T,
                   labels_col = tmp.stage,
                   angle_col = 45)

plot.dims <- get_plot_dims(ph)


jpeg(file = myFileName(prefix = paste0("res/fig/eLR_neg_",cutoff,"_",nrow(tmp.mat.plot)),suffix = ".jpg"),
     width = plot.dims$width+0.5,height = plot.dims$height+0.5,units = "in",res = 300)
ph
dev.off()

#####-----GO----------

####LR.pos
tmp.mat.plot <- IS.mat[tmp.res.LR.pos$LRpairs,Stage.levels]
tmp.mat.plot <- log10(tmp.mat.plot+1)
tmp.mat.plot <- myRemoveNA(tmp.mat.plot)
num_clusters <- 3
ph <- pheatmap(tmp.mat.plot,
               color = rdbu(100), 
               cluster_cols = F,
               cutree_rows = num_clusters ,
               #scale = "row",
               show_rownames = F,silent = T)
annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                     num_clusters)))

for(ii in 1:num_clusters){
  tmp_eLR <- annotation_row %>%
       rownames_to_column("Gene") %>%
       filter(Cluster==ii) %>%
       pull(Gene)
  eLR_L <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[1])) 
  eLR_R <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[2]))
  eLR_gene <- c(eLR_L,eLR_R)
  myPerformDavidGO(gene.lits = eLR_gene,
                     output.table = myFileName(prefix = paste0("res/txt/eLR_pos_",cutoff,"_","cluster",ii),suffix = ".txt"),
                     output.fig = myFileName(prefix = paste0("res/fig/eLR_pos_",cutoff,"_","cluster",ii),suffix = ".jpg"),
                     plot.title = paste0("eLR_pos_",'cluster',ii),
                     annotation = "GOTERM_BP_DIRECT")

}

####LR.neg
tmp.mat.plot <- IS.mat[tmp.res.LR.neg$LRpairs,Stage.levels]
tmp.mat.plot <- log10(tmp.mat.plot+1)
tmp.mat.plot <- myRemoveNA(tmp.mat.plot)

num_clusters <- 3
ph <- pheatmap(tmp.mat.plot,
               color = rdbu(100), 
               cluster_cols = F,
               cutree_rows = num_clusters ,
               #scale = "row",
               show_rownames = F,silent = T)
annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, 
                                                     num_clusters)))

for(ii in 1:num_clusters){
  tmp_eLR <- annotation_row %>%
    rownames_to_column("Gene") %>%
    filter(Cluster==ii) %>%
    pull(Gene)
  eLR_L <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[1])) 
  eLR_R <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[2]))
  eLR_gene <- c(eLR_L,eLR_R)
  myPerformDavidGO(gene.lits = eLR_gene,
                   output.table = myFileName(prefix = paste0("res/txt/eLR_neg_",cutoff,"_","cluster",ii),suffix = ".txt"),
                   output.fig = myFileName(prefix = paste0("res/fig/eLR_neg_",cutoff,"_","cluster",ii),suffix = ".jpg"),
                   plot.title = paste0("eLR_neg_",'cluster',ii),
                   annotation = "GOTERM_BP_DIRECT")
  
}

####-------update GO plot------
tmp.df <- read.delim(file = "res/txt/eLR_pos_0.2_cluster1_20201215.txt")

p1 <- myDavid_GOplot(df = tmp.df,
                    fill.color = "#C6DBEFFF",
                    show_number = 30,
                    title = "p-eLR cluster 1",
                    font.size = 18,
                    term.pos.adjust = 0)+
  theme(axis.line.y = element_blank())

myggsave(p = p1,prefix = "res/fig/eLR_pos_cluster1_GO",suffix = ".jpg",width = 8,height = 8,dpi=350)


tmp.df <- read.delim(file = "res/txt/eLR_pos_0.2_cluster2_20201215.txt")

p1 <- myDavid_GOplot(df = tmp.df,
                     fill.color = "#C6DBEFFF",
                     show_number = 30,
                     title = "p-eLR cluster 2",
                     font.size = 18,
                     term.pos.adjust = 0)+
  theme(axis.line.y = element_blank())

myggsave(p = p1,prefix = "res/fig/eLR_pos_cluster2_GO",suffix = ".jpg",width = 8,height = 8,dpi=350)

tmp.df <- read.delim(file = "res/txt/eLR_pos_0.2_cluster3_20201215.txt")
p1 <- myDavid_GOplot(df = tmp.df,
                     fill.color = "#C6DBEFFF",
                     show_number = 30,
                     title = "p-eLR cluster 3",
                     font.size = 18,
                     term.pos.adjust = 0)+
  theme(axis.line.y = element_blank())

myggsave(p = p1,prefix = "res/fig/eLR_pos_cluster3_GO",suffix = ".jpg",width = 8,height = 8,dpi=350)


tmp.df <- read.delim(file = "res/txt/eLR_neg_0.2_cluster1_20201215.txt")
p2 <- myDavid_GOplot(df = tmp.df,
                     fill.color = "#C6DBEFFF",
                     show_number = 30,
                     title = "n-eLR cluster 1",
                     font.size = 18,
                     term.pos.adjust = 0)+
  theme(axis.line.y = element_blank())
p2
myggsave(p = p2,prefix = "res/fig/eLR_neg_cluster1_GO",suffix = ".jpg",width = 8,height = 6,dpi=350)



tmp.df <- read.delim(file = "res/txt/eLR_neg_0.2_cluster2_20201215.txt")
p2 <- myDavid_GOplot(df = tmp.df,
                     fill.color = "#C6DBEFFF",
                     show_number = 30,
                     title = "n-eLR cluster 2",
                     font.size = 18,
                     term.pos.adjust = 0)+
  theme(axis.line.y = element_blank())
p2
myggsave(p = p2,prefix = "res/fig/eLR_neg_cluster2_GO",suffix = ".jpg",width = 8,height = 8,dpi=350)



tmp.df <- read.delim(file = "res/txt/eLR_neg_0.2_cluster3_20201215.txt")

p2 <- myDavid_GOplot(df = tmp.df,
                    fill.color = "#C6DBEFFF",
                    show_number = 30,
                    title = "n-eLR cluster 3",
                    font.size = 18,
                    term.pos.adjust = 0)+
  theme(axis.line.y = element_blank())
p2
myggsave(p = p2,prefix = "res/fig/eLR_neg_cluster3_GO",suffix = ".jpg",width = 8,height = 4,dpi=350)


#####-----4.2.2 count Stage specific eLR --------

####load IS mat
Stage.levels <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
           "4cell","8cell","16cell","earlyblast","midblast",  
           "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")

IS.mat <- readRDS(file = "res/R/early.scRNAseq_cytotalk_LR_IS_20201215.rds")
####load eLR
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
tmp.res.LR.pos <- tmp.res.cor.1 %>%
  filter(PCC > cutoff)
tmp.res.LR.neg <- tmp.res.cor.1 %>%
  filter(PCC < -cutoff)
tmp.res.LR <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff)

candidate_eLR <- IS.mat[tmp.res.LR$LRpairs,Stage.levels]
head(candidate_eLR)
candidate_eLR <- myRemoveNA(candidate_eLR)


test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))
res <- apply(candidate_eLR,2,function(x) sum(x!=0))

data.plot <- as.data.frame(res) %>%
  rownames_to_column("Stage") %>%
  mutate(Stage = factor(Stage,levels = Stage.levels))

tmp.stage <- c("MII oocyte","zygote",
               "early 2-cell","mid 2-cell","late 2-cell",
               "4-cell","8-cell","16-cell",
               "early blastocyst","mid blastocyst","late blastocyst",
               "E5.25","E5.5","E6.25","E6.5","fibroblast")

p <- ggplot(data = data.plot,aes(Stage,res,fill=Stage))+
  geom_bar(stat = "identity")+
  geom_text(aes(label = res),vjust = -0.2,size = 5)+
  ylab("count")+
  scale_x_discrete(label=tmp.stage)+
  scale_fill_manual(values = test.color.2(length(Stage.levels)))+
  theme_cowplot(font_size = 18)+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.text.x = element_text(angle = 30,hjust = 0.8),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ))
p
ggsave(p,filename =  myFileName(prefix = "res/fig/eLR_stage_specific_",suffix = ".jpg"),
       width = 8,height = 8)
  

eLR_putative <- rownames(candidate_eLR)
saveRDS(eLR_putative,file = "res/R/eLR_putative_20201215.rds")

####--------4.2.3 intersetct-----
Stage.levels <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                  "4cell","8cell","16cell","earlyblast","midblast",  
                  "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")

eLR_stage_list.df <- candidate_eLR %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "Stage",value = "strength",-LRpairs) %>%
  filter(strength != 0)

eLR_stage_list <- eLR_stage_list.df %>%
  dplyr::select(LRpairs,Stage)

eLR_stage_list <- sapply(Stage.levels,
                         FUN = function(x){
                           eLR_stage_list %>%
                             dplyr::filter(Stage == x) %>%
                             pull(LRpairs)
                           },USE.NAMES = T)



#### please not run
for(ii in 1:length(Stage.levels)){
  tmp.index <- ii:(ii+4)
  tmp.subset <- Stage.levels[tmp.index]
  tmp.plot.list <- sapply(tmp.subset,function(x) eLR_stage_list[[x]],USE.NAMES = T)
  venn.diagram(tmp.plot.list,
               filename = myFileName(prefix = paste0("res/fig/venn","_",
                                                     Stage.levels[ii],"_",Stage.levels[ii+4]),
                                     suffix = ".jpg"),
               height = 8,
               width = 8,
               units = "in",
               col = "transparent",
               fill = c("#fec04e","#cba6d1","#74a5cc",
                        "#84c683","#f18d8c"),
               cat.dist = rep(0.08,length(tmp.subset)),
               cex = 2,
               cat.cex = 2,
               resolution = 350)
}
?venn.diagram
tmp.subset <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                 "4cell","8cell","16cell","earlyblast","midblast",  
                 "lateblast","E5.25","E5.5","E6.25","E6.5")
tmp.plot.list <- sapply(tmp.subset,function(x) eLR_stage_list[[x]],USE.NAMES = T)
all.LR.pairs <- Reduce(intersect,tmp.plot.list)


tmp_eLR <- all.LR.pairs
eLR_L <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[1])) 
eLR_R <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[2]))
eLR_gene <- c(eLR_L,eLR_R)

myPerformDavidGO(gene.lits = eLR_gene,
                 output.table = "res/txt/all_stage_LRpairs_intersect_20201102.txt",
                 output.fig = "res/fig/all_stage_LRpairs_intersect_20201102.jpg",
                 plot.title = "all_stage_intersect_LRpairs",
                 annotation = "GOTERM_BP_DIRECT")

#####-------4.3 gain lost--------

Stage.levels <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                  "4cell","8cell","16cell","earlyblast","midblast",  
                  "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")

eLR_stage_list.df <- candidate_eLR %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "Stage",value = "strength",-LRpairs) %>%
  filter(strength != 0)

eLR_stage_list <- eLR_stage_list.df %>%
  dplyr::select(LRpairs,Stage)

eLR_stage_list <- sapply(Stage.levels,
                         FUN = function(x){
                           eLR_stage_list %>%
                             dplyr::filter(Stage == x) %>%
                             pull(LRpairs)
                         },USE.NAMES = T)

lost <- unlist(lapply(1:(length(Stage.levels)-1),
                      FUN = function(x) length(setdiff(eLR_stage_list[[x]],
                                                       eLR_stage_list[[x+1]]))))
lost <- c(lost,0)

gain <- unlist(lapply(1:(length(Stage.levels)-1),
                      FUN = function(x) length(setdiff(eLR_stage_list[[x+1]],
                                                       eLR_stage_list[[x]]))))
gain <- c(0,gain)


data.plot <- data.frame(Stage=Stage.levels,
                        lost = lost, 
                        gain = gain,
                        stringsAsFactors = F) %>%
  gather(key = "group",value = "count",-Stage)
tmp.1 <- eLR_stage_list.df %>%
  group_by(Stage) %>%
  summarise(count = n()) %>%
  ungroup()
tmp.1$group = "expressed"
tmp.1 <- tmp.1[,c("Stage","group","count")]

data.plot <- rbind(data.plot,tmp.1)

data.plot$Stage <- factor(data.plot$Stage,levels = Stage.levels)

head(data.plot)
tmp.stage <- c("MII oocyte","zygote",
               "early 2-cell","mid 2-cell","late 2-cell",
               "4-cell","8-cell","16-cell",
               "early blastocyst","mid blastocyst","late blastocyst",
               "E5.25","E5.5","E6.25","E6.5","fibroblast")
p <- ggplot(subset(data.plot,group!="expressed"),aes(Stage,count,fill=group,label=count))+
  geom_col(position = "dodge")+
  geom_text(position = position_dodge(0.9),vjust = -0.2,size=5)+
  scale_x_discrete(label=tmp.stage)+
  theme_cowplot(font_size = 15)+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p = p,prefix = "res/fig/eLR_gain_lost",suffix = ".jpg",width = 8 ,height = 8, dpi = 350)


#####--------4.5 ZGA gene-----------

####------4.5.1 load data--------

###------zhangyi zga------

zga_gene_zhangyi <- read_delim(file = "database/other_ZGA/Nfya_KD_zhangyi.txt",delim = "\t") %>%
  pull(gene)

tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
tmp.res.LR.pos <- tmp.res.cor.1 %>%
  filter(PCC > cutoff)
tmp.res.LR.neg <- tmp.res.cor.1 %>%
  filter(PCC < -cutoff)
tmp.res.LR <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff)

tmp_eLR <- tmp.res.LR$LRpairs
eLR_L <- tmp.res.LR$Lgene
eLR_R <- tmp.res.LR$Rgene
####choose zhangyi
intersect(eLR_L,zga_gene_zhangyi)
intersect(eLR_R,zga_gene_zhangyi)


####-----4.5.2 ZGA: intersect with zhangyi ZGA gene--------

tmp.idx.1 <- which(eLR_L %in% zga_gene_zhangyi)
tmp.idx.2 <- which(eLR_R %in% zga_gene_zhangyi)
tmp.idx <- union(tmp.idx.1,tmp.idx.2)

tmp.LR.pairs <- tmp_eLR[tmp.idx]

tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
length(intersect(eLR_L,zga_gene_zhangyi))
length(intersect(eLR_R,zga_gene_zhangyi))

tmp.data.plot <- tmp.res.cor.1 %>%
  filter(LRpairs %in% tmp.LR.pairs) %>%
  arrange(-PCC)


#### using beforeEPI data
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")

beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)

Stage.level <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                  "4cell","8cell","16cell","earlyblast","midblast",  
                  "lateblast")


tmp.cell_id <- beforeEPI.meta %>%
  filter(Stage %in% Stage.level) %>%
  pull(cell_id)

myCorplot <- function(mat=beforeEPI.exp[,tmp.cell_id],
                      mat.meta=beforeEPI.meta,
                      tmp.Lgene=tmp.res.cor.1$Lgene[1],
                      tmp.Rgene=tmp.res.cor.1$Rgene[1],
                      PCC=tmp.res.cor.1$PCC[1],
                      Stage.level=Stage.level,
                      color=test.color(length(Stage.level))){
  
  mat.L <- mat[tmp.Lgene,]
  mat.R <- mat[tmp.Rgene,]
  
  res.df <- data.frame(t(mat.L),
                       t(mat.R),
                       stringsAsFactors = F) %>%
    rownames_to_column("cell_id") %>%
    left_join(mat.meta,by = "cell_id") %>%
    mutate(Stage=factor(Stage,levels = Stage.level))
  
  p<- ggplot(data = res.df,
             aes_string(x=tmp.Lgene,y=tmp.Rgene,color="Stage"))+
    geom_point(size=8)+
    scale_color_manual(values = color)+
    ggtitle(paste0(tmp.Lgene,"_",tmp.Rgene,",cor:",PCC))+
    theme_cowplot(font_size = 15)+
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text( size = 15 ) ,
          axis.text = element_text(size = 15),
          axis.line.x = element_line( size = 1 ),
          axis.line.y = element_line( size = 1 ))
  return(p)
}


p <- myCorplot(mat=beforeEPI.exp[,tmp.cell_id],
               mat.meta=beforeEPI.meta,
               tmp.Lgene=tmp.data.plot$Lgene[1],
               tmp.Rgene=tmp.data.plot$Rgene[1],
               PCC=tmp.data.plot$PCC[1],
               Stage.level=Stage.level,
               color=test.color.2(length(Stage.level)))+
  scale_x_discrete(label=tmp.stage)+
  theme_cowplot(font_size = 18)+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
      axis.line.x = element_line( size = 1 ),
      axis.line.y = element_line( size = 1 ))
p
myggsave(p = p,prefix = paste0("res/fig/eLR_",tmp.data.plot$LRpairs[1]),suffix = ".jpg",width = 8,height = 8,dpi=300)



p <- myCorplot(mat=beforeEPI.exp[,tmp.cell_id],
               mat.meta=beforeEPI.meta,
               tmp.Lgene=tmp.data.plot$Lgene[nrow(tmp.data.plot)],
               tmp.Rgene=tmp.data.plot$Rgene[nrow(tmp.data.plot)],
               PCC=tmp.data.plot$PCC[nrow(tmp.data.plot)],
               Stage.level=Stage.level,
               color=test.color.2(length(Stage.level)))+
  theme_cowplot(font_size = 18)+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ))
p
myggsave(p = p,prefix = paste0("res/fig/eLR_",tmp.data.plot$LRpairs[nrow(tmp.data.plot)]),suffix = ".jpg",width = 8,height = 8,dpi=300)

####-----4.5.3 point plot--------

data.merge.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20201028.rds")
data.merge.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20201028.rds")

Stage.level <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                  "4cell","8cell","16cell","earlyblast","midblast",  
                  "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")


tmp.cell_id <- data.merge.meta %>%
  filter(Stage %in% Stage.level) %>%
  pull(cell_id)


myScatterPlot <- function(mat = data.merge.exp[,tmp.cell_id],
                          mat.meta = data.merge.meta,
                          tmp.Lgene = tmp.data.plot$Lgene[nrow(tmp.data.plot)],
                          tmp.Rgene = tmp.data.plot$Rgene[nrow(tmp.data.plot)],
                          PCC = tmp.data.plot$PCC[nrow(tmp.data.plot)]){
  mat.L <- mat[tmp.Lgene,]
  mat.R <- mat[tmp.Rgene,]
  
  res.df <- data.frame(t(mat.L),
                       t(mat.R),
                       stringsAsFactors = F) %>%
    rownames_to_column("cell_id") %>%
    left_join(mat.meta,by = "cell_id") %>%
    mutate(Stage=factor(Stage,levels = Stage.level)) %>%
    gather(key = "gene",value = "gene.exp",-cell_id,-Stage) %>%
    mutate(gene=factor(gene,levels = c(tmp.Lgene,tmp.Rgene)))
  
  
  
  p <- ggplot(data = res.df,aes_string(x="Stage",y="gene.exp",colour="gene"))+
    geom_point(alpha=0.1,size=1)+
    geom_smooth(method = "gam",aes_string(group="gene"),se = F)+
    scale_color_manual(values = c("#E41A1C","#377EB8"))+
    ylab("gene.exp")+
    ggtitle(paste0(tmp.Lgene,"_",tmp.Rgene,",cor:",PCC))+
    theme_cowplot(font_size = 18)+
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text( size = 18 ) ,
          axis.text = element_text(size = 18),
          axis.line.x = element_line( size = 1 ),
          axis.line.y = element_line( size = 1 ),
          axis.text.x = element_text(hjust = 0.8,angle = 30))
  return(p)
}

p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.data.plot$Lgene[nrow(tmp.data.plot)],
                   tmp.Rgene = tmp.data.plot$Rgene[nrow(tmp.data.plot)],
                   PCC = tmp.data.plot$PCC[nrow(tmp.data.plot)])+
  scale_x_discrete(label=tmp.stage)

p
myggsave(p = p,prefix = paste0("res/fig/eLR_ZGA_scatter_",tmp.data.plot$LRpairs[nrow(tmp.data.plot)]),suffix = ".jpg",width = 8,height = 8,dpi=350)


p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.data.plot$Lgene[1],
                   tmp.Rgene = tmp.data.plot$Rgene[1],
                   PCC = tmp.data.plot$PCC[1])+
  scale_x_discrete(label=tmp.stage)
p
myggsave(p = p,prefix = paste0("res/fig/eLR_ZGA_scatter_",tmp.data.plot$LRpairs[1]),suffix = ".jpg",width = 8,height = 8,dpi=350)


tmp.data.plot$Lgene[nrow(tmp.data.plot)] %in% zga_gene_zhangyi
tmp.data.plot$Rgene[nrow(tmp.data.plot)] %in% zga_gene_zhangyi

tmp.data.plot$Lgene[1] %in% zga_gene_zhangyi
tmp.data.plot$Rgene[1] %in% zga_gene_zhangyi





####--------5.coordination-----------

####--------5.1 load data-----------

####--------5.1.1 single-cell--------
####load data
data.plot <- readRDS(file = "res/R/data_merge_plot_20201214.rds")
#####Background
data.plot.mean.Background <- data.plot %>%
  filter(group=="Background") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp),
            gene.exp.sd=sd(gene.exp),
            gene.exp.cv=ifelse(mean(gene.exp)!=0,sd(gene.exp)/mean(gene.exp),0)) %>%
  ungroup()
data.plot.mean.mat <- data.plot.mean.Background %>%
  dplyr::select(Gene,Stage,gene.exp.mean) %>%
  spread(key = "Stage",value = "gene.exp.mean") %>%
  column_to_rownames("Gene")
####test order
boxplot(data.plot.mean.mat)
saveRDS(data.plot.mean.mat,
        file = myFileName(prefix = "res/R/data.plot.mean.mat",suffix = ".rds"))
saveRDS(data.plot.mean.Background,file = myFileName(prefix = "res/R/data.plot.mean.Background",suffix = ".rds"))




####----5.1.2 xiewei data--------

xiewei.data.exp <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
xiewei.data.exp <- myRemoveNA(xiewei.data.exp)
xiewei.data.exp <- myNormalize(xiewei.data.exp)
boxplot(xiewei.data.exp,col="skyblue")

####----5.1.3 load LR pairs-------
####load LRpairs 1744
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20201214.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)


####load peLR 266 pairs
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)



#####-----------5.2 LR strength--------------

#####-----------5.2.1 load data-----------

####load LRpairs 1744
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20201214.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)

####load mean mat
data.mean.plot.mat <- readRDS(file = "res/R/data.plot.mean.mat_20201214.rds")

#####----------5.2.2 calculate IS---------

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")

tmp.mat.L <- data.mean.plot.mat[Lgene.list,] 
tmp.mat.R <- data.mean.plot.mat[Rgene.list,]

tmp.res <- (tmp.mat.L * tmp.mat.R) %>%
  rownames_to_column("LRid") %>%
  mutate(LRid = LRpairs) %>%
  gather(key = "Stage",value = "IS",-LRid) %>%
  mutate(log10IS=log10(IS),
         Stage=factor(Stage,levels = Stage.levels),
         group="LR")

myBoxPlot(data.plot = tmp.res,
          x = "Stage",
          y = "log10IS",
          xlab = "Stage",
          ylab = "Interaction Strength")


#####-------5.2.3 generate randomLR-----------
gene_symbols <- rownames(data.plot.mean.mat)
tmp.res.list <- list()
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20201214.rds")
for (ii in 1:100) {
  cat("simu:",ii,",start\n")
  gene_symbols <- rownames(data.plot.mean.mat)
  tmp.LR.pairs <- NULL
  while (length(tmp.LR.pairs) < length(LRpairs)) {
    tmp.Lgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.Rgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.LR.pairs <- paste0(tmp.Lgene.list,"_",tmp.Rgene.list)
  }
  tmp.res.df <- data.frame(stringsAsFactors = F)
  
  tmp.mat.L <- data.mean.plot.mat[tmp.Lgene.list,] 
  tmp.mat.R <- data.mean.plot.mat[tmp.Rgene.list,]
  
  tmp.res.df <- (tmp.mat.L * tmp.mat.R) %>%
    rownames_to_column("LRid") %>%
    mutate(LRid = LRpairs) %>%
    gather(key = "Stage",value = "IS",-LRid) %>%
    mutate(log10IS=log10(IS),
           Stage=factor(Stage,levels = Stage.levels),
           group="randomLR")
  tmp.res.df$simu_count <- paste0("simu_",ii)
  tmp.res.list <- c(tmp.res.list,list(tmp.res.df))
  cat("simu:",ii,",end\n")
}

tmp.data.plot <- tmp.res.list %>%
  purrr::reduce(rbind) %>%
  dplyr::select(-simu_count)

####---------5.2.4 plot boxplot -------------
tmp.data.plot <- rbind(tmp.res,tmp.data.plot)
mycolor <- c("#73c0dc","#f2f2f2")
tmp.stage <- c("MII oocyte","zygote",
               "early 2-cell","mid 2-cell","late 2-cell",
               "4-cell","8-cell","16-cell",
               "early blastocyst","mid blastocyst","late blastocyst",
               "E5.25","E5.5","E6.25","E6.5","fibroblast")

p <- ggplot(tmp.data.plot,aes(Stage,log10IS,fill=group,color=group))+
  geom_boxplot(linetype = "dashed", 
               outlier.shape = NA,width = 0.8 ,
               position = position_dodge(0.9)) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),
               outlier.shape=NA,width = 0.8 ,
               position = position_dodge(0.9)) +
  stat_boxplot(geom = "errorbar", 
               aes(ymin = ..ymax..),width = 0.8 ,
               position = position_dodge(0.9)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..),
               position = position_dodge(0.9)) +
  stat_boxplot(data = tmp.data.plot,
               geom = "errorbar", 
               aes(x = Stage,
                   y = log10IS,
                   ymin = ..middle..,
                   ymax = ..middle..),
               size = 1,
               col="black",
               width = 0.8 ,
               position = position_dodge(0.9))+
  theme_cowplot(font_size = 18)+
  ylab("Interaction Strength")+
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  scale_x_discrete(labels=tmp.stage)+
  scale_y_continuous(limits = c(-8,3),
                     breaks = seq(-8,3,by = 1),
                     labels = parse(text=paste0(10,"^",seq(-8,3,by = 1))))+
  theme(axis.line = element_line(size = 1),
        axis.text.x = element_text(hjust = 0.8,angle = 30))

myggsave(p,prefix = "res/fig/Fig5_IS_groupby",suffix = ".jpg",width = 8,height = 8,dpi=350)


####-----------5.3 eLR strength--------------

#####-----------5.2.1 load data-----------

####load LRpairs 239
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)

####load mean mat
data.mean.plot.mat <- readRDS(file = "res/R/data.plot.mean.mat_20201214.rds")

#####----------5.2.2 calculate IS---------

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")

tmp.mat.L <- data.mean.plot.mat[Lgene.list,] 
tmp.mat.R <- data.mean.plot.mat[Rgene.list,]

tmp.res <- (tmp.mat.L * tmp.mat.R) %>%
  rownames_to_column("LRid") %>%
  mutate(LRid = LRpairs) %>%
  gather(key = "Stage",value = "IS",-LRid) %>%
  mutate(log10IS=log10(IS),
         Stage=factor(Stage,levels = Stage.levels),
         group="eLR")

myBoxPlot(data.plot = tmp.res,
          x = "Stage",
          y = "log10IS",
          xlab = "Stage",
          ylab = "Interaction Strength")


#####-------5.2.3 generate randomLR-----------
gene_symbols <- rownames(data.plot.mean.mat)
tmp.res.list <- list()

for (ii in 1:100) {
  cat("simu:",ii,",start\n")
  gene_symbols <- rownames(data.plot.mean.mat)
  tmp.LR.pairs <- NULL
  while (length(tmp.LR.pairs) < length(LRpairs)) {
    tmp.Lgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.Rgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.LR.pairs <- paste0(tmp.Lgene.list,"_",tmp.Rgene.list)
  }
  tmp.res.df <- data.frame(stringsAsFactors = F)
  
  tmp.mat.L <- data.mean.plot.mat[tmp.Lgene.list,] 
  tmp.mat.R <- data.mean.plot.mat[tmp.Rgene.list,]
  
  tmp.res.df <- (tmp.mat.L * tmp.mat.R) %>%
    rownames_to_column("LRid") %>%
    mutate(LRid = LRpairs) %>%
    gather(key = "Stage",value = "IS",-LRid) %>%
    mutate(log10IS=log10(IS),
           Stage=factor(Stage,levels = Stage.levels),
           group="randomLR")
  tmp.res.df$simu_count <- paste0("simu_",ii)
  tmp.res.list <- c(tmp.res.list,list(tmp.res.df))
  cat("simu:",ii,",end\n")
}

tmp.data.plot <- tmp.res.list %>%
  purrr::reduce(rbind) %>%
  dplyr::select(-simu_count)

####---------5.2.4 plot boxplot -------------
tmp.data.plot <- rbind(tmp.res,tmp.data.plot)
mycolor <- c("#73c0dc","#f2f2f2")
tmp.stage <- c("MII oocyte","zygote",
               "early 2-cell","mid 2-cell","late 2-cell",
               "4-cell","8-cell","16-cell",
               "early blastocyst","mid blastocyst","late blastocyst",
               "E5.25","E5.5","E6.25","E6.5","fibroblast")

p <- ggplot(tmp.data.plot,aes(Stage,log10IS,fill=group,color=group))+
  geom_boxplot(linetype = "dashed", 
               outlier.shape = NA,width = 0.8 ,
               position = position_dodge(0.9)) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),
               outlier.shape=NA,width = 0.8 ,
               position = position_dodge(0.9)) +
  stat_boxplot(geom = "errorbar", 
               aes(ymin = ..ymax..),width = 0.8 ,
               position = position_dodge(0.9)) +
  stat_boxplot(geom = "errorbar", aes(ymax = ..ymin..),
               position = position_dodge(0.9)) +
  stat_boxplot(data = tmp.data.plot,
               geom = "errorbar", 
               aes(x = Stage,
                   y = log10IS,
                   ymin = ..middle..,
                   ymax = ..middle..),
               size = 1,
               col="black",
               width = 0.8 ,
               position = position_dodge(0.9))+
  theme_cowplot(font_size = 18)+
  ylab("Interaction Strength")+
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  scale_x_discrete(labels=tmp.stage)+
  scale_y_continuous(limits = c(-8,3),
                     breaks = seq(-8,3,by = 1),
                     labels = parse(text=paste0(10,"^",seq(-8,3,by = 1))))+
  theme(axis.line = element_line(size = 1),
        axis.text.x = element_text(hjust = 0.8,angle = 30))

myggsave(p,prefix = "res/fig/Fig5_eLR_IS_groupby",suffix = ".jpg",width = 8,height = 8,dpi=350)


#####-----------5.4 LR difference-------------------

#####------5.4(a) difference single-cell pseduobulk---------

#####------5.2.1 all LR pairs-----------
tmp.mat <- readRDS(file = "res/R/data.plot.mean.mat_20201214.rds")

LRpairs <- readRDS(file = "res/R/LRpairs_1619_20201117.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 


# Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
#                    "4cell","8cell","16cell","earlyblast","midblast",  
#                    "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast")


tmp.mat.L <- tmp.mat[Lgene.list,Stage.levels]
tmp.mat.R <- tmp.mat[Rgene.list,Stage.levels]
tmp.LR.pairs <- LRpairs

LR.diffference <- abs(tmp.mat.L-tmp.mat.R) %>%
  rownames_to_column("LRpairs") %>%
  dplyr::mutate(LRpairs= tmp.LR.pairs) %>%
  gather(key = "Stage",value=difference,-LRpairs) %>%
  dplyr::mutate(Group = "LR",
                Stage = factor(Stage,Stage.levels))
head(LR.diffference)
ggplot(LR.diffference,aes(Stage,difference))+
  geom_boxplot()

myBoxplot_advanced(LR.diffference,
                   x = "Stage",
                   y = "difference",
                   fill = "Stage",
                   color = "Stage",
                   mycolor = rep("#73bfdc",16),
                   ylimts = c(0,16),
                   ybreaks = seq(0,16,2),
                   outlier.shape = NA,
                   outlier.alpha = 0.1)

### add background

gene_symbols <- rownames(tmp.mat)

tmp.res.list <- list()
for(ii in 1:100){
  tmp.LR.pairs <- NULL
  while (length(tmp.LR.pairs) < length(LRpairs)) {
    tmp.Lgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.Rgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.LR.pairs <- paste0(tmp.Lgene.list,"_",tmp.Rgene.list)
  }
  tmp.mat.L <- tmp.mat[tmp.Lgene.list,Stage.levels]
  tmp.mat.R <- tmp.mat[tmp.Rgene.list,Stage.levels]
  
  
  tmp.LR.diffference <- abs(tmp.mat.L-tmp.mat.R) %>%
    rownames_to_column("LRpairs") %>%
    dplyr::mutate(LRpairs= tmp.LR.pairs) %>%
    gather(key = "Stage",value=difference,-LRpairs) %>%
    dplyr::mutate(Group = "randomLR",
                  Stage = factor(Stage,Stage.levels))
  
  tmp.res.list <- c(tmp.res.list,list(tmp.LR.diffference))
}

randomLR.difference <- tmp.res.list %>%
  purrr::reduce(rbind)


head(LR.diffference)
head(randomLR.difference)
difference.data.plot <- rbind(LR.diffference,
                              randomLR.difference)

mycolor <- c("#73bfdc","#f2f2f2")
p <- myBoxplot_advanced(difference.data.plot,
                        x = "Stage",
                        y = "difference",
                        fill = "Group",
                        color = "Group",
                        mycolor = mycolor,
                        ylimts = c(0,10),
                        ybreaks = seq(0,10,2),
                        outlier.shape = NA,
                        outlier.alpha = 0.1)
p
myggsave(p,prefix = "res/fig/difference_1619_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


test <- difference.data.plot %>%
  group_by(Stage,Group) %>%
  summarise(test = median(difference)) %>%
  ungroup()
test.1 <- test %>%
  dplyr::filter(Group=="LR")
test.2 <- test %>%
  dplyr::filter(Group=="randomLR")

test.res <- data.frame(Stage = test.1$Stage,
                       delta.diff = test.1$test-test.2$test)

p <- ggplot(data = test.res,
       aes(Stage,delta.diff,fill=Stage,color=Stage))+
  geom_point(stat = "identity",size=5)+
  geom_line(group="")+
  theme_cowplot(font_size = 15)+
  scale_color_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/difference_delta_1619_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


####---------5.2.2 266 peLR-------

tmp.mat <- readRDS(file = "res/R/data.plot.mean.mat.rds")

####load peLR 266 pairs
LRpairs <- readRDS(file = "res/R/eLR_pairs_266_20201120.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 


# Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
#                    "4cell","8cell","16cell","earlyblast","midblast",  
#                    "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast")


tmp.mat.L <- tmp.mat[Lgene.list,Stage.levels]
tmp.mat.R <- tmp.mat[Rgene.list,Stage.levels]
tmp.LR.pairs <- LRpairs

LR.diffference <- abs(tmp.mat.L-tmp.mat.R) %>%
  rownames_to_column("LRpairs") %>%
  dplyr::mutate(LRpairs= tmp.LR.pairs) %>%
  gather(key = "Stage",value=difference,-LRpairs) %>%
  dplyr::mutate(Group = "peLR",
                Stage = factor(Stage,Stage.levels))
head(LR.diffference)
ggplot(LR.diffference,aes(Stage,difference))+
  geom_boxplot()

myBoxplot_advanced(LR.diffference,
                   x = "Stage",
                   y = "difference",
                   fill = "Stage",
                   color = "Stage",
                   mycolor = rep("#73bfdc",16),
                   ylimts = c(0,16),
                   ybreaks = seq(0,16,2),
                   outlier.shape = NA,
                   outlier.alpha = 0.1)


gene_symbols <- rownames(tmp.mat)

tmp.res.list <- list()
for(ii in 1:100){
  tmp.LR.pairs <- NULL
  while (length(tmp.LR.pairs) < length(LRpairs)) {
    tmp.Lgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.Rgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.LR.pairs <- paste0(tmp.Lgene.list,"_",tmp.Rgene.list)
  }
  tmp.mat.L <- tmp.mat[tmp.Lgene.list,Stage.levels]
  tmp.mat.R <- tmp.mat[tmp.Rgene.list,Stage.levels]
  
  
  tmp.LR.diffference <- abs(tmp.mat.L-tmp.mat.R) %>%
    rownames_to_column("LRpairs") %>%
    dplyr::mutate(LRpairs= tmp.LR.pairs) %>%
    gather(key = "Stage",value=difference,-LRpairs) %>%
    dplyr::mutate(Group = "randomLR",
                  Stage = factor(Stage,Stage.levels))
  
  tmp.res.list <- c(tmp.res.list,list(tmp.LR.diffference))
}

randomLR.difference <- tmp.res.list %>%
  purrr::reduce(rbind)


difference.data.plot <- rbind(LR.diffference,
                              randomLR.difference)

mycolor <- c("#73bfdc","#f2f2f2")
p <- myBoxplot_advanced(difference.data.plot,
                        x = "Stage",
                        y = "difference",
                        fill = "Group",
                        color = "Group",
                        mycolor = mycolor,
                        ylimts = c(0,16),
                        ybreaks = seq(0,16,2),
                        outlier.shape = NA,
                        outlier.alpha = 0.1)
myggsave(p,prefix = "res/fig/difference_266_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)

####---------5.2.3  eLR-------

tmp.mat <- readRDS(file = "res/R/data.plot.mean.mat_20201214.rds")

tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")

# Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
#                    "4cell","8cell","16cell","earlyblast","midblast",  
#                    "lateblast")


tmp.mat.L <- tmp.mat[Lgene.list,Stage.levels]
tmp.mat.R <- tmp.mat[Rgene.list,Stage.levels]
tmp.LR.pairs <- LRpairs

LR.diffference <- abs(tmp.mat.L-tmp.mat.R) %>%
  rownames_to_column("LRpairs") %>%
  dplyr::mutate(LRpairs= tmp.LR.pairs) %>%
  gather(key = "Stage",value=difference,-LRpairs) %>%
  dplyr::mutate(Group = "eLR",
                Stage = factor(Stage,Stage.levels))
head(LR.diffference)
ggplot(LR.diffference,aes(Stage,difference))+
  geom_boxplot()

myBoxplot_advanced(LR.diffference,
                   x = "Stage",
                   y = "difference",
                   fill = "Stage",
                   color = "Stage",
                   mycolor = rep("#73bfdc",16),
                   ylimts = c(0,16),
                   ybreaks = seq(0,16,2),
                   outlier.shape = NA,
                   outlier.alpha = 0.1)


gene_symbols <- rownames(tmp.mat)

tmp.res.list <- list()
for(ii in 1:100){
  tmp.LR.pairs <- NULL
  while (length(tmp.LR.pairs) < length(LRpairs)) {
    tmp.Lgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.Rgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.LR.pairs <- paste0(tmp.Lgene.list,"_",tmp.Rgene.list)
  }
  tmp.mat.L <- tmp.mat[tmp.Lgene.list,Stage.levels]
  tmp.mat.R <- tmp.mat[tmp.Rgene.list,Stage.levels]
  
  
  tmp.LR.diffference <- abs(tmp.mat.L-tmp.mat.R) %>%
    rownames_to_column("LRpairs") %>%
    dplyr::mutate(LRpairs= tmp.LR.pairs) %>%
    gather(key = "Stage",value=difference,-LRpairs) %>%
    dplyr::mutate(Group = "randomLR",
                  Stage = factor(Stage,Stage.levels))
  
  tmp.res.list <- c(tmp.res.list,list(tmp.LR.diffference))
}

randomLR.difference <- tmp.res.list %>%
  purrr::reduce(rbind)


difference.data.plot <- rbind(LR.diffference,
                              randomLR.difference)

mycolor <- c("#73bfdc","#f2f2f2")
tmp.stage <- c("MII oocyte","zygote",
               "early 2-cell","mid 2-cell","late 2-cell",
               "4-cell","8-cell","16-cell",
               "early blastocyst","mid blastocyst","late blastocyst",
               "E5.25","E5.5","E6.25","E6.5","fibroblast")
p <- myBoxplot_advanced(difference.data.plot,
                        x = "Stage",
                        y = "difference",
                        fill = "Group",
                        color = "Group",
                        mycolor = mycolor,
                        ylimts = c(0,10),
                        ybreaks = seq(0,10,2),
                        outlier.shape = NA,
                        outlier.alpha = 0.1)+
  scale_x_discrete(label=tmp.stage)
p
myggsave(p,prefix = "res/fig/difference_239_eLR",suffix = ".jpg",width = 8,height = 8,dpi=300)


# test <- difference.data.plot %>%
#   group_by(Stage,Group) %>%
#   summarise(test = median(difference)) %>%
#   ungroup()
# test.1 <- test %>%
#   dplyr::filter(Group=="eLR")
# test.2 <- test %>%
#   dplyr::filter(Group=="randomLR")
# 
# test.res <- data.frame(Stage = test.1$Stage,
#                        delta.diff = test.1$test-test.2$test)
# 
# p <- ggplot(data = test.res,
#             aes(Stage,delta.diff,fill=Stage,color=Stage))+
#   geom_point(stat = "identity",size=5)+
#   geom_line(group="")+
#   theme_cowplot(font_size = 15)+
#   scale_color_manual(values = test.color.2(16))+
#   theme(legend.position = "none",
#         plot.title = element_text(hjust = 0.5),
#         legend.text = element_text( size = 15 ) ,
#         axis.text = element_text(size = 15),
#         axis.line.x = element_line( size = 1 ),
#         axis.line.y = element_line( size = 1 ),
#         axis.text.x = element_text(hjust = 0.8,angle = 30))
# p
# myggsave(p,prefix = "res/fig/difference_delta_210_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)



#####------5.2(b) difference single-cell pseduobulk---------

#####------5.2.1 1619 LR pairs-----------
tmp.mat <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
tmp.mat <- myNormalize(myRemoveNA(tmp.mat))
colnames(tmp.mat)[1] <- "Oocyte"
boxplot(tmp.mat)

LRpairs <- readRDS(file = "res/R/LRpairs_1619_20201117.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 


Stage.levels <-  c("Oocyte","zygote","early_2cell","2cell",
                   "4cell","8cell","ICM")


tmp.mat.L <- tmp.mat[Lgene.list,Stage.levels]
tmp.mat.R <- tmp.mat[Rgene.list,Stage.levels]
tmp.LR.pairs <- LRpairs

LR.diffference <- abs(tmp.mat.L-tmp.mat.R) %>%
  rownames_to_column("LRpairs") %>%
  dplyr::mutate(LRpairs= tmp.LR.pairs) %>%
  gather(key = "Stage",value=difference,-LRpairs) %>%
  dplyr::mutate(Group = "LR",
                Stage = factor(Stage,Stage.levels))
head(LR.diffference)
ggplot(LR.diffference,aes(Stage,difference))+
  geom_boxplot()

myBoxplot_advanced(LR.diffference,
                   x = "Stage",
                   y = "difference",
                   fill = "Stage",
                   color = "Stage",
                   mycolor = rep("#73bfdc",16),
                   ylimts = c(0,16),
                   ybreaks = seq(0,16,2),
                   outlier.shape = NA,
                   outlier.alpha = 0.1)

### add background

gene_symbols <- rownames(tmp.mat)

tmp.res.list <- list()
for(ii in 1:100){
  tmp.LR.pairs <- NULL
  while (length(tmp.LR.pairs) < length(LRpairs)) {
    tmp.Lgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.Rgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.LR.pairs <- paste0(tmp.Lgene.list,"_",tmp.Rgene.list)
  }
  tmp.mat.L <- tmp.mat[tmp.Lgene.list,Stage.levels]
  tmp.mat.R <- tmp.mat[tmp.Rgene.list,Stage.levels]
  
  tmp.LR.diffference <- abs(tmp.mat.L-tmp.mat.R) %>%
    rownames_to_column("LRpairs") %>%
    dplyr::mutate(LRpairs= tmp.LR.pairs) %>%
    gather(key = "Stage",value=difference,-LRpairs) %>%
    dplyr::mutate(Group = "randomLR",
                  Stage = factor(Stage,Stage.levels))
  
  tmp.res.list <- c(tmp.res.list,list(tmp.LR.diffference))
}

randomLR.difference <- tmp.res.list %>%
  purrr::reduce(rbind)


head(LR.diffference)
head(randomLR.difference)
difference.data.plot <- rbind(LR.diffference,
                              randomLR.difference)

mycolor <- c("#73bfdc","#f2f2f2")
p <- myBoxplot_advanced(difference.data.plot,
                        x = "Stage",
                        y = "difference",
                        fill = "Group",
                        color = "Group",
                        mycolor = mycolor,
                        ylimts = c(0,10),
                        ybreaks = seq(0,10,2),
                        outlier.shape = NA,
                        outlier.alpha = 0.1)
p
myggsave(p,prefix = "res/fig/xiweidata_difference_1619_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)

test <- difference.data.plot %>%
  group_by(Stage,Group) %>%
  summarise(test = median(difference)) %>%
  ungroup()
test.1 <- test %>%
  dplyr::filter(Group=="LR")
test.2 <- test %>%
  dplyr::filter(Group=="randomLR")

test.res <- data.frame(Stage = test.1$Stage,
                       delta.diff = test.1$test-test.2$test)

p <- ggplot(data = test.res,
            aes(Stage,delta.diff,fill=Stage,color=Stage))+
  geom_point(stat = "identity",size=5)+
  geom_line(group="")+
  theme_cowplot(font_size = 15)+
  scale_color_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/xieweidata_difference_delta_1619_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)




####---------5.2.2 266 peLR-------

tmp.mat <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
tmp.mat <- myNormalize(myRemoveNA(tmp.mat))
colnames(tmp.mat)[1] <- "Oocyte"
boxplot(tmp.mat)

####load peLR 266 pairs
LRpairs <- readRDS(file = "res/R/eLR_pairs_266_20201120.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early_2cell","2cell",
                   "4cell","8cell","ICM")


tmp.mat.L <- tmp.mat[Lgene.list,Stage.levels]
tmp.mat.R <- tmp.mat[Rgene.list,Stage.levels]
tmp.LR.pairs <- LRpairs

LR.diffference <- abs(tmp.mat.L-tmp.mat.R) %>%
  rownames_to_column("LRpairs") %>%
  dplyr::mutate(LRpairs= tmp.LR.pairs) %>%
  gather(key = "Stage",value=difference,-LRpairs) %>%
  dplyr::mutate(Group = "peLR",
                Stage = factor(Stage,Stage.levels))
head(LR.diffference)
ggplot(LR.diffference,aes(Stage,difference))+
  geom_boxplot()

myBoxplot_advanced(LR.diffference,
                   x = "Stage",
                   y = "difference",
                   fill = "Stage",
                   color = "Stage",
                   mycolor = rep("#73bfdc",16),
                   ylimts = c(0,16),
                   ybreaks = seq(0,16,2),
                   outlier.shape = NA,
                   outlier.alpha = 0.1)


gene_symbols <- rownames(tmp.mat)

tmp.res.list <- list()
for(ii in 1:100){
  tmp.LR.pairs <- NULL
  while (length(tmp.LR.pairs) < length(LRpairs)) {
    tmp.Lgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.Rgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.LR.pairs <- paste0(tmp.Lgene.list,"_",tmp.Rgene.list)
  }
  tmp.mat.L <- tmp.mat[tmp.Lgene.list,Stage.levels]
  tmp.mat.R <- tmp.mat[tmp.Rgene.list,Stage.levels]
  
  
  tmp.LR.diffference <- abs(tmp.mat.L-tmp.mat.R) %>%
    rownames_to_column("LRpairs") %>%
    dplyr::mutate(LRpairs= tmp.LR.pairs) %>%
    gather(key = "Stage",value=difference,-LRpairs) %>%
    dplyr::mutate(Group = "randomLR",
                  Stage = factor(Stage,Stage.levels))
  
  tmp.res.list <- c(tmp.res.list,list(tmp.LR.diffference))
}

randomLR.difference <- tmp.res.list %>%
  purrr::reduce(rbind)


difference.data.plot <- rbind(LR.diffference,
                              randomLR.difference)

mycolor <- c("#73bfdc","#f2f2f2")
p <- myBoxplot_advanced(difference.data.plot,
                        x = "Stage",
                        y = "difference",
                        fill = "Group",
                        color = "Group",
                        mycolor = mycolor,
                        ylimts = c(0,10),
                        ybreaks = seq(0,10,2),
                        outlier.shape = NA,
                        outlier.alpha = 0.1)
myggsave(p,prefix = "res/fig/xieweidata_difference_266_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)

####---------5.2.3 210 eLR-------


tmp.mat <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
tmp.mat <- myNormalize(myRemoveNA(tmp.mat))
colnames(tmp.mat)[1] <- "Oocyte"
boxplot(tmp.mat)

LRpairs <- readRDS(file = "res/R/eLR_correlation.Rds") %>%
  dplyr::filter(abs(PCC) > 0.2) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early_2cell","2cell",
                   "4cell","8cell","ICM")


tmp.mat.L <- tmp.mat[Lgene.list,Stage.levels]
tmp.mat.R <- tmp.mat[Rgene.list,Stage.levels]
tmp.LR.pairs <- LRpairs

LR.diffference <- abs(tmp.mat.L-tmp.mat.R) %>%
  rownames_to_column("LRpairs") %>%
  dplyr::mutate(LRpairs= tmp.LR.pairs) %>%
  gather(key = "Stage",value=difference,-LRpairs) %>%
  dplyr::mutate(Group = "eLR",
                Stage = factor(Stage,Stage.levels))
head(LR.diffference)
ggplot(LR.diffference,aes(Stage,difference))+
  geom_boxplot()

myBoxplot_advanced(LR.diffference,
                   x = "Stage",
                   y = "difference",
                   fill = "Stage",
                   color = "Stage",
                   mycolor = rep("#73bfdc",16),
                   ylimts = c(0,16),
                   ybreaks = seq(0,16,2),
                   outlier.shape = NA,
                   outlier.alpha = 0.1)


gene_symbols <- rownames(tmp.mat)

tmp.res.list <- list()
for(ii in 1:100){
  tmp.LR.pairs <- NULL
  while (length(tmp.LR.pairs) < length(LRpairs)) {
    tmp.Lgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.Rgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.LR.pairs <- paste0(tmp.Lgene.list,"_",tmp.Rgene.list)
  }
  tmp.mat.L <- tmp.mat[tmp.Lgene.list,Stage.levels]
  tmp.mat.R <- tmp.mat[tmp.Rgene.list,Stage.levels]
  
  
  tmp.LR.diffference <- abs(tmp.mat.L-tmp.mat.R) %>%
    rownames_to_column("LRpairs") %>%
    dplyr::mutate(LRpairs= tmp.LR.pairs) %>%
    gather(key = "Stage",value=difference,-LRpairs) %>%
    dplyr::mutate(Group = "randomLR",
                  Stage = factor(Stage,Stage.levels))
  
  tmp.res.list <- c(tmp.res.list,list(tmp.LR.diffference))
}

randomLR.difference <- tmp.res.list %>%
  purrr::reduce(rbind)


difference.data.plot <- rbind(LR.diffference,
                              randomLR.difference)

mycolor <- c("#73bfdc","#f2f2f2")
p <- myBoxplot_advanced(difference.data.plot,
                        x = "Stage",
                        y = "difference",
                        fill = "Group",
                        color = "Group",
                        mycolor = mycolor,
                        ylimts = c(0,10),
                        ybreaks = seq(0,10,2),
                        outlier.shape = NA,
                        outlier.alpha = 0.1)
p
myggsave(p,prefix = "res/fig/xieweidata_difference_210_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


test <- difference.data.plot %>%
  group_by(Stage,Group) %>%
  summarise(test = median(difference)) %>%
  ungroup()
test.1 <- test %>%
  dplyr::filter(Group=="eLR")
test.2 <- test %>%
  dplyr::filter(Group=="randomLR")

test.res <- data.frame(Stage = test.1$Stage,
                       delta.diff = test.1$test-test.2$test)

p <- ggplot(data = test.res,
            aes(Stage,delta.diff,fill=Stage,color=Stage))+
  geom_point(stat = "identity",size=5)+
  geom_line(group="")+
  theme_cowplot(font_size = 15)+
  scale_color_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/xieweidata_difference_delta_210_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


####-------5.3 low high level expression (a) singlecell----------

####-------5.3.0 prepare state mat---------

tmp.mat <- readRDS(file = "res/R/data.plot.mean.mat.rds")
boxplot.stats(tmp.mat)
?boxplot.stats

tmp.tttt <- boxplot(tmp.mat)$stats
apply(tmp.tttt, 1, min)

tmp.threshold <- c(0.10,5.15)


tmp.mat.long <- tmp.mat %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = gene.exp.mean,-Gene) %>%
  mutate(Gene_state = ifelse(gene.exp.mean < tmp.threshold[1],
                             "low",
                             ifelse((gene.exp.mean >= tmp.threshold[1]) &  (gene.exp.mean < tmp.threshold[2]),
                                    "middle","high")))
# head(tmp.mat.long)
Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")

tmp.mat.long$Gene_state <- factor(tmp.mat.long$Gene_state,
                                  levels = c("low","middle","high"))
tmp.mat.long$Stage <- factor(tmp.mat.long$Stage,levels = Stage.levels)

tmp.summary <- tmp.mat.long %>%
  group_by(Stage,Gene_state) %>%
  summarise(n=n()) %>%
  ungroup()
tmp.summary$Gene_state <- factor(tmp.summary$Gene_state,
                                 levels = rev(c("low","middle","high")))

mycolor <- c("#DEEBF7","#9ECAE1","#3182BD")
p <- myBoxplot_advanced(tmp.mat.long,
                        x = "Stage",
                        y = "gene.exp.mean",
                        fill = "Gene_state",
                        color = "Gene_state",
                        mycolor = mycolor,
                        ylimts = c(0,16),
                        ybreaks = seq(0,16,2))

myggsave(p,prefix = "res/fig/gene_expresion_levels",suffix = ".jpg",width = 8,height = 8,dpi=300)

p <- ggplot(tmp.summary,
            aes(Stage,n,fill=Gene_state,label=n))+
  geom_col(position = "stack")+
  #geom_text()+
  ylab("count")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = rev(mycolor))+
  scale_y_continuous(limits = c(0,20000),
                     breaks = seq(0,20000,1000),
                     expand = c(0,0.1))+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p

myggsave(p,prefix = "res/fig/gene_expresion_levels_stat",suffix = ".jpg",width = 8,height = 8,dpi=300)



###LR coordination
head(tmp.mat.long)

tmp.state.mat <- tmp.mat.long %>%
  dplyr::select(-gene.exp.mean) %>%
  spread(key = "Stage",value = Gene_state) %>%
  column_to_rownames("Gene")

head(tmp.state.mat[,1:3])


saveRDS(tmp.state.mat,file = "res/R/early_exp_state_mat_20201130.rds")

####---------5.3.1 LR 1619 ---------
tmp.mat <- readRDS(file = "res/R/early_exp_state_mat_20201130.rds")

LRpairs <- readRDS(file = "res/R/NATMI_LRpairs_1619LR_20201103.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")


tmp.mat.L <- tmp.mat[Lgene.list,]
tmp.mat.R <- tmp.mat[Rgene.list,]
tmp.LR.pairs <- LRpairs

head(tmp.mat.L)
head(tmp.mat.R)

tmp.res.list <- sapply(Stage.levels,function(x) paste0(tmp.mat.L[,x],"+",tmp.mat.R[,x]),USE.NAMES = T)
tmp.res.list <- as.data.frame(tmp.res.list)
rownames(tmp.res.list) <- LRpairs


tmp.levels <- c()
tmp.temp <- c("low","middle","high")
for (ii in tmp.temp) {
  tmp.res <- paste0(ii,"+",tmp.temp)
  tmp.levels <- c(tmp.levels,tmp.res)
}
tmp.levels


tmp.res.list.df <- tmp.res.list %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "Stage",value = "State",-LRpairs) %>%
  mutate(Stage=factor(Stage,levels = Stage.levels),
         State=factor(State,levels = tmp.levels))

saveRDS(tmp.res.list.df,file = "res/R/LR_expression_state_20201130.rds")

tmp.res.list.df.summ <- tmp.res.list.df %>%
  group_by(Stage,State) %>%
  summarise(count=n()) %>%
  ungroup()

tmp.res.list.df.summ$State <- factor(tmp.res.list.df.summ$State,levels = rev(tmp.levels))
tmp.res.list.df.summ$count <- tmp.res.list.df.summ$count / length(LRpairs)

mycolor <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues")) 
p <- ggplot(tmp.res.list.df.summ,
            aes(Stage,
                count,
                fill=State,
                label=count))+
  geom_bar(stat = "identity",position = "fill")+
  #geom_text()+
  ylab("Percent")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = rev(mycolor(9)))+
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.1),
                     expand = c(0,0))+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/gene_expresion_levels_state_coord",suffix = ".jpg",width = 8,height = 8,dpi=300)

####-----5.3.2 LR 266----
tmp.mat <- readRDS(file = "res/R/early_exp_state_mat_20201130.rds")

####load peLR 266 pairs
LRpairs <- readRDS(file = "res/R/eLR_pairs_266_20201120.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")


tmp.mat.L <- tmp.mat[Lgene.list,]
tmp.mat.R <- tmp.mat[Rgene.list,]
tmp.LR.pairs <- LRpairs

head(tmp.mat.L)
head(tmp.mat.R)

tmp.res.list <- sapply(Stage.levels,function(x) paste0(tmp.mat.L[,x],"+",tmp.mat.R[,x]),USE.NAMES = T)
tmp.res.list <- as.data.frame(tmp.res.list)
rownames(tmp.res.list) <- LRpairs


tmp.levels <- c()
tmp.temp <- c("low","middle","high")
for (ii in tmp.temp) {
  tmp.res <- paste0(ii,"+",tmp.temp)
  tmp.levels <- c(tmp.levels,tmp.res)
}
tmp.levels


tmp.res.list.df <- tmp.res.list %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "Stage",value = "State",-LRpairs) %>%
  mutate(Stage=factor(Stage,levels = Stage.levels),
         State=factor(State,levels = tmp.levels))

saveRDS(tmp.res.list.df,file = "res/R/LR_266_expression_state_20201130.rds")

tmp.res.list.df.summ <- tmp.res.list.df %>%
  group_by(Stage,State) %>%
  summarise(count=n()) %>%
  ungroup()

tmp.res.list.df.summ$State <- factor(tmp.res.list.df.summ$State,levels = rev(tmp.levels))
tmp.res.list.df.summ$count <- tmp.res.list.df.summ$count / length(LRpairs)

mycolor <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues")) 
p <- ggplot(tmp.res.list.df.summ,
            aes(Stage,
                count,
                fill=State,
                label=count))+
  geom_bar(stat = "identity",position = "fill")+
  #geom_text()+
  ylab("Percent")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = rev(mycolor(9)))+
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.1),
                     expand = c(0,0))+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/peLR_266_expresion_levels_state_coord",suffix = ".jpg",width = 8,height = 8,dpi=300)

####-----5.3.3 LR 210----
tmp.mat <- readRDS(file = "res/R/early_exp_state_mat_20201130.rds")

LRpairs <- readRDS(file = "res/R/eLR_correlation.Rds") %>%
  dplyr::filter(abs(PCC) > 0.2) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")


tmp.mat.L <- tmp.mat[Lgene.list,]
tmp.mat.R <- tmp.mat[Rgene.list,]
tmp.LR.pairs <- LRpairs

head(tmp.mat.L)
head(tmp.mat.R)

tmp.res.list <- sapply(Stage.levels,function(x) paste0(tmp.mat.L[,x],"+",tmp.mat.R[,x]),USE.NAMES = T)
tmp.res.list <- as.data.frame(tmp.res.list)
rownames(tmp.res.list) <- LRpairs


tmp.levels <- c()
tmp.temp <- c("low","middle","high")
for (ii in tmp.temp) {
  tmp.res <- paste0(ii,"+",tmp.temp)
  tmp.levels <- c(tmp.levels,tmp.res)
}
tmp.levels


tmp.res.list.df <- tmp.res.list %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "Stage",value = "State",-LRpairs) %>%
  mutate(Stage=factor(Stage,levels = Stage.levels),
         State=factor(State,levels = tmp.levels))

saveRDS(tmp.res.list.df,file = "res/R/LR_210_expression_state_20201130.rds")

tmp.res.list.df.summ <- tmp.res.list.df %>%
  group_by(Stage,State) %>%
  summarise(count=n()) %>%
  ungroup()

tmp.res.list.df.summ$State <- factor(tmp.res.list.df.summ$State,levels = rev(tmp.levels))
tmp.res.list.df.summ$count <- tmp.res.list.df.summ$count / length(LRpairs)

mycolor <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues")) 
p <- ggplot(tmp.res.list.df.summ,
            aes(Stage,
                count,
                fill=State,
                label=count))+
  geom_bar(stat = "identity",position = "fill")+
  #geom_text()+
  ylab("Percent")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = rev(mycolor(9)))+
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.1),
                     expand = c(0,0))+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/peLR_210_expresion_levels_state_coord",suffix = ".jpg",width = 8,height = 8,dpi=300)


####------5.3 stat summary--------

####------5.3.1 1619 LR pairs-----------

tmp.res.list.df <- readRDS(file = "res/R/LR_expression_state_20201130.rds")
LRpairs <- readRDS(file = "res/R/NATMI_LRpairs_1619LR_20201103.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 



tmp.res.list.df.summ <- tmp.res.list.df %>%
  group_by(Stage,State) %>%
  summarise(count=n()) %>%
  ungroup()

tmp.res.list.df.summ$State <- factor(tmp.res.list.df.summ$State,levels = rev(tmp.levels))
tmp.res.list.df.summ$count <- tmp.res.list.df.summ$count / length(LRpairs)

tmp.state <- tmp.res.list.df.summ$State
tmp.state <- strsplit(as.character(tmp.state),split = "\\+")
tmp.L.state <- unlist(lapply(tmp.state, function(x) x[1]))
tmp.R.state <- unlist(lapply(tmp.state, function(x) x[2]))
tmp.res.list.df.summ$L.state <- tmp.L.state
tmp.res.list.df.summ$R.state <- tmp.R.state
tmp.res.list.df.summ <- tmp.res.list.df.summ %>%
  mutate(Group=ifelse(tmp.L.state==tmp.R.state,"consistent","inconsistent"))

data.plot <- tmp.res.list.df.summ %>%
  group_by(Stage,Group) %>%
  summarise(tmp.sum=sum(count)) %>%
  ungroup()

mycolor <- c("#a6cee3","#1f78b4")
p <- ggplot(data.plot,aes(Stage,tmp.sum,fill=Group,label=round(tmp.sum,2)))+
  geom_bar(stat = "identity",position = position_dodge(0.8))+
  geom_text(vjust = -0.2,position = position_dodge(0.8),size = 4)+
  ylab("Percent")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = mycolor)+
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.1),
                     expand = c(0,0))+
  scale_x_discrete(expand = c(0.05,0))+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/consistence_LR_1619_expresion_levels_state",suffix = ".jpg",width = 8,height = 8,dpi=300)

####----------5.3.2 210 LR ----------------

tmp.res.list.df <- readRDS(file = "res/R/LR_210_expression_state_20201130.rds")

LRpairs <- readRDS(file = "res/R/eLR_correlation.Rds") %>%
  dplyr::filter(abs(PCC) > 0.2) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 


tmp.res.list.df.summ <- tmp.res.list.df %>%
  group_by(Stage,State) %>%
  summarise(count=n()) %>%
  ungroup()

tmp.res.list.df.summ$State <- factor(tmp.res.list.df.summ$State,levels = rev(tmp.levels))
tmp.res.list.df.summ$count <- tmp.res.list.df.summ$count / length(LRpairs)

tmp.state <- tmp.res.list.df.summ$State
tmp.state <- strsplit(as.character(tmp.state),split = "\\+")
tmp.L.state <- unlist(lapply(tmp.state, function(x) x[1]))
tmp.R.state <- unlist(lapply(tmp.state, function(x) x[2]))
tmp.res.list.df.summ$L.state <- tmp.L.state
tmp.res.list.df.summ$R.state <- tmp.R.state
tmp.res.list.df.summ <- tmp.res.list.df.summ %>%
  mutate(Group=ifelse(tmp.L.state==tmp.R.state,"consistent","inconsistent"))

data.plot <- tmp.res.list.df.summ %>%
  group_by(Stage,Group) %>%
  summarise(tmp.sum=sum(count)) %>%
  ungroup()

mycolor <- c("#a6cee3","#1f78b4")
p <- ggplot(data.plot,aes(Stage,tmp.sum,fill=Group,label=round(tmp.sum,2)))+
  geom_bar(stat = "identity",position = position_dodge(0.8))+
  geom_text(vjust = -0.2,position = position_dodge(0.8),size = 4)+
  ylab("Percent")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = mycolor)+
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.1),
                     expand = c(0,0))+
  scale_x_discrete(expand = c(0.05,0))+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/consistence_LR_210_expresion_levels_state",suffix = ".jpg",width = 8,height = 8,dpi=300)

##### sub by group-------
mycolor <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues")) 

tmp.levels <- c()
tmp.temp <- c("low","middle","high")
for (ii in tmp.temp) {
  tmp.res <- paste0(ii,"+",tmp.temp)
  tmp.levels <- c(tmp.levels,tmp.res)
}

for(ii in 1:length(tmp.levels)){
  p <- ggplot(subset(tmp.res.list.df.summ,
                     State==tmp.levels[ii]),
              aes(Stage,count,fill=State,label=round(count,2)))+
    geom_bar(stat = "identity",position = position_dodge(0.8))+
    geom_text(vjust = -0.2,position = position_dodge(0.8),size = 4)+
    ylab("Percent")+
    theme_cowplot(font_size = 15)+
    scale_fill_manual(values = mycolor(9)[ii])+
    scale_y_continuous(limits = c(0,1),
                       breaks = seq(0,1,0.1),
                       expand = c(0,0))+
    scale_x_discrete(expand = c(0.05,0))+
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text( size = 15 ) ,
          axis.text = element_text(size = 15),
          axis.line.x = element_line( size = 1 ),
          axis.line.y = element_line( size = 1 ),
          axis.text.x = element_text(hjust = 0.8,angle = 30))
  p
  myggsave(p,prefix = paste0("res/fig/consistence_LR_210_expresion_levels_state_",tmp.levels[ii]),suffix = ".jpg",width = 8,height = 8,dpi=300)
}



####----------5.3 (b) xieweidata ----------


####---------5.3.0 load data-----------
tmp.mat <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
tmp.mat <- myNormalize(myRemoveNA(tmp.mat))
colnames(tmp.mat)[1] <- "Oocyte"
boxplot(tmp.mat)$stat

tmp.tttt <- boxplot(tmp.mat)$stats
apply(tmp.tttt, 1, min)

tmp.threshold <- c(0.10,5.15)

Stage.levels <-  c("Oocyte","zygote","early_2cell","2cell",
                   "4cell","8cell","ICM")

tmp.mat.long <- tmp.mat[,Stage.levels] %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = gene.exp.mean,-Gene) %>%
  mutate(Gene_state = ifelse(gene.exp.mean < tmp.threshold[1],
                             "low",
                             ifelse((gene.exp.mean >= tmp.threshold[1]) &  (gene.exp.mean < tmp.threshold[2]),
                                    "middle","high")))
# head(tmp.mat.long)


tmp.mat.long$Gene_state <- factor(tmp.mat.long$Gene_state,
                                  levels = c("low","middle","high"))
tmp.mat.long$Stage <- factor(tmp.mat.long$Stage,levels = Stage.levels)

tmp.summary <- tmp.mat.long %>%
  group_by(Stage,Gene_state) %>%
  summarise(n=n()) %>%
  ungroup()
tmp.summary$Gene_state <- factor(tmp.summary$Gene_state,
                                 levels = rev(c("low","middle","high")))

mycolor <- c("#DEEBF7","#9ECAE1","#3182BD")
p <- myBoxplot_advanced(tmp.mat.long,
                        x = "Stage",
                        y = "gene.exp.mean",
                        fill = "Gene_state",
                        color = "Gene_state",
                        mycolor = mycolor,
                        ylimts = c(0,16),
                        ybreaks = seq(0,16,2))
p
myggsave(p,prefix = "res/fig/xieweidata_gene_expresion_levels",suffix = ".jpg",width = 8,height = 8,dpi=300)

p <- ggplot(tmp.summary,
            aes(Stage,n,fill=Gene_state,label=n))+
  geom_col(position = "stack")+
  #geom_text()+
  ylab("count")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = rev(mycolor))+
  scale_y_continuous(limits = c(0,30000),
                     breaks = seq(0,30000,2000),
                     expand = c(0,0.1))+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p

myggsave(p,prefix = "res/fig/xieweidata_gene_expresion_levels_stat",suffix = ".jpg",width = 8,height = 8,dpi=300)



###LR coordination
head(tmp.mat.long)
tmp.state.mat <- tmp.mat.long %>%
  dplyr::select(-gene.exp.mean) %>%
  spread(key = "Stage",value = Gene_state) %>%
  column_to_rownames("Gene")
head(tmp.state.mat[,1:3])
saveRDS(tmp.state.mat,file = "res/R/xieweidata_early_exp_state_mat_20201130.rds")

######--------5.3.1 1619 LR-------------

tmp.mat <- readRDS(file = "res/R/xieweidata_early_exp_state_mat_20201130.rds")

LRpairs <- readRDS(file = "res/R/NATMI_LRpairs_1619LR_20201103.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early_2cell","2cell",
                   "4cell","8cell","ICM")


tmp.mat.L <- tmp.mat[Lgene.list,]
tmp.mat.R <- tmp.mat[Rgene.list,]
tmp.LR.pairs <- LRpairs

head(tmp.mat.L)
head(tmp.mat.R)

tmp.res.list <- sapply(Stage.levels,function(x) paste0(tmp.mat.L[,x],"+",tmp.mat.R[,x]),USE.NAMES = T)
tmp.res.list <- as.data.frame(tmp.res.list)
rownames(tmp.res.list) <- LRpairs


tmp.levels <- c()
tmp.temp <- c("low","middle","high")
for (ii in tmp.temp) {
  tmp.res <- paste0(ii,"+",tmp.temp)
  tmp.levels <- c(tmp.levels,tmp.res)
}
tmp.levels


tmp.res.list.df <- tmp.res.list %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "Stage",value = "State",-LRpairs) %>%
  mutate(Stage=factor(Stage,levels = Stage.levels),
         State=factor(State,levels = tmp.levels))

saveRDS(tmp.res.list.df,file = "res/R/xieweidata_LR_expression_state_20201130.rds")

tmp.res.list.df.summ <- tmp.res.list.df %>%
  group_by(Stage,State) %>%
  summarise(count=n()) %>%
  ungroup()

tmp.res.list.df.summ$State <- factor(tmp.res.list.df.summ$State,levels = rev(tmp.levels))
tmp.res.list.df.summ$count <- tmp.res.list.df.summ$count / length(LRpairs)

mycolor <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues")) 
p <- ggplot(tmp.res.list.df.summ,
            aes(Stage,
                count,
                fill=State,
                label=count))+
  geom_bar(stat = "identity",position = "fill")+
  #geom_text()+
  ylab("Percent")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = rev(mycolor(9)))+
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.1),
                     expand = c(0,0))+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/xieweidata_gene_expresion_levels_state_coord",suffix = ".jpg",width = 8,height = 8,dpi=300)


######--------5.3.2 210 LR-------------

tmp.mat <- readRDS(file = "res/R/xieweidata_early_exp_state_mat_20201130.rds")

LRpairs <- readRDS(file = "res/R/eLR_correlation.Rds") %>%
  dplyr::filter(abs(PCC) > 0.2) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 


Stage.levels <-  c("Oocyte","zygote","early_2cell","2cell",
                   "4cell","8cell","ICM")


tmp.mat.L <- tmp.mat[Lgene.list,]
tmp.mat.R <- tmp.mat[Rgene.list,]
tmp.LR.pairs <- LRpairs

head(tmp.mat.L)
head(tmp.mat.R)

tmp.res.list <- sapply(Stage.levels,function(x) paste0(tmp.mat.L[,x],"+",tmp.mat.R[,x]),USE.NAMES = T)
tmp.res.list <- as.data.frame(tmp.res.list)
rownames(tmp.res.list) <- LRpairs


tmp.levels <- c()
tmp.temp <- c("low","middle","high")
for (ii in tmp.temp) {
  tmp.res <- paste0(ii,"+",tmp.temp)
  tmp.levels <- c(tmp.levels,tmp.res)
}
tmp.levels


tmp.res.list.df <- tmp.res.list %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "Stage",value = "State",-LRpairs) %>%
  mutate(Stage=factor(Stage,levels = Stage.levels),
         State=factor(State,levels = tmp.levels))

saveRDS(tmp.res.list.df,file = "res/R/xieweidata_eLR_210_expression_state_20201130.rds")

tmp.res.list.df.summ <- tmp.res.list.df %>%
  group_by(Stage,State) %>%
  summarise(count=n()) %>%
  ungroup()

tmp.res.list.df.summ$State <- factor(tmp.res.list.df.summ$State,levels = rev(tmp.levels))
tmp.res.list.df.summ$count <- tmp.res.list.df.summ$count / length(LRpairs)

mycolor <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues")) 
p <- ggplot(tmp.res.list.df.summ,
            aes(Stage,
                count,
                fill=State,
                label=count))+
  geom_bar(stat = "identity",position = "fill")+
  #geom_text()+
  ylab("Percent")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = rev(mycolor(9)))+
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.1),
                     expand = c(0,0))+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/xieweidata_210_gene_expresion_levels_state_coord",suffix = ".jpg",width = 8,height = 8,dpi=300)



tmp.state <- tmp.res.list.df.summ$State
tmp.state <- strsplit(as.character(tmp.state),split = "\\+")
tmp.L.state <- unlist(lapply(tmp.state, function(x) x[1]))
tmp.R.state <- unlist(lapply(tmp.state, function(x) x[2]))
tmp.res.list.df.summ$L.state <- tmp.L.state
tmp.res.list.df.summ$R.state <- tmp.R.state
tmp.res.list.df.summ <- tmp.res.list.df.summ %>%
  mutate(Group=ifelse(tmp.L.state==tmp.R.state,"consistent","inconsistent"))

data.plot <- tmp.res.list.df.summ %>%
  group_by(Stage,Group) %>%
  summarise(tmp.sum=sum(count)) %>%
  ungroup()

mycolor <- c("#a6cee3","#1f78b4")
p <- ggplot(data.plot,aes(Stage,tmp.sum,fill=Group,label=round(tmp.sum,2)))+
  geom_bar(stat = "identity",position = position_dodge(0.8))+
  geom_text(vjust = -0.2,position = position_dodge(0.8),size = 4)+
  ylab("Percent")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = mycolor)+
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.1),
                     expand = c(0,0))+
  scale_x_discrete(expand = c(0.05,0))+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/xieweidata_consistence_LR_210_expresion_levels_state",suffix = ".jpg",width = 8,height = 8,dpi=300)



### consitant

#####-------sub by group-------
mycolor <- colorRampPalette(RColorBrewer::brewer.pal(9,"Blues")) 

tmp.levels <- c()
tmp.temp <- c("low","middle","high")
for (ii in tmp.temp) {
  tmp.res <- paste0(ii,"+",tmp.temp)
  tmp.levels <- c(tmp.levels,tmp.res)
}

for(ii in 1:length(tmp.levels)){
  p <- ggplot(subset(tmp.res.list.df.summ,
                     State==tmp.levels[ii]),
              aes(Stage,count,fill=State,label=round(count,2)))+
    geom_bar(stat = "identity",position = position_dodge(0.8))+
    geom_text(vjust = -0.2,position = position_dodge(0.8),size = 4)+
    xlab("")+
    ylab("Percent")+
    theme_cowplot(font_size = 15)+
    scale_fill_manual(values = mycolor(9)[ii])+
    scale_y_continuous(limits = c(0,1),
                       breaks = seq(0,1,0.1),
                       expand = c(0,0))+
    scale_x_discrete(expand = c(0.05,0))+
    theme(legend.position = "top",
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text( size = 15 ) ,
          axis.text = element_text(size = 15),
          axis.line.x = element_line( size = 1 ),
          axis.line.y = element_line( size = 1 ),
          axis.text.x = element_text(hjust = 0.8,angle = 30))
  p
  myggsave(p,prefix = paste0("res/fig/xieweidata_consistence_LR_210_expresion_levels_state_",tmp.levels[ii]),suffix = ".jpg",width = 8,height = 8,dpi=300)
}


#####-------x.5.4 coordination score ?????----------

####-------x.5.4.1 load data-----------
tmp.mat <- readRDS(file = "res/R/data.plot.mean.mat_20201214.rds")
tmp.mat.long <- tmp.mat %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = gene.exp.mean,-Gene) %>%
  mutate(Gene_state = cut_interval(gene.exp.mean,10))
levels(tmp.mat.long$Gene_state) <- 1:10
tmp.mat.long$Gene_state <- as.integer(tmp.mat.long$Gene_state)
head(tmp.mat.long)

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")

tmp.mat.long$Stage <- factor(tmp.mat.long$Stage,Stage.levels)

tmp.state.mat <- tmp.mat.long %>%
  dplyr::select(-gene.exp.mean) %>%
  spread(key = "Stage",value = Gene_state) %>%
  column_to_rownames("Gene")

head(tmp.state.mat[,1:3])

saveRDS(tmp.state.mat,file = "res/R/early_exp_state_mat_20201216.rds")

####------x.5.4.2 calculate consistent score--------------

####---------x5.4.2.1 eLR---------
tmp.mat <- readRDS(file = "res/R/early_exp_state_mat_20201216.rds")
head(tmp.mat)

tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")


tmp.mat.L <- tmp.mat[Lgene.list,]
tmp.mat.R <- tmp.mat[Rgene.list,]
tmp.LR.pairs <- LRpairs

mat <- data.frame(stringsAsFactors = F)
for (ii in 1:10) {
  for (jj in 1:10) {
    mat[ii,jj] <- ii * jj / 100 * exp(-abs(ii-jj))
  }
}

mycolor <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")) 
tmp.cutoff <- 0.01
mat[mat>tmp.cutoff] <- 1
mat[mat<=tmp.cutoff] <- 0
jpeg(filename = myFileName(prefix = "res/fig/consistent_score_levels",suffix = ".jpg"),width = 8,height = 8,units = "in",res = 300)
pheatmap(mat = mat,
         color = c("#f2f2f2","#74bfdc"),
         labels_col = 1:10,
         angle_col = 0,fontsize = 18,
         legend_labels = NA,
         legend_breaks = c(0,1),
         show_rownames = T,
         show_colnames = T,
         cluster_rows = F,
         cluster_cols = F)
dev.off()

tmp.res.list <- sapply(Stage.levels,function(x) {
  tmp.res <- (tmp.mat.L[,x] * tmp.mat.R[,x])/100 * exp(-abs(tmp.mat.L[,x]-tmp.mat.R[,x]))
  tmp.res <- ifelse(tmp.res > 0.01,"consistent","inconsistent")
  },USE.NAMES = T)
tmp.res.list <- as.data.frame(tmp.res.list)
rownames(tmp.res.list) <- LRpairs

data.plot <- tmp.res.list %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "Stage",value = "state",-LRpairs) %>%
  mutate(Stage=factor(Stage,levels = Stage.levels)) %>%
  group_by(Stage,state) %>%
  summarise(percent=n()) %>%
  mutate(percent = percent/length(LRpairs)) %>%
  ungroup()
tmp.stage <- c("MII oocyte","zygote",
               "early 2-cell","mid 2-cell","late 2-cell",
               "4-cell","8-cell","16-cell",
               "early blastocyst","mid blastocyst","late blastocyst",
               "E5.25","E5.5","E6.25","E6.5","fibroblast")
p <- ggplot(subset(data.plot,state=="consistent"),aes(Stage,percent,fill=state,label=round(percent,2))) +
  geom_bar(stat="identity")+
  geom_text(vjust = -0.2,position = position_dodge(0.8),size = 4)+
  ylab("Consistent pairs Frequency")+
  theme_cowplot(font_size = 18)+
  #scale_fill_manual(values = mycolor(9)[ii])+
  scale_y_continuous(limits = c(0,0.3),
                     breaks = seq(0,0.3,0.1),
                     expand = c(0,0))+
  scale_x_discrete(label = tmp.stage,expand = c(0.05,0))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 18 ) ,
        axis.text = element_text(size = 18),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/consistence_score_eLR_expresion_levels_state_",suffix = ".jpg",width = 8,height = 8,dpi=300)


####---------x5.4.2.1 LR 210 ---------
tmp.mat <- readRDS(file = "res/R/early_exp_state_mat_20201203.rds")
head(tmp.mat)

LRpairs <- readRDS(file = "res/R/eLR_correlation.Rds") %>%
  dplyr::filter(abs(PCC) > 0.2) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 


Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")


tmp.mat.L <- tmp.mat[Lgene.list,]
tmp.mat.R <- tmp.mat[Rgene.list,]
tmp.LR.pairs <- LRpairs

mat <- data.frame(stringsAsFactors = F)
for (ii in 1:10) {
  for (jj in 1:10) {
    mat[ii,jj] <- ii * jj / 100 * exp(-abs(ii-jj))
  }
}

mycolor <- colorRampPalette(RColorBrewer::brewer.pal(9,"YlOrRd")) 
jpeg(filename = myFileName(prefix = "res/fig/consistent_score_levels",suffix = ".jpg"),width = 8,height = 8,units = "in",res = 300)
pheatmap(mat = mat,
         color = mycolor(100),
         labels_col = 1:10,
         angle_col = 0,
         fontsize = 15,
         show_rownames = T,
         show_colnames = T,
         cluster_rows = F,
         cluster_cols = F)
dev.off()

tmp.res.list <- sapply(Stage.levels,function(x) {
  tmp.res <- (tmp.mat.L[,x] * tmp.mat.R[,x])/100 * exp(-abs(tmp.mat.L[,x]-tmp.mat.R[,x]))
  tmp.res <- ifelse(tmp.res > 0.01,"consistent","inconsistent")
},USE.NAMES = T)
tmp.res.list <- as.data.frame(tmp.res.list)
rownames(tmp.res.list) <- LRpairs

data.plot <- tmp.res.list %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "Stage",value = "state",-LRpairs) %>%
  mutate(Stage=factor(Stage,levels = Stage.levels)) %>%
  group_by(Stage,state) %>%
  summarise(percent=n()) %>%
  mutate(percent = percent/length(LRpairs)) %>%
  ungroup()

p <- ggplot(subset(data.plot,state=="consistent"),aes(Stage,percent,fill=state,label=round(percent,2))) +
  geom_bar(stat="identity")+
  geom_text(vjust = -0.2,position = position_dodge(0.8),size = 4)+
  ylab("Percent")+
  theme_cowplot(font_size = 15)+
  #scale_fill_manual(values = mycolor(9)[ii])+
  scale_y_continuous(limits = c(0,1),
                     breaks = seq(0,1,0.1),
                     expand = c(0,0))+
  scale_x_discrete(expand = c(0.05,0))+
  theme(legend.position = "right",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/consistence_score_LR_210_expresion_levels_state_",suffix = ".jpg",width = 8,height = 8,dpi=300)



####----------5.4 correlation(a)--------------

####----------5.4.1 LR  1619--------------

tmp.mat <- readRDS(file = "res/R/data.plot.mean.mat.rds")

LRpairs <- readRDS(file = "res/R/NATMI_LRpairs_1619LR_20201103.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")


tmp.mat.L <- tmp.mat[Lgene.list,]
tmp.mat.R <- tmp.mat[Rgene.list,]
tmp.LR.pairs <- LRpairs

tmp.res.list <- sapply(Stage.levels,
                       function(x) cor(tmp.mat.L[,x],tmp.mat.R[,x]),
                       USE.NAMES = T)
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))


tmp.res.list.df <- data.frame(Stage=names(tmp.res.list),
                              cor = tmp.res.list,
                              stringsAsFactors = F)
tmp.res.list.df$Stage <- factor(tmp.res.list.df$Stage,levels = Stage.levels)
p <- ggplot(tmp.res.list.df,aes(Stage,cor,fill=Stage))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))

myggsave(p,prefix = "res/fig/LR_1619_correlation",suffix = ".jpg",width = 8,height = 8,dpi=300)

####----------5.4.2 LR 266--------------

tmp.mat <- readRDS(file = "res/R/data.plot.mean.mat.rds")

LRpairs <- readRDS(file = "res/R/peLR_266.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")


tmp.mat.L <- tmp.mat[Lgene.list,]
tmp.mat.R <- tmp.mat[Rgene.list,]
tmp.LR.pairs <- LRpairs

tmp.res.list <- sapply(Stage.levels,
                       function(x) cor(tmp.mat.L[,x],tmp.mat.R[,x]),
                       USE.NAMES = T)
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))


tmp.res.list.df <- data.frame(Stage=names(tmp.res.list),
                              cor = tmp.res.list,
                              stringsAsFactors = F)
tmp.res.list.df$Stage <- factor(tmp.res.list.df$Stage,levels = Stage.levels)
p <- ggplot(tmp.res.list.df,aes(Stage,cor,fill=Stage))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))

myggsave(p,prefix = "res/fig/LR_266_correlation",suffix = ".jpg",width = 8,height = 8,dpi=300)

####----------5.4.3 LR 210--------------

tmp.mat <- readRDS(file = "res/R/data.plot.mean.mat.rds")

LRpairs <- readRDS(file = "res/R/eLR_correlation.Rds") %>%
  dplyr::filter(abs(PCC) > 0.2) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 


Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")


tmp.mat.L <- tmp.mat[Lgene.list,]
tmp.mat.R <- tmp.mat[Rgene.list,]
tmp.LR.pairs <- LRpairs

tmp.res.list <- sapply(Stage.levels,
                       function(x) cor(tmp.mat.L[,x],tmp.mat.R[,x]),
                       USE.NAMES = T)
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))


tmp.res.list.df <- data.frame(Stage=names(tmp.res.list),
                              cor = tmp.res.list,
                              stringsAsFactors = F)
tmp.res.list.df$Stage <- factor(tmp.res.list.df$Stage,levels = Stage.levels)
p <- ggplot(tmp.res.list.df,aes(Stage,cor,fill=Stage))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))

myggsave(p,prefix = "res/fig/LR_210_correlation",suffix = ".jpg",width = 8,height = 8,dpi=300)


####20201201
####----------5.4.4 pos--------------
tmp.mat <- readRDS(file = "res/R/data.plot.mean.mat.rds")

LRpairs <- readRDS(file = "res/R/eLR_correlation.Rds") %>%
  dplyr::filter(PCC > 0.2) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 


Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")


tmp.mat.L <- tmp.mat[Lgene.list,]
tmp.mat.R <- tmp.mat[Rgene.list,]
tmp.LR.pairs <- LRpairs

tmp.res.list <- sapply(Stage.levels,
                       function(x) cor(tmp.mat.L[,x],tmp.mat.R[,x]),
                       USE.NAMES = T)
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))


tmp.res.list.df <- data.frame(Stage=names(tmp.res.list),
                              cor = tmp.res.list,
                              stringsAsFactors = F)
tmp.res.list.df$Stage <- factor(tmp.res.list.df$Stage,levels = Stage.levels)
p <- ggplot(tmp.res.list.df,aes(Stage,cor,fill=Stage))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/LR_210_cor_160_pos_correlation",suffix = ".jpg",width = 8,height = 8,dpi=300)


####----------5.4.4 neg--------------
tmp.mat <- readRDS(file = "res/R/data.plot.mean.mat.rds")

LRpairs <- readRDS(file = "res/R/eLR_correlation.Rds") %>%
  dplyr::filter(PCC <  -0.2) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 


Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")


tmp.mat.L <- tmp.mat[Lgene.list,]
tmp.mat.R <- tmp.mat[Rgene.list,]
tmp.LR.pairs <- LRpairs

tmp.res.list <- sapply(Stage.levels,
                       function(x) cor(tmp.mat.L[,x],tmp.mat.R[,x]),
                       USE.NAMES = T)
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))


tmp.res.list.df <- data.frame(Stage=names(tmp.res.list),
                              cor = tmp.res.list,
                              stringsAsFactors = F)
tmp.res.list.df$Stage <- factor(tmp.res.list.df$Stage,levels = Stage.levels)
p <- ggplot(tmp.res.list.df,aes(Stage,cor,fill=Stage))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/LR_210_cor_50_neg_correlation",suffix = ".jpg",width = 8,height = 8,dpi=300)


####-------5.4 correlation (b) ---------

####-------5.4.1 LR 1619---------

tmp.mat <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
tmp.mat <- myNormalize(myRemoveNA(tmp.mat))
colnames(tmp.mat)[1] <- "Oocyte"
boxplot(tmp.mat)

LRpairs <- readRDS(file = "res/R/NATMI_LRpairs_1619LR_20201103.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early_2cell","2cell",
                   "4cell","8cell","ICM")


tmp.mat.L <- tmp.mat[Lgene.list,Stage.levels]
tmp.mat.R <- tmp.mat[Rgene.list,Stage.levels]
tmp.LR.pairs <- LRpairs

tmp.res.list <- sapply(Stage.levels,
                       function(x) cor(tmp.mat.L[,x],tmp.mat.R[,x]),
                       USE.NAMES = T)
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))


tmp.res.list.df <- data.frame(Stage=names(tmp.res.list),
                              cor = tmp.res.list,
                              stringsAsFactors = F)
tmp.res.list.df$Stage <- factor(tmp.res.list.df$Stage,levels = Stage.levels)
p <- ggplot(tmp.res.list.df,aes(Stage,cor,fill=Stage))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/xiewei_data_LR_1619_correlation",suffix = ".jpg",width = 8,height = 8,dpi=300)

####-------5.4.2 LR 266---------

tmp.mat <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
tmp.mat <- myNormalize(myRemoveNA(tmp.mat))
colnames(tmp.mat)[1] <- "Oocyte"
boxplot(tmp.mat)

LRpairs <- readRDS(file = "res/R/peLR_266.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early_2cell","2cell",
                   "4cell","8cell","ICM")


tmp.mat.L <- tmp.mat[Lgene.list,Stage.levels]
tmp.mat.R <- tmp.mat[Rgene.list,Stage.levels]
tmp.LR.pairs <- LRpairs

tmp.res.list <- sapply(Stage.levels,
                       function(x) cor(tmp.mat.L[,x],tmp.mat.R[,x]),
                       USE.NAMES = T)
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))


tmp.res.list.df <- data.frame(Stage=names(tmp.res.list),
                              cor = tmp.res.list,
                              stringsAsFactors = F)
tmp.res.list.df$Stage <- factor(tmp.res.list.df$Stage,levels = Stage.levels)
p <- ggplot(tmp.res.list.df,aes(Stage,cor,fill=Stage))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/xiewei_data_LR_266_correlation",suffix = ".jpg",width = 8,height = 8,dpi=300)

####-------5.4.3 LR 210---------

tmp.mat <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
tmp.mat <- myNormalize(myRemoveNA(tmp.mat))
colnames(tmp.mat)[1] <- "Oocyte"
boxplot(tmp.mat)

LRpairs <- readRDS(file = "res/R/eLR_correlation.Rds") %>%
  dplyr::filter(abs(PCC) > 0.2) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early_2cell","2cell",
                   "4cell","8cell","ICM")


tmp.mat.L <- tmp.mat[Lgene.list,Stage.levels]
tmp.mat.R <- tmp.mat[Rgene.list,Stage.levels]
tmp.LR.pairs <- LRpairs

tmp.res.list <- sapply(Stage.levels,
                       function(x) cor(tmp.mat.L[,x],tmp.mat.R[,x]),
                       USE.NAMES = T)
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))


tmp.res.list.df <- data.frame(Stage=names(tmp.res.list),
                              cor = tmp.res.list,
                              stringsAsFactors = F)
tmp.res.list.df$Stage <- factor(tmp.res.list.df$Stage,levels = Stage.levels)
p <- ggplot(tmp.res.list.df,aes(Stage,cor,fill=Stage))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/xiewei_data_LR_210_correlation",suffix = ".jpg",width = 8,height = 8,dpi=300)


####-------5.4.4 pos---------

tmp.mat <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
tmp.mat <- myNormalize(myRemoveNA(tmp.mat))
colnames(tmp.mat)[1] <- "Oocyte"
boxplot(tmp.mat)

LRpairs <- readRDS(file = "res/R/eLR_correlation.Rds") %>%
  dplyr::filter(PCC > 0.2) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early_2cell","2cell",
                   "4cell","8cell","ICM")


tmp.mat.L <- tmp.mat[Lgene.list,Stage.levels]
tmp.mat.R <- tmp.mat[Rgene.list,Stage.levels]
tmp.LR.pairs <- LRpairs

tmp.res.list <- sapply(Stage.levels,
                       function(x) cor(tmp.mat.L[,x],tmp.mat.R[,x]),
                       USE.NAMES = T)
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))


tmp.res.list.df <- data.frame(Stage=names(tmp.res.list),
                              cor = tmp.res.list,
                              stringsAsFactors = F)
tmp.res.list.df$Stage <- factor(tmp.res.list.df$Stage,levels = Stage.levels)
p <- ggplot(tmp.res.list.df,aes(Stage,cor,fill=Stage))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/xiewei_data_LR_160_pos_correlation",suffix = ".jpg",width = 8,height = 8,dpi=300)


####-------5.4.4 neg---------

tmp.mat <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
tmp.mat <- myNormalize(myRemoveNA(tmp.mat))
colnames(tmp.mat)[1] <- "Oocyte"
boxplot(tmp.mat)

LRpairs <- readRDS(file = "res/R/eLR_correlation.Rds") %>%
  dplyr::filter(PCC < -0.2) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early_2cell","2cell",
                   "4cell","8cell","ICM")


tmp.mat.L <- tmp.mat[Lgene.list,Stage.levels]
tmp.mat.R <- tmp.mat[Rgene.list,Stage.levels]
tmp.LR.pairs <- LRpairs

tmp.res.list <- sapply(Stage.levels,
                       function(x) cor(tmp.mat.L[,x],tmp.mat.R[,x]),
                       USE.NAMES = T)
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))


tmp.res.list.df <- data.frame(Stage=names(tmp.res.list),
                              cor = tmp.res.list,
                              stringsAsFactors = F)
tmp.res.list.df$Stage <- factor(tmp.res.list.df$Stage,levels = Stage.levels)
p <- ggplot(tmp.res.list.df,aes(Stage,cor,fill=Stage))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/xiewei_data_LR_50_neg_correlation",suffix = ".jpg",width = 8,height = 8,dpi=300)

####-------5.5 bcdcor-----------

####-------5.5.1 load data-------

data.merge.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20201028.rds")
data.merge.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20201028.rds")

####------5.5.2 1619 LR-----------
LRpairs <- readRDS(file = "res/R/NATMI_LRpairs_1619LR_20201103.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")
res <- c()
for(ii in Stage.levels){
  cat("Stage:",ii)
  tmp.stage <- ii
  tmp.cell_id <- data.merge.meta %>%
    filter(Stage == tmp.stage) %>%
    pull(cell_id)
  X <- t(data.merge.exp[Lgene.list,tmp.cell_id])
  Y <- t(data.merge.exp[Rgene.list,tmp.cell_id])
  tmp.res <- bcdcor(X,Y)
  res <- c(res,tmp.res)
}

res.df <- data.frame(Stage=Stage.levels,
                     bcdccor=res,
                     stringsAsFactors = F)


res.df$Stage <- factor(res.df$Stage,levels = Stage.levels)


p <- ggplot(res.df,aes(Stage,bcdccor,fill=Stage))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/LR_1619_bcdcor",suffix = ".jpg",width = 8,height = 8,dpi=300)

####--------5.5.3 210 LR-----------
LRpairs <- readRDS(file = "res/R/eLR_correlation.Rds") %>%
  dplyr::filter(abs(PCC) > 0.2) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.levels <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                   "4cell","8cell","16cell","earlyblast","midblast",  
                   "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")
res <- c()
for(ii in Stage.levels){
  cat("Stage:",ii)
  tmp.stage <- ii
  tmp.cell_id <- data.merge.meta %>%
    filter(Stage == tmp.stage) %>%
    pull(cell_id)
  X <- t(data.merge.exp[Lgene.list,tmp.cell_id])
  Y <- t(data.merge.exp[Rgene.list,tmp.cell_id])
  tmp.res <- bcdcor(X,Y)
  res <- c(res,tmp.res)
}

res.df <- data.frame(Stage=Stage.levels,
                     bcdccor=res,
                     stringsAsFactors = F)


res.df$Stage <- factor(res.df$Stage,levels = Stage.levels)


p <- ggplot(res.df,aes(Stage,bcdccor,fill=Stage))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 15)+
  scale_fill_manual(values = test.color.2(16))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text( size = 15 ) ,
        axis.text = element_text(size = 15),
        axis.line.x = element_line( size = 1 ),
        axis.line.y = element_line( size = 1 ),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
myggsave(p,prefix = "res/fig/eLR_210_bcdcor",suffix = ".jpg",width = 8,height = 8,dpi=300)


####----6.assymetry----------


####----6.1 by pca-------

####----6.1.1 load data----
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20201019.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20201019.rds")

beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)

Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                 "4cell","8cell","16cell","earlyblast","midblast",  
                 "lateblast","fibroblast")

head(beforeEPI.meta)


####------6.1.2 1619 LR--------

LRpairs <- readRDS(file = "res/R/NATMI_LRpairs_1619LR_20201103.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

###----early2cell----

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="early2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="early2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_early2cell_1619_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)

###----mid2cell----

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="mid2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="mid2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_mid2cell_1619_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


###----late2cell----

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="late2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="late2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_late2cell_1619_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


####------6.1.2 266 LR--------

LRpairs <- readRDS(file = "res/R/peLR_266.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

###----early2cell-----

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="early2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="early2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_early2cell_266_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


###----mid2cell----

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="mid2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="mid2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_mid2cell_266_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


###----late2cell----
data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="late2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="late2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_late2cell_266_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


####------6.1.2 210 LR--------
LRpairs <- readRDS(file = "res/R/eLR_160_pos_50_neg_20201130.rds") %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 


####-----early2cell---------
data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="early2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="early2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_early2cell_210_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)

###----mid2cell----

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="mid2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="mid2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_mid2cell_210_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


###----late2cell----
data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="late2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="late2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_late2cell_210_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)

####----6.assymetry----


####----6.1 by pca-------

####----6.1.1 load data----
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20201019.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20201019.rds")

beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)

Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                 "4cell","8cell","16cell","earlyblast","midblast",  
                 "lateblast","fibroblast")

head(beforeEPI.meta)


####------6.1.2 1619 LR--------

LRpairs <- readRDS(file = "res/R/NATMI_LRpairs_1619LR_20201103.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

###----early2cell----

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="early2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="early2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_early2cell_1619_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)

###----mid2cell----

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="mid2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="mid2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_mid2cell_1619_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


###----late2cell----

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="late2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="late2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_late2cell_1619_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


####------6.1.2 266 LR--------

LRpairs <- readRDS(file = "res/R/peLR_266.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

###----early2cell-----

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="early2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="early2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_early2cell_266_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


###----mid2cell----

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="mid2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="mid2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_mid2cell_266_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


###----late2cell----
data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="late2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="late2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_late2cell_266_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


####------6.1.2 210 LR--------
LRpairs <- readRDS(file = "res/R/eLR_160_pos_50_neg_20201130.rds") %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 


####-----early2cell---------
data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="early2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="early2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_early2cell_210_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)

###----mid2cell----

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="mid2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="mid2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_mid2cell_210_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)


###----late2cell----
data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  filter(Stage=="late2cell") %>%
  pull(cell_id)

data.plot.L <- data.plot[Lgene.list,tmp.select]
data.plot.R <- data.plot[Rgene.list,tmp.select]

data.plot <- abs(data.plot.L-data.plot.R)
boxplot(data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####get meta data
tmp.meta <- beforeEPI.meta %>%
  filter(Stage=="late2cell")


data.plot <- pca.res$x %>%
  as.data.frame() %>%
  dplyr::select(PC1,PC2) %>%
  rownames_to_column("cell_id") %>%
  left_join(tmp.meta,by = "cell_id")


p <- ggplot(data.plot,aes(PC1,PC2,color=embryo_id))+
  geom_point(size=5)+
  xlab(paste0("PC1:",ext_labels[1]))+
  ylab(paste0("PC2:",ext_labels[2]))+
  scale_color_manual(values =c("#FCB03B","#0F6937","#67BD45",
                               "#169192","#B09DBA","#954496",
                               "#4E1550","#8EAF3D","#F05153",
                               "#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 18,font_family = "sans")
p
myggsave(p,prefix = "res/fig/Assymetry_late2cell_210_LR",suffix = ".jpg",width = 8,height = 8,dpi=300)




####----6.2 by density-------

####----6.2.1 load data----
data.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20201028.rds")
data.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20201028.rds")


Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                  "4cell","8cell","16cell","earlyblast","midblast",  
                  "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")


####------6.2.2 1619 LR--------

LRpairs <- readRDS(file = "res/R/LRpairs_1619_20201117.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

####------LR------
tmp.res.df <- data.frame(stringsAsFactors = F)
for (ii in Stage.level) {
  data.plot <- data.exp
  tmp.stage <- ii
  tmp.select <- data.meta %>%
    filter(Stage==tmp.stage) %>%
    pull(cell_id)
  data.plot.L <- data.plot[Lgene.list,tmp.select]
  data.plot.R <- data.plot[Rgene.list,tmp.select]
  data.plot <- colMeans(data.plot.L+data.plot.R)
  tmp <- data.frame(Stage=ii,value=data.plot,cell_id=names(data.plot),stringsAsFactors = F)
  tmp.res.df <- rbind(tmp.res.df,tmp)
}

tmp.data.plot <- tmp.res.df

plot.list <- list()
for(ii in Stage.level){
  cat("Stage:",ii,"\n")
  tmp.stage <- ii
  p <- ggplot()+
    geom_density(data = subset(tmp.data.plot,Stage==tmp.stage),
                 aes(value),alpha=0.01)+
    scale_x_continuous(limits = c(0,5),expand = c(0,0.1))+
    #scale_y_continuous(expand = c(0,0.1))+
    ggtitle(tmp.stage)+
    theme_cowplot(font_size = 18)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
  
  plot.list[[ii]] <- p
}
p <- plot_grid(plotlist = plot.list,nrow = 4,ncol = 4)

myggsave(p,prefix = "res/fig/Assymetry_density_1744_LR",suffix = ".jpg",width = 12,height = 16,dpi=300)


####-----background sample 100 times----
tmp.res.list <- list()
for(ii in 1:100){
  cat("simu:",ii,",start\n")
  gene_symbols <- rownames(data.exp)
  tmp.LR.pairs <- NULL
  while (length(tmp.LR.pairs) < length(LRpairs)) {
    tmp.Lgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.Rgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.LR.pairs <- paste0(tmp.Lgene.list,"_",tmp.Rgene.list)
  }
  tmp.res.df <- data.frame(stringsAsFactors = F)
  for (jj in Stage.level) {
    data.plot <- data.exp
    tmp.stage <- jj
    tmp.select <- data.meta %>%
      filter(Stage==tmp.stage) %>%
      pull(cell_id)
    data.plot.L <- data.plot[tmp.Lgene.list,tmp.select]
    data.plot.R <- data.plot[tmp.Rgene.list,tmp.select]
    data.plot <- colMeans(data.plot.L+data.plot.R)
    tmp <- data.frame(Stage=jj,value=data.plot,cell_id=names(data.plot),stringsAsFactors = F)
    tmp.res.df <- rbind(tmp.res.df,tmp)
  }
  tmp.res.df$simu_count <- paste0("simu_",ii)
  tmp.res.list <- c(tmp.res.list,list(tmp.res.df))
  cat("simu:",ii,",end\n")
}

tmp.data.plot <- tmp.res.list %>%
  purrr::reduce(rbind)

plot.list <- list()
for(ii in Stage.level){
  cat("Stage:",ii,"\n")
  tmp.stage <- ii
  p <- ggplot()+
    geom_density(data = subset(tmp.data.plot,Stage==tmp.stage),
                 aes(value,color=simu_count),alpha=0.01)+
    scale_color_manual(values = rep("#dbdbdb",100))+
    geom_density(data = tmp.data.plot %>% 
                   dplyr::filter(Stage==tmp.stage) %>%
                   group_by(cell_id) %>%
                   summarise(tt.res = mean(value)) %>%
                   ungroup(),aes(tt.res))+
    scale_x_continuous(limits = c(0,5),expand = c(0,0.1))+
    ggtitle(tmp.stage)+
    scale_y_continuous(expand = c(0,0.1))+
    theme_cowplot(font_size = 18)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
  
  plot.list[[ii]] <- p
}
p <- plot_grid(plotlist = plot.list,nrow = 4,ncol = 4)
myggsave(p,prefix = "res/fig/Assymetry_density_1744_randomLR_simu100",suffix = ".jpg",width = 12,height = 16,dpi=300)




####-------6.2.3 210 eLR---------
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)

Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                 "4cell","8cell","16cell","earlyblast","midblast",  
                 "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")

####------LR------
tmp.res.df <- data.frame(stringsAsFactors = F)
for (ii in Stage.level) {
  data.plot <- data.exp
  tmp.stage <- ii
  tmp.select <- data.meta %>%
    filter(Stage==tmp.stage) %>%
    pull(cell_id)
  data.plot.L <- data.plot[Lgene.list,tmp.select]
  data.plot.R <- data.plot[Rgene.list,tmp.select]
  data.plot <- colMeans(data.plot.L+data.plot.R)
  tmp <- data.frame(Stage=ii,value=data.plot,cell_id=names(data.plot),stringsAsFactors = F)
  tmp.res.df <- rbind(tmp.res.df,tmp)
}

tmp.data.plot <- tmp.res.df

plot.list <- list()
for(ii in Stage.level){
  cat("Stage:",ii,"\n")
  tmp.stage <- ii
  p <- ggplot()+
    geom_density(data = subset(tmp.data.plot,Stage==tmp.stage),
                 aes(value),alpha=0.01)+
    scale_x_continuous(limits = c(0,7),expand = c(0,0.1))+
    #scale_y_continuous(expand = c(0,0.1))+
    ggtitle(tmp.stage)+
    theme_cowplot(font_size = 18)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
  
  plot.list[[ii]] <- p
}
p <- plot_grid(plotlist = plot.list,nrow = 4,ncol = 4)

myggsave(p,prefix = "res/fig/Assymetry_density_239_eLR",suffix = ".jpg",width = 12,height = 16,dpi=300)


####-----background sample 100 times----
tmp.res.list <- list()
for(ii in 1:100){
  cat("simu:",ii,",start\n")
  gene_symbols <- rownames(data.exp)
  tmp.LR.pairs <- NULL
  while (length(tmp.LR.pairs) < length(LRpairs)) {
    tmp.Lgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.Rgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.LR.pairs <- paste0(tmp.Lgene.list,"_",tmp.Rgene.list)
  }
  tmp.res.df <- data.frame(stringsAsFactors = F)
  for (jj in Stage.level) {
    data.plot <- data.exp
    tmp.stage <- jj
    tmp.select <- data.meta %>%
      filter(Stage==tmp.stage) %>%
      pull(cell_id)
    data.plot.L <- data.plot[tmp.Lgene.list,tmp.select]
    data.plot.R <- data.plot[tmp.Rgene.list,tmp.select]
    data.plot <- colMeans(data.plot.L+data.plot.R)
    tmp <- data.frame(Stage=jj,value=data.plot,cell_id=names(data.plot),stringsAsFactors = F)
    tmp.res.df <- rbind(tmp.res.df,tmp)
  }
  tmp.res.df$simu_count <- paste0("simu_",ii)
  tmp.res.list <- c(tmp.res.list,list(tmp.res.df))
  cat("simu:",ii,",end\n")
}

tmp.data.plot <- tmp.res.list %>%
  purrr::reduce(rbind)

plot.list <- list()
for(ii in Stage.level){
  cat("Stage:",ii,"\n")
  tmp.stage <- ii
  p <- ggplot()+
    geom_density(data = subset(tmp.data.plot,Stage==tmp.stage),
                 aes(value,color=simu_count),alpha=0.01)+
    scale_color_manual(values = rep("#dbdbdb",100))+
    geom_density(data = tmp.data.plot %>% 
                   dplyr::filter(Stage==tmp.stage) %>%
                   group_by(cell_id) %>%
                   summarise(tt.res = mean(value)) %>%
                   ungroup(),aes(tt.res))+
    scale_x_continuous(limits = c(0,7),expand = c(0,0.1))+
    ggtitle(tmp.stage)+
    scale_y_continuous(expand = c(0,0.1))+
    theme_cowplot(font_size = 18)+
    theme(legend.position = "none",plot.title = element_text(hjust = 0.5))
  
  plot.list[[ii]] <- p
}
p <- plot_grid(plotlist = plot.list,nrow = 4,ncol = 4)
myggsave(p,prefix = "res/fig/Assymetry_density_239_randomLR_simu100",suffix = ".jpg",width = 12,height = 16,dpi=300)



####-------xxxx. what is function of eLR---------

ttt <- readRDS("res/R/eLR_correlation.Rds")




#####------7.integrated with atac-seq----------

#####------7.1 atac density in gene------------

#####preprocess the bedgraph file in terminal


###----7.1.1 gene reduced length------
tmp.df <- read.delim(file = "res/atac-seq-process/mm9_gene.bed",header = F,stringsAsFactors = F)
tttt <- tmp.df %>%
  group_by(V4) %>%
  arrange(desc(V5)) %>%
  dplyr::slice(1) %>%
  arrange(V1,V2)
write.table(tttt,file = "res/atac-seq-process/mm9_gene_reduced.bed",quote = F,col.names = F,row.names = F,sep = "\t")

length(unique(tttt$V4)) == nrow(tttt)


####------7.1.2 read bedtools map no rep------
tmp.dir <- "res/atac-seq-process/density/"
tmp.files <- list.files(path = tmp.dir)

tmp.res <- lapply(tmp.files, function(x){
  tmp.df <- read.delim(file = paste0(tmp.dir,x),
                       header = F,
                       stringsAsFactors = F)
  colnames(tmp.df) <- c("chr","start","end","gene","gene.length","strand","density")
  tmp.df$density <- as.numeric(tmp.df$density)
  tmp.df$density <- ifelse(is.na(tmp.df$density),0,tmp.df$density)
  tmp.df.res <- tmp.df %>%
    mutate(stage=gsub(pattern = "\\.bed",replacement="",x)) %>%
    dplyr::select(gene,density,stage)
  return(tmp.df.res)
}) %>%
  purrr::reduce(rbind) %>%
  spread(key = "stage",value = "density") %>%
  column_to_rownames("gene")

saveRDS(tmp.res,file = "res/R/xiewei.data.atac.norep_20201127.rds")


####-----7.1.3 gene promoter atac density correlation------

####------7.1.3.1 load data-------

xiewei.data.atac <- readRDS(file = "res/R/xiewei.data.atac.norep_20201127.rds")
xiewei.data.atac <- myNormalize(xiewei.data.atac)

xiewei.data.exp <- readRDS("res/R/xiewei.data.exp_20201012.rds")
xiewei.data.exp <- myNormalize(xiewei.data.exp)

colnames(xiewei.data.exp)[3] <- "early 2-cell"
colnames(xiewei.data.exp)[4] <- "late 2-cell"
colnames(xiewei.data.atac)[2] <- "early 2-cell"
colnames(xiewei.data.atac)[1] <- "late 2-cell"


####-------7.1.3.2 peLR gene-----------
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)

#LRpairs <- readRDS(file = "res/R/NATMI_LRpairs_20201028.rds")
####
tmp_eLR <- LRpairs
tmp.LR.pairs <- tmp_eLR
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)


###Lgene
tmp.gene.use <- intersect(tmp.gene.symbol,tmp.L.gene)
tmp.atac <- xiewei.data.atac[tmp.gene.use,]
tmp.exp <- xiewei.data.exp[tmp.gene.use,]

colnames(tmp.atac)
colnames(tmp.exp)

tmp.stage <- c("early 2-cell","late 2-cell","4cell","8cell","ICM")
tmp.res <- c()
tmp.res.p <- c()
for (ii in tmp.stage) {
  x <- tmp.atac[,ii]
  y <- tmp.exp[,ii]
  tmp.res <- c(tmp.res,cor(x,y))
  tmp.res.p <- c(tmp.res.p,cor.test(x,y)$p.value)
}

tmp.df <- data.frame(Stage = tmp.stage,
                     PCC = tmp.res,
                     P = tmp.res.p,
                     stringsAsFactors = F)
tmp.df$Stage <- factor(tmp.df$Stage,levels = tmp.stage)
ggplot(data = tmp.df,aes(Stage,PCC))+
  geom_point(shape=15,color="skyblue",size=5)+
  geom_line(group="",size=1,color="blue")+
  geom_text(aes(label=ifelse(P < 1e-2,"**","")),vjust=-1)+
  ggtitle(label = "eLR Ligand gene")+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size  = 18,hjust = 0.8,angle = 30),
        axis.text.y = element_text(size = 18))
ggsave(filename = myFileName(prefix = "res/fig/atac_eLR_L_cor",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)

###Rgene
tmp.gene.use <- intersect(tmp.gene.symbol,tmp.R.gene)
tmp.atac <- xiewei.data.atac[tmp.gene.use,]
tmp.exp <- xiewei.data.exp[tmp.gene.use,]

colnames(tmp.atac)
colnames(tmp.exp)

tmp.stage <- c("early 2-cell","late 2-cell","4cell","8cell","ICM")
tmp.res <- c()
tmp.res.p <- c()
for (ii in tmp.stage) {
  x <- tmp.atac[,ii]
  y <- tmp.exp[,ii]
  tmp.res <- c(tmp.res,cor(x,y))
  tmp.res.p <- c(tmp.res.p,cor.test(x,y)$p.value)
}

tmp.df <- data.frame(Stage = tmp.stage,
                     PCC = tmp.res,
                     P = tmp.res.p,
                     stringsAsFactors = F)
tmp.df$Stage <- factor(tmp.df$Stage,levels = tmp.stage)
ggplot(data = tmp.df,aes(Stage,PCC))+
  geom_point(shape=15,color="skyblue",size=5)+
  geom_line(group="",size=1,color="blue")+
  geom_text(aes(label=ifelse(P < 1e-2,"**","")),vjust=-1)+
  ggtitle(label = "eLR Receptor gene")+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size  = 18,hjust = 0.8,angle = 30),
        axis.text.y = element_text(size = 18))
ggsave(filename = myFileName(prefix = "res/fig/atac_eLR_R_cor",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)




####-------7.1.3.3 correlation of eLR-----------
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff) %>%
  pull(LRpairs)
tmp_eLR <- LRpairs
tmp.LR.pairs <- tmp_eLR
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)

x <- xiewei.data.atac[tmp.L,]
y <- xiewei.data.atac[tmp.R,]
tmp.stage <- c("early 2-cell","late 2-cell","4cell","8cell","ICM")
tmp.res <- c()
tmp.res.p <- c()
for (ii in tmp.stage) {
  xx <- x[,ii]
  yy <- y[,ii]
  tmp.res <- c(tmp.res,cor(xx,yy))
  tmp.res.p <- c(tmp.res.p,cor.test(xx,yy)$p.value)
}
tmp.df <- data.frame(Stage = tmp.stage,
                     PCC = tmp.res,
                     P = tmp.res.p,
                     stringsAsFactors = F)
tmp.df$Stage <- factor(tmp.df$Stage,levels = tmp.stage)
ggplot(data = tmp.df,aes(Stage,PCC))+
  geom_point(shape=15,color="skyblue",size=5)+
  geom_line(group="",size=1,color="blue")+
  geom_text(aes(label=ifelse(P < 1e-2,"**","")),vjust=-1)+
  ggtitle(label = "eLR")+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size  = 18,hjust = 0.8,angle = 30),
        axis.text.y = element_text(size = 18))
ggsave(filename = myFileName(prefix = "res/fig/atac_eLR_cor",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


####------7.1.3.4 correlation of pos eLR----------
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(PCC > cutoff) %>%
  pull(LRpairs)
tmp_eLR <- LRpairs
tmp.LR.pairs <- tmp_eLR
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)

x <- xiewei.data.atac[tmp.L,]
y <- xiewei.data.atac[tmp.R,]
tmp.stage <- c("early 2-cell","late 2-cell","4cell","8cell","ICM")
tmp.res <- c()
tmp.res.p <- c()
for (ii in tmp.stage) {
  xx <- x[,ii]
  yy <- y[,ii]
  tmp.res <- c(tmp.res,cor(xx,yy))
  tmp.res.p <- c(tmp.res.p,cor.test(xx,yy)$p.value)
}
tmp.df <- data.frame(Stage = tmp.stage,
                     PCC = tmp.res,
                     P = tmp.res.p,
                     stringsAsFactors = F)
tmp.df$Stage <- factor(tmp.df$Stage,levels = tmp.stage)
ggplot(data = tmp.df,aes(Stage,PCC))+
  geom_point(shape=15,color="skyblue",size=5)+
  geom_line(group="",size=1,color="blue")+
  geom_text(aes(label=ifelse(P < 1e-2,"**","")),vjust=-1)+
  ggtitle(label = "pos correlated pos eLR")+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size  = 18,hjust = 0.8,angle = 30),
        axis.text.y = element_text(size = 18))
ggsave(filename = myFileName(prefix = "res/fig/atac_eLR_LR_pos_cor",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


###Lgene
tmp.gene.use <- intersect(tmp.gene.symbol,tmp.L.gene)
tmp.atac <- xiewei.data.atac[tmp.gene.use,]
tmp.exp <- xiewei.data.exp[tmp.gene.use,]

colnames(tmp.atac)
colnames(tmp.exp)

tmp.stage <- c("early 2-cell","late 2-cell","4cell","8cell","ICM")
tmp.res <- c()
tmp.res.p <- c()
for (ii in tmp.stage) {
  x <- tmp.atac[,ii]
  y <- tmp.exp[,ii]
  tmp.res <- c(tmp.res,cor(x,y))
  tmp.res.p <- c(tmp.res.p,cor.test(x,y)$p.value)
}

tmp.df <- data.frame(Stage = tmp.stage,
                     PCC = tmp.res,
                     P = tmp.res.p,
                     stringsAsFactors = F)
tmp.df$Stage <- factor(tmp.df$Stage,levels = tmp.stage)
ggplot(data = tmp.df,aes(Stage,PCC))+
  geom_point(shape=15,color="skyblue",size=5)+
  geom_line(group="",size=1,color="blue")+
  geom_text(aes(label=ifelse(P < 1e-2,"**","")),vjust=-1)+
  ggtitle(label = "eLR Ligand gene")+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size  = 18,hjust = 0.8,angle = 30),
        axis.text.y = element_text(size = 18))
ggsave(filename = myFileName(prefix = "res/fig/atac_eLR_pos_L_cor",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)

###Rgene
tmp.gene.use <- intersect(tmp.gene.symbol,tmp.R.gene)
tmp.atac <- xiewei.data.atac[tmp.gene.use,]
tmp.exp <- xiewei.data.exp[tmp.gene.use,]

colnames(tmp.atac)
colnames(tmp.exp)

tmp.stage <- c("early 2-cell","late 2-cell","4cell","8cell","ICM")
tmp.res <- c()
tmp.res.p <- c()
for (ii in tmp.stage) {
  x <- tmp.atac[,ii]
  y <- tmp.exp[,ii]
  tmp.res <- c(tmp.res,cor(x,y))
  tmp.res.p <- c(tmp.res.p,cor.test(x,y)$p.value)
}

tmp.df <- data.frame(Stage = tmp.stage,
                     PCC = tmp.res,
                     P = tmp.res.p,
                     stringsAsFactors = F)
tmp.df$Stage <- factor(tmp.df$Stage,levels = tmp.stage)
ggplot(data = tmp.df,aes(Stage,PCC))+
  geom_point(shape=15,color="skyblue",size=5)+
  geom_line(group="",size=1,color="blue")+
  geom_text(aes(label=ifelse(P < 1e-2,"**","")),vjust=-1)+
  ggtitle(label = "eLR Receptor gene")+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size  = 18,hjust = 0.8,angle = 30),
        axis.text.y = element_text(size = 18))
ggsave(filename = myFileName(prefix = "res/fig/atac_eLR_pos_R_cor",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)



####-----------7.1.3.5 correlation of eLR neg--------------
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(PCC < -cutoff) %>%
  pull(LRpairs)
tmp_eLR <- LRpairs
tmp.LR.pairs <- tmp_eLR
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)

x <- xiewei.data.atac[tmp.L,]
y <- xiewei.data.atac[tmp.R,]
tmp.stage <- c("early 2-cell","late 2-cell","4cell","8cell","ICM")
tmp.res <- c()
tmp.res.p <- c()
for (ii in tmp.stage) {
  xx <- x[,ii]
  yy <- y[,ii]
  tmp.res <- c(tmp.res,cor(xx,yy))
  tmp.res.p <- c(tmp.res.p,cor.test(xx,yy)$p.value)
}


tmp.df <- data.frame(Stage = tmp.stage,
                     PCC = tmp.res,
                     P = tmp.res.p,
                     stringsAsFactors = F)
tmp.df$Stage <- factor(tmp.df$Stage,levels = tmp.stage)
ggplot(data = tmp.df,aes(Stage,PCC))+
  geom_point(shape=15,color="skyblue",size=5)+
  geom_line(group="",size=1,color="blue")+
  geom_text(aes(label=ifelse(P < 1e-2,"**","")),vjust=-1)+
  ggtitle(label = "neg correlated eLR")+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size  = 18,hjust = 0.8,angle = 30),
        axis.text.y = element_text(size = 18))
ggsave(filename = myFileName(prefix = "res/fig/atac_eLR_LR_neg_cor",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)

###Lgene
tmp.gene.use <- intersect(tmp.gene.symbol,tmp.L.gene)
tmp.atac <- xiewei.data.atac[tmp.gene.use,]
tmp.exp <- xiewei.data.exp[tmp.gene.use,]

colnames(tmp.atac)
colnames(tmp.exp)

tmp.stage <- c("early 2-cell","late 2-cell","4cell","8cell","ICM")
tmp.res <- c()
tmp.res.p <- c()
for (ii in tmp.stage) {
  x <- tmp.atac[,ii]
  y <- tmp.exp[,ii]
  tmp.res <- c(tmp.res,cor(x,y))
  tmp.res.p <- c(tmp.res.p,cor.test(x,y)$p.value)
}

tmp.df <- data.frame(Stage = tmp.stage,
                     PCC = tmp.res,
                     P = tmp.res.p,
                     stringsAsFactors = F)
tmp.df$Stage <- factor(tmp.df$Stage,levels = tmp.stage)
ggplot(data = tmp.df,aes(Stage,PCC))+
  geom_point(shape=15,color="skyblue",size=5)+
  geom_line(group="",size=1,color="blue")+
  geom_text(aes(label=ifelse(P < 1e-2,"**","")),vjust=-1)+
  ggtitle(label = "eLR Ligand gene")+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size  = 18,hjust = 0.8,angle = 30),
        axis.text.y = element_text(size = 18))
ggsave(filename = myFileName(prefix = "res/fig/atac_eLR_neg_L_cor",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)

###Rgene
tmp.gene.use <- intersect(tmp.gene.symbol,tmp.R.gene)
tmp.atac <- xiewei.data.atac[tmp.gene.use,]
tmp.exp <- xiewei.data.exp[tmp.gene.use,]

colnames(tmp.atac)
colnames(tmp.exp)

tmp.stage <- c("early 2-cell","late 2-cell","4cell","8cell","ICM")
tmp.res <- c()
tmp.res.p <- c()
for (ii in tmp.stage) {
  x <- tmp.atac[,ii]
  y <- tmp.exp[,ii]
  tmp.res <- c(tmp.res,cor(x,y))
  tmp.res.p <- c(tmp.res.p,cor.test(x,y)$p.value)
}

tmp.df <- data.frame(Stage = tmp.stage,
                     PCC = tmp.res,
                     P = tmp.res.p,
                     stringsAsFactors = F)
tmp.df$Stage <- factor(tmp.df$Stage,levels = tmp.stage)
ggplot(data = tmp.df,aes(Stage,PCC))+
  geom_point(shape=15,color="skyblue",size=5)+
  geom_line(group="",size=1,color="blue")+
  geom_text(aes(label=ifelse(P < 1e-2,"**","")),vjust=-1)+
  ggtitle(label = "eLR Receptor gene")+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(size  = 18,hjust = 0.8,angle = 30),
        axis.text.y = element_text(size = 18))
ggsave(filename = myFileName(prefix = "res/fig/atac_eLR_neg_R_cor",suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)

####---------7.2 Intersect with open gene-------

####---------7.2.1 load data------
tmp.res.list <- list()
tmp.files <- list.files("res/atac-seq-process/open_gene/")
for(ii in tmp.files){
  tmp.res <- readLines(paste0("res/atac-seq-process/open_gene/",ii))
  tmp.res.list <- c(tmp.res.list,list(tmp.res))
}
names(tmp.res.list) <- gsub(tmp.files,pattern = "*.txt",replacement = "")

####eLR
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff) %>%
  pull(LRpairs)
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)

tmp.stage <- c("2cell_early","2cell","4cell","8cell","icm")
tmp.res.list <- sapply(tmp.stage,FUN = function(x) tmp.res.list[[x]],USE.NAMES = T)


###intersect test

tmp.res.L <- sapply(tmp.res.list,function(x) intersect(x,Lgene))
tmp.res.R <- sapply(tmp.res.list,function(x) intersect(x,Rgene))

####------7.3 help track visual--------

xiewei.data.exp <- readRDS("res/R/xiewei.data.exp_20201012.rds")
xiewei.data.exp <- myNormalize(xiewei.data.exp)
xiewei.data.exp <- myRemoveNA(xiewei.data.exp)
xiewei.data.exp <- pheatmap:::scale_mat(xiewei.data.exp,scale = "row")
xiewei.data.exp <- myRemoveNA(xiewei.data.exp)


tmp.stage <- c("early_2cell","2cell","4cell","8cell","ICM")
tmp.mat <- xiewei.data.exp[,tmp.stage]
tmp.gene <- c("Adam12","Sdc4")



tmp.limits <- c(-3,3)
tmp.gene <- c("Adam12","Sdc4")
pdf(file = "res/tmp/Adam12-Sdc4_exp.pdf",width = 8,height = 8)
par(mfrow=c(2,1))
barplot(rep(1,length(tmp.stage)),
        col = map2color(x = tmp.mat[tmp.gene[1],],
                        limits = tmp.limits,
                        pal = rdbu(100)),
        border = NA,
        axes = F)
barplot(rep(1,length(tmp.stage)),
        col = map2color(x = tmp.mat[tmp.gene[2],],
                        limits = tmp.limits,
                        pal = rdbu(100)),
        border = NA,
        axes = F)
dev.off()

tmp.limits <- c(-3,3)
tmp.gene <- c("Pyy","Dpp4")
pdf(file = "res/tmp/Pyy-Dpp4_exp.pdf",width = 8,height = 8)
par(mfrow=c(2,1))
barplot(rep(1,length(tmp.stage)),
        col = map2color(x = tmp.mat[tmp.gene[1],],
                        limits = tmp.limits,
                        pal = rdbu(100)),
        border = NA,
        axes = F)
barplot(rep(1,length(tmp.stage)),
        col = map2color(x = tmp.mat[tmp.gene[2],],
                        limits = tmp.limits,
                        pal = rdbu(100)),
        border = NA,
        axes = F)
dev.off()



#####---------7.4 atac-seq correlation changes ------------
####------7.4.1 load data-------

xiewei.data.atac <- readRDS(file = "res/R/xiewei.data.atac.norep_20201127.rds")
xiewei.data.atac <- myNormalize(xiewei.data.atac)

xiewei.data.exp <- readRDS("res/R/xiewei.data.exp_20201012.rds")
xiewei.data.exp <- myNormalize(xiewei.data.exp)
xiewei.data.exp <- myRemoveNA(xiewei.data.exp)

colnames(xiewei.data.exp)[3:6] <- c("early 2-cell","late 2-cell","4-cell","8-cell")
colnames(xiewei.data.atac)[1:4] <- c("late 2-cell","early 2-cell","4-cell","8-cell")


range(xiewei.data.atac)
range(xiewei.data.exp)
colnames(xiewei.data.atac)
colnames(xiewei.data.exp)
tmp.stage <- c("early 2-cell","late 2-cell","4-cell","8-cell","ICM")


####-----7.4.2 atac-coordination -----------

####-----7.4.2.1 pos eLR -----

tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(PCC > cutoff) %>%
  pull(LRpairs)
tmp_eLR <- LRpairs
tmp.LR.pairs <- tmp_eLR
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)

x <- xiewei.data.atac[tmp.L,]
y <- xiewei.data.atac[tmp.R,]
tmp.stage <- c("early 2-cell","late 2-cell","4-cell","8-cell","ICM")
tmp.res <- c()
tmp.res.p <- c()
for (ii in tmp.stage) {
  xx <- x[,ii]
  yy <- y[,ii]
  tmp.res <- c(tmp.res,cor(xx,yy))
  tmp.res.p <- c(tmp.res.p,cor.test(xx,yy)$p.value)
}
tmp.df <- data.frame(Stage = tmp.stage,
                     PCC = tmp.res,
                     P = tmp.res.p,
                     stringsAsFactors = F)
tmp.df$Stage <- factor(tmp.df$Stage,levels = tmp.stage)
p <- ggplot(data = tmp.df,aes(Stage,PCC))+
  geom_line(group="",size=1,color="blue")+
  geom_point(shape=15,color="skyblue",size=5)+
  geom_text(aes(label=ifelse(P < 1e-2,"**","")),size=6,vjust=-1)+
  ggtitle(label = "pos correlated eLR")+
  scale_y_continuous(expand = c(0,0.1))+
  theme_cowplot()+
  theme(axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5,size = 28),
        axis.text.x = element_text(size  = 28,hjust = 0.8,angle = 30),
        axis.text.y = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28))
p
ggsave(plot = p,
       filename = myFileName(prefix = "res/fig/atac_eLR_LR_pos_cor",suffix = ".jpg"),
       width = 8,height = 8,
       dpi = 350)
####-----7.4.2.2 neg eLR -----

tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(PCC < -cutoff) %>%
  pull(LRpairs)
tmp_eLR <- LRpairs
tmp.LR.pairs <- tmp_eLR
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)

x <- xiewei.data.atac[tmp.L,]
y <- xiewei.data.atac[tmp.R,]
tmp.stage <- c("early 2-cell","late 2-cell","4-cell","8-cell","ICM")
tmp.res <- c()
tmp.res.p <- c()
for (ii in tmp.stage) {
  xx <- x[,ii]
  yy <- y[,ii]
  tmp.res <- c(tmp.res,cor(xx,yy))
  tmp.res.p <- c(tmp.res.p,cor.test(xx,yy)$p.value)
}
tmp.df <- data.frame(Stage = tmp.stage,
                     PCC = tmp.res,
                     P = tmp.res.p,
                     stringsAsFactors = F)
tmp.df$Stage <- factor(tmp.df$Stage,levels = tmp.stage)
p <- ggplot(data = tmp.df,aes(Stage,PCC))+
  geom_line(group="",size=1,color="blue")+
  geom_point(shape=15,color="skyblue",size=5)+
  geom_text(aes(label=ifelse(P < 1e-2,"**","")),size=6,vjust=-1)+
  ggtitle(label = "negative correlated eLR")+
  scale_y_continuous(expand = c(0,0.1))+
  theme_cowplot()+
  theme(axis.line = element_line(size = 1),
        plot.title = element_text(hjust = 0.5,size = 28),
        axis.text.x = element_text(size  = 28,hjust = 0.8,angle = 30),
        axis.text.y = element_text(size = 28),
        axis.title.x = element_text(size = 28),
        axis.title.y = element_text(size = 28))
p
ggsave(plot = p,
       filename = myFileName(prefix = "res/fig/atac_eLR_LR_neg_cor",suffix = ".jpg"),
       width = 8,height = 8,
       dpi = 350)



####-----7.4.3 atac-coordination score -----------



####-----7.4.3.1 pos eLR -----

####------pos eLR--------
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(PCC > cutoff) %>%
  pull(LRpairs)
tmp_eLR <- LRpairs
tmp.LR.pairs <- tmp_eLR
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)



####save to entrez id
gene.list <- tmp.L.gene
eg = bitr(gene.list, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)
gene.list  <- eg$ENTREZID
cat(tmp.L.gene,sep = "\n",file = myFileName(prefix = "res/txt/eLR_L_gene_entrez",suffix = ".txt"))


gene.list <- tmp.R.gene
eg = bitr(gene.list, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)
gene.list  <- eg$ENTREZID
cat(tmp.R.gene,sep = "\n",file = myFileName(prefix = "res/txt/eLR_R_gene_entrez",suffix = ".txt"))


x <- xiewei.data.atac[tmp.L,]
y <- xiewei.data.atac[tmp.R,]
tmp.stage <- c("early 2-cell","late 2-cell","4-cell","8-cell","ICM")
tmp.res <- list()

for (ii in tmp.stage) {
  xx <- x[,ii]
  yy <- y[,ii]
  tmp.res.df <- data.frame(xx*yy*exp(-abs(xx-yy)))
  colnames(tmp.res.df) <- ii
  tmp.res <- c(tmp.res,list(tmp.res.df))
}

tmp.data.plot <- purrr::reduce(tmp.res,cbind) %>%
  gather(key = "Stage",value = "Score") %>%
  mutate(group = "p-eLR")

tmp.data.plot.1 <- tmp.data.plot
###-------generate random eLR----------
gene_symbols <- rownames(xiewei.data.atac)
tmp.res.list <- list()
for (ii in 1:100) {
  cat("simu:",ii,",start\n")
  tmp.LR.pairs <- NULL
  while (length(tmp.LR.pairs) < length(LRpairs)) {
    tmp.Lgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.Rgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.LR.pairs <- paste0(tmp.Lgene.list,"_",tmp.Rgene.list)
  }
  tmp.res <- list()
  
  x <- xiewei.data.atac[tmp.Lgene.list,]
  y <- xiewei.data.atac[tmp.Rgene.list,]
  tmp.stage <- c("early 2-cell","late 2-cell","4-cell","8-cell","ICM")
  
  for (jj in tmp.stage) {
    xx <- x[,jj]
    yy <- y[,jj]
    tmp.res.df <- data.frame(xx*yy*exp(-abs(xx-yy)))
    colnames(tmp.res.df) <- jj
    tmp.res <- c(tmp.res,list(tmp.res.df))
  }
  
  tmp.res.df <- purrr::reduce(tmp.res ,cbind) %>%
    gather(key = "Stage",value = "Score") %>%
    mutate(group = "randomLR",simu_count = ii)
  ### save data frame to each list
  tmp.res.list <- c(tmp.res.list,list(tmp.res.df))
  cat("simu:",ii,",end\n")
}

tmp.data.plot.2 <- tmp.res.list %>%
  purrr::reduce(rbind) %>%
  dplyr::select(-simu_count)

tmp.data.plot <- rbind(tmp.data.plot.1,tmp.data.plot.2)
tmp.stage <- c("early 2-cell","late 2-cell","4-cell","8-cell","ICM")
tmp.data.plot$Stage <- factor(tmp.data.plot$Stage,levels = tmp.stage )
tmp.data.plot$group <- factor(tmp.data.plot$group,levels = c("p-eLR","randomLR"))


#####-------plot----------
mycolor <- c("#73c0dc","#f2f2f2")

p <- ggplot(tmp.data.plot,aes(Stage,Score,fill=group,color=group))+
  geom_boxplot(linetype = "dashed", 
               outlier.shape = NA,
               width = 0.8 ,
               position = position_dodge(0.9)) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),
               outlier.shape=NA,
               width = 0.8 ,
               position = position_dodge(0.9)) +
  stat_boxplot(geom = "errorbar", 
               aes(ymin = ..ymax..),
               width = 0.3 ,
               position = position_dodge(0.9)) +
  stat_boxplot(geom = "errorbar", 
               aes(ymax = ..ymin..),
               width = 0.3,
               position = position_dodge(0.9)) +
  stat_boxplot(data = tmp.data.plot,
               geom = "errorbar", 
               aes(x = Stage,
                   y = Score,
                   ymin = ..middle..,
                   ymax = ..middle..),
               size = 1,
               col="black",
               width = 0.8 ,
               position = position_dodge(0.9))+
  theme_cowplot(font_size = 28)+
  ylab("Open Chromatin Coordination Score")+
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  scale_x_discrete(labels=tmp.stage)+
  scale_y_continuous(limits = c(0,75),
                     breaks = seq(0,75,by = 25))+
  theme(axis.line = element_line(size = 1),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
ggsave(plot = p,
       filename = myFileName(prefix = "res/fig/atac_eLR_pos_OCCS",suffix = ".jpg"),
       width = 8,height = 8,
       dpi = 350)

####-----7.4.3.2 neg eLR -----

####------neg eLR--------
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(PCC < -cutoff) %>%
  pull(LRpairs)
tmp_eLR <- LRpairs
tmp.LR.pairs <- tmp_eLR
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)

x <- xiewei.data.atac[tmp.L,]
y <- xiewei.data.atac[tmp.R,]
tmp.stage <- c("early 2-cell","late 2-cell","4-cell","8-cell","ICM")
tmp.res <- list()

for (ii in tmp.stage) {
  xx <- x[,ii]
  yy <- y[,ii]
  tmp.res.df <- data.frame(xx*yy*exp(-abs(xx-yy)))
  colnames(tmp.res.df) <- ii
  tmp.res <- c(tmp.res,list(tmp.res.df))
}

tmp.data.plot <- purrr::reduce(tmp.res,cbind)
boxplot(tmp.data.plot)
tmp.data.plot <- purrr::reduce(tmp.res,cbind) %>%
  gather(key = "Stage",value = "Score") %>%
  mutate(group = "n-eLR")
tmp.data.plot.1 <- tmp.data.plot
###-------generate random eLR----------
gene_symbols <- rownames(xiewei.data.atac)
tmp.res.list <- list()
for (ii in 1:100) {
  cat("simu:",ii,",start\n")
  tmp.LR.pairs <- NULL
  while (length(tmp.LR.pairs) < length(LRpairs)) {
    tmp.Lgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.Rgene.list <- sample(gene_symbols,length(LRpairs),replace = T)
    tmp.LR.pairs <- paste0(tmp.Lgene.list,"_",tmp.Rgene.list)
  }
  tmp.res <- list()
  
  x <- xiewei.data.atac[tmp.Lgene.list,]
  y <- xiewei.data.atac[tmp.Rgene.list,]
  tmp.stage <- c("early 2-cell","late 2-cell","4-cell","8-cell","ICM")
  
  for (jj in tmp.stage) {
    xx <- x[,jj]
    yy <- y[,jj]
    tmp.res.df <- data.frame(xx*yy*exp(-abs(xx-yy)))
    colnames(tmp.res.df) <- jj
    tmp.res <- c(tmp.res,list(tmp.res.df))
  }
  
  tmp.res.df <- purrr::reduce(tmp.res ,cbind) %>%
    gather(key = "Stage",value = "Score") %>%
    mutate(group = "randomLR",simu_count = ii)
  ### save data frame to each list
  tmp.res.list <- c(tmp.res.list,list(tmp.res.df))
  cat("simu:",ii,",end\n")
}

tmp.data.plot.2 <- tmp.res.list %>%
  purrr::reduce(rbind) %>%
  dplyr::select(-simu_count)

tmp.data.plot <- rbind(tmp.data.plot.1,tmp.data.plot.2)
tmp.stage <- c("early 2-cell","late 2-cell","4-cell","8-cell","ICM")
tmp.data.plot$Stage <- factor(tmp.data.plot$Stage,levels = tmp.stage )
tmp.data.plot$group <- factor(tmp.data.plot$group,levels = c("n-eLR","randomLR"))


#####-------plot----------
mycolor <- c("#73c0dc","#f2f2f2")

p <- ggplot(tmp.data.plot,aes(Stage,Score,fill=group,color=group))+
  geom_boxplot(linetype = "dashed", 
               outlier.shape = NA,
               width = 0.8 ,
               position = position_dodge(0.9)) +
  stat_boxplot(aes(ymin = ..lower.., ymax = ..upper..),
               outlier.shape=NA,
               width = 0.8 ,
               position = position_dodge(0.9)) +
  stat_boxplot(geom = "errorbar", 
               aes(ymin = ..ymax..),
               width = 0.3 ,
               position = position_dodge(0.9)) +
  stat_boxplot(geom = "errorbar", 
               aes(ymax = ..ymin..),
               width = 0.3,
               position = position_dodge(0.9)) +
  stat_boxplot(data = tmp.data.plot,
               geom = "errorbar", 
               aes(x = Stage,
                   y = Score,
                   ymin = ..middle..,
                   ymax = ..middle..),
               size = 1,
               col="black",
               width = 0.8 ,
               position = position_dodge(0.9))+
  theme_cowplot(font_size = 28)+
  ylab("Open Chromatin Coordination Score")+
  scale_fill_manual(values = mycolor)+
  scale_color_manual(values = mycolor)+
  scale_x_discrete(labels=tmp.stage)+
  scale_y_continuous(limits = c(0,75),
                     breaks = seq(0,75,by = 25))+
  theme(axis.line = element_line(size = 1),
        axis.text.x = element_text(hjust = 0.8,angle = 30))
p
ggsave(plot = p,
       filename = myFileName(prefix = "res/fig/atac_eLR_neg_OCCS",suffix = ".jpg"),
       width = 8,height = 8,
       dpi = 350)


####---------7.5 TF network by ATAC-seq-----------


####---------7.5.0 load datasets------

####------load eLR--------
p_eLR <- readRDS(file = "res/R/p_eLR_20201230.rds")
n_eLR <- readRDS(file = "res/R/n_eLR_20201230.rds")


####-----load exp data-------

data.mean.plot.mat <- readRDS(file = "res/R/data.plot.mean.mat_20201214.rds")




####------load motif and tss data-----

data("motifAnnotations_mgi")
motifRankings <- importRankings("E://wlt/project/data/cistarget/mm9-tss-centered-10kb-10species.mc9nr.feather")


####---------7.5.1 p_eLR -----

geneLists <- list()
for (ii in 1:3) {
  tmp.LR.pairs <- p_eLR %>%
    filter(Cluster == ii) %>%
    pull(LRpairs)
  tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
  tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
  tmp.L.gene  <- unique(tmp.L)
  tmp.R.gene  <- unique(tmp.R)
  geneList1 <- c(tmp.L.gene,tmp.R.gene)
  geneLists <- c(geneLists,list(geneList1))
}
names(geneLists) <- paste0("p-eLR_",1:3)

####-------7.5.1.2 run cisTarget-----------
motifEnrichmentTable_wGenes <- cisTarget(geneLists, 
                                         motifRankings = motifRankings,
                                         motifAnnot = motifAnnotations_mgi,
                                         nesThreshold = 0, 
                                         geneErnMethod="aprox")






####-----7.5.3 plot ---


####prepare vertices
tmp.res <- motifEnrichmentTable_wGenes %>% 
  filter(NES >= 3,TF_highConf != "") %>%
  dplyr::select(-motif, -TF_lowConf )
tmp.res$TF_highConf <- gsub(" \\(.*\\). ", ";", tmp.res$TF_highConf)
data.plot <- tmp.res %>%
  filter(geneSet == "p-eLR_1") %>%
  dplyr::select(source = TF_highConf,target = enrichedGenes,NES)
data.plot$source <- strsplit(data.plot$source,split = ";")
data.plot$target <- strsplit(data.plot$target,split = ";")
data.plot <- apply(data.plot, 1, 
                   FUN = function(x){
                     res <- data.frame(source = rep(x$source,each = length(x$target)),
                                       target = x$target,
                                       NES = x$NES,
                                       stringsAsFactors = F)
                     return(res)
                   })
data.plot <- data.plot %>%
  purrr::reduce(rbind)
data.plot$source <- gsub(data.plot$source,pattern = " ",replacement = "")

tmp.LR.pairs <- p_eLR %>%
  filter(Cluster == 1) %>%
  pull(LRpairs)
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)


###nodes
data.plot.vertice <- data.frame(Gene=c(tmp.L.gene,
                                       tmp.R.gene,
                                       unique(data.plot$source)),
                                Type=c(rep("Ligand",length(tmp.L.gene)),
                                       rep("Receptor",length(tmp.R.gene)),
                                       rep("TF",length(unique(data.plot$source)))),
                                stringsAsFactors = F)

tmp.exp <- data.mean.plot.mat[data.plot.vertice$Gene,]
tmp.exp <- myRemoveNA(tmp.exp)
tmp.exp <- tmp.exp %>%
  rownames_to_column("Gene")
tmp.exp$Gene <- data.plot.vertice$Gene

data.plot.vertice <- merge(data.plot.vertice,tmp.exp,by="Gene")

###edages
data.plot.adj <- data.plot %>%
  dplyr::select(-NES) %>%
  mutate(type="regulation")

tmp.LR.network <- data.frame(source=tmp.L,
                             target=tmp.R,
                             type = "communication",
                             stringsAsFactors = F)

data.plot.adj <- rbind(data.plot.adj,tmp.LR.network)


write.table(data.plot.adj,sep = "\t",file = myFileName(prefix = "res/txt/peLR_cluster1_edges",
                                            suffix = ".txt"),row.names = F,quote = F)
write.table(data.plot.vertice,sep = "\t",file = myFileName(prefix = "res/txt/peLR_cluster1_nodess",
                                                       suffix = ".txt"),row.names = F,quote = F)

#### work in batch
for(ii in 1:3){
  ####prepare vertices
  tmp.res <- motifEnrichmentTable_wGenes %>% 
    filter(NES >= 3,TF_highConf != "") %>%
    dplyr::select(-motif, -TF_lowConf )
  tmp.res$TF_highConf <- gsub(" \\(.*\\). ", ";", tmp.res$TF_highConf)
  data.plot <- tmp.res %>%
    filter(geneSet == paste0("p-eLR_",ii)) %>%
    dplyr::select(source = TF_highConf,target = enrichedGenes,NES)
  data.plot$source <- strsplit(data.plot$source,split = ";")
  data.plot$target <- strsplit(data.plot$target,split = ";")
  data.plot <- apply(data.plot, 1, 
                     FUN = function(x){
                       res <- data.frame(source = rep(x$source,each = length(x$target)),
                                         target = x$target,
                                         NES = x$NES,
                                         stringsAsFactors = F)
                       return(res)
                     })
  data.plot <- data.plot %>%
    purrr::reduce(rbind)
  data.plot$source <- gsub(data.plot$source,pattern = " ",replacement = "")
  
  tmp.LR.pairs <- p_eLR %>%
    filter(Cluster == ii) %>%
    pull(LRpairs)
  tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
  tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
  tmp.L.gene  <- unique(tmp.L)
  tmp.R.gene  <- unique(tmp.R)
  
  
  ###nodes
  data.plot.vertice <- data.frame(Gene=c(tmp.L.gene,
                                         tmp.R.gene,
                                         unique(data.plot$source)),
                                  Type=c(rep("Ligand",length(tmp.L.gene)),
                                         rep("Receptor",length(tmp.R.gene)),
                                         rep("TF",length(unique(data.plot$source)))),
                                  stringsAsFactors = F)
  
  tmp.exp <- data.mean.plot.mat[data.plot.vertice$Gene,]
  tmp.exp <- myRemoveNA(tmp.exp)
  tmp.exp <- tmp.exp %>%
    rownames_to_column("Gene")
  tmp.exp$Gene <- data.plot.vertice$Gene
  
  data.plot.vertice <- merge(data.plot.vertice,tmp.exp,by="Gene")
  
  ###edages
  data.plot.adj <- data.plot %>%
    dplyr::select(-NES) %>%
    mutate(type="regulation")
  
  tmp.LR.network <- data.frame(source=tmp.L,
                               target=tmp.R,
                               type = "communication",
                               stringsAsFactors = F)
  
  data.plot.adj <- rbind(data.plot.adj,tmp.LR.network)
  
  
  write.table(data.plot.adj,sep = "\t",file = myFileName(prefix = paste0("res/txt/peLR_cluster",ii,"_edges"),
                                                         suffix = ".txt"),row.names = F,quote = F)
  write.table(data.plot.vertice,sep = "\t",file = myFileName(prefix = paste0("res/txt/peLR_cluster",ii,"_nodes"),
                                                             suffix = ".txt"),row.names = F,quote = F)
  
}


####---------7.5.2 n_eLR -----

geneLists <- list()
for (ii in 1:3) {
  tmp.LR.pairs <- n_eLR %>%
    filter(Cluster == ii) %>%
    pull(LRpairs)
  tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
  tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
  tmp.L.gene  <- unique(tmp.L)
  tmp.R.gene  <- unique(tmp.R)
  geneList1 <- c(tmp.L.gene,tmp.R.gene)
  geneLists <- c(geneLists,list(geneList1))
}
names(geneLists) <- paste0("n-eLR_",1:3)

####-------7.5.2 run cisTarget-----------
motifEnrichmentTable_wGenes <- cisTarget(geneLists, 
                                         motifRankings = motifRankings,
                                         motifAnnot = motifAnnotations_mgi,
                                         nesThreshold = 0, 
                                         geneErnMethod="aprox")






####-----7.5.3 plot ---


####prepare vertices
tmp.res <- motifEnrichmentTable_wGenes %>% 
  filter(NES >= 3,TF_highConf != "") %>%
  dplyr::select(-motif, -TF_lowConf )
tmp.res$TF_highConf <- gsub(" \\(.*\\). ", ";", tmp.res$TF_highConf)

#### work in batch
for(ii in 1:3){
  ####prepare vertices
  tmp.res <- motifEnrichmentTable_wGenes %>% 
    filter(NES >= 3,TF_highConf != "") %>%
    dplyr::select(-motif, -TF_lowConf )
  tmp.res$TF_highConf <- gsub(" \\(.*\\). ", ";", tmp.res$TF_highConf)
  data.plot <- tmp.res %>%
    filter(geneSet == paste0("n-eLR_",ii)) %>%
    dplyr::select(source = TF_highConf,target = enrichedGenes,NES)
  data.plot$source <- strsplit(data.plot$source,split = ";")
  data.plot$target <- strsplit(data.plot$target,split = ";")
  data.plot <- apply(data.plot, 1, 
                     FUN = function(x){
                       res <- data.frame(source = rep(x$source,each = length(x$target)),
                                         target = x$target,
                                         NES = x$NES,
                                         stringsAsFactors = F)
                       return(res)
                     })
  data.plot <- data.plot %>%
    purrr::reduce(rbind)
  data.plot$source <- gsub(data.plot$source,pattern = " ",replacement = "")
  
  tmp.LR.pairs <- n_eLR %>%
    filter(Cluster == ii) %>%
    pull(LRpairs)
  tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
  tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
  tmp.L.gene  <- unique(tmp.L)
  tmp.R.gene  <- unique(tmp.R)
  
  
  ###nodes
  data.plot.vertice <- data.frame(Gene=c(tmp.L.gene,
                                         tmp.R.gene,
                                         unique(data.plot$source)),
                                  Type=c(rep("Ligand",length(tmp.L.gene)),
                                         rep("Receptor",length(tmp.R.gene)),
                                         rep("TF",length(unique(data.plot$source)))),
                                  stringsAsFactors = F)
  
  tmp.exp <- data.mean.plot.mat[data.plot.vertice$Gene,]
  tmp.exp <- myRemoveNA(tmp.exp)
  tmp.exp <- tmp.exp %>%
    rownames_to_column("Gene")
  tmp.exp$Gene <- data.plot.vertice$Gene
  
  data.plot.vertice <- merge(data.plot.vertice,tmp.exp,by="Gene")
  
  ###edages
  data.plot.adj <- data.plot %>%
    dplyr::select(-NES) %>%
    mutate(type="regulation")
  
  tmp.LR.network <- data.frame(source=tmp.L,
                               target=tmp.R,
                               type = "communication",
                               stringsAsFactors = F)
  
  data.plot.adj <- rbind(data.plot.adj,tmp.LR.network)
  
  
  write.table(data.plot.adj,sep = "\t",file = myFileName(prefix = paste0("res/txt/neLR_cluster",ii,"_edges"),
                                                         suffix = ".txt"),row.names = F,quote = F)
  write.table(data.plot.vertice,sep = "\t",file = myFileName(prefix = paste0("res/txt/neLR_cluster",ii,"_nodes"),
                                                             suffix = ".txt"),row.names = F,quote = F)
  
}


### using cytoscape for beautilayout

####------8. eLR function annotation-----------

####--------8.1 total eLR--------
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff) %>%
  pull(LRpairs)
tmp_eLR <- LRpairs
tmp.LR.pairs <- tmp_eLR
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)

cat(tmp.L.gene,sep = "\n",file = myFileName(prefix = "res/txt/eLR_L_gene",suffix = ".txt"))
cat(tmp.R.gene,sep = "\n",file = myFileName(prefix = "res/txt/eLR_R_gene",suffix = ".txt"))


tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
cutoff <- 0.2
LRpairs <- tmp.res.cor.1 %>%
  dplyr::filter(abs(PCC) > cutoff)
write.csv(LRpairs,file = "res/txt/eLR_LR_info.csv",quote = F,row.names = F)


###----------8.1.1 overlap with all candiate eLR-----------
cutoff <- 0.2
###p-elR
tmp.LR.pairs <- tmp.res.cor.1 %>%
  dplyr::filter(PCC > cutoff) %>%
  pull(LRpairs)
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)
p_eLR_gene <- c(tmp.L.gene,tmp.R.gene)
###n-eLR
tmp.LR.pairs <- tmp.res.cor.1 %>%
  dplyr::filter(PCC < -cutoff) %>%
  pull(LRpairs)
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)
n_eLR_gene <- c(tmp.L.gene,tmp.R.gene)


reprogramming_factor_candidate <- read.delim(file = "database/genes_365_for_iPSCs_202001303_DBTMEE.tsv",
                                             sep = "\t",stringsAsFactors = F) %>%
  pull("Gene")
zga_gene_zhangyi <- read_delim(file = "database/other_ZGA/Nfya_KD_zhangyi.txt",delim = "\t") %>%
  pull(gene)
genelist <- list(p_eLR_gene=p_eLR_gene,
                 n_eLR_gene=n_eLR_gene,
                 reprogramming_factor=reprogramming_factor_candidate,
                 zga_gene=zga_gene_zhangyi)
###draw venn diagram
p <- venn.diagram(genelist,
                  filename = NULL,
                  height = 8,
                  width = 8,
                  units = "in",
                  col = "transparent",
                  #category.names = c("p-eLR","n-eLR","R factor"),
                  fill = c("#84c683","#fec04e","#cba6d1","#74a5cc"),
                  cex = 2,
                  cat.cex = 2)


pdf(myFileName(prefix = paste0("res/fig/venn_reprogramming_factor_all_eLR"),
               suffix = ".pdf"),
    width = 8,height = 8)
grid.draw(p)
dev.off()


intersect(genelist$p_eLR_gene,genelist$reprogramming_factor)


#####---------8.2 compare to reprogramming factor----------------

reprogramming_factor_candidate <- read.delim(file = "database/genes_365_for_iPSCs_202001303_DBTMEE.tsv",
                                   sep = "\t",stringsAsFactors = F) %>%
  pull("Gene")

zga_gene_zhangyi <- read_delim(file = "database/other_ZGA/Nfya_KD_zhangyi.txt",delim = "\t") %>%
  pull(gene)


tmp.LR.pairs <- p_eLR$LRpairs
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)

p_eLR_gene <- c(tmp.L.gene,tmp.R.gene)

tmp.LR.pairs <- n_eLR$LRpairs
tmp.L <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[1])) 
tmp.R <- unlist(lapply(strsplit(tmp.LR.pairs,split = "_"),FUN = function(x) x[2]))
tmp.L.gene  <- unique(tmp.L)
tmp.R.gene  <- unique(tmp.R)

n_eLR_gene <- c(tmp.L.gene,tmp.R.gene)

genelist <- list(p_eLR_gene=p_eLR_gene,
                 n_eLR_gene=n_eLR_gene,
                 reprogramming_factor_candidate=reprogramming_factor_candidate,
                 zga_gene=zga_gene_zhangyi)

###draw venn diagram
p <- venn.diagram(genelist,
             filename = NULL,
             height = 8,
             width = 8,
             units = "in",
             col = "transparent",
             #category.names = c("p-eLR","n-eLR","R factor"),
             fill = c("#84c683","#fec04e","#cba6d1","#74a5cc"),
             cex = 2,
             cat.cex = 2)
dev.off()
grid.draw(p)

pdf(myFileName(prefix = paste0("res/fig/venn_reprogramming_factor"),
               suffix = ".pdf"),
    width = 8,height = 8)
grid.draw(p)
dev.off()


tmp.gene.list <- intersect(genelist$p_eLR_gene,genelist$reprogramming_factor_candidate)
cat(tmp.gene.list,sep = "\n",file = myFileName(prefix = "res/txt/reprogramming_factor_overlap",
                                               suffix = ".txt"))



####-----8.3 test database--------

####test escape


tmp.txt <- read.delim(file = "database/ESCAPE/proteomicsESC.txt/proteomicsESC.txt",stringsAsFactors = F)

proteomicsESC <- tmp.txt %>%
  pull(geneName) %>%
  stringr::str_to_title()

intersect(proteomicsESC,genelist$p_eLR_gene)
intersect(proteomicsESC,genelist$n_eLR_gene)



###-------8.4 visualize gene expression---------

###load data
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20201019.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20201019.rds")
beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)
Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                 "4cell","8cell","16cell","earlyblast","midblast",  
                 "lateblast","fibroblast")

LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20201214.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

Stage.level <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                 "4cell","8cell","16cell","earlyblast","midblast","lateblast")

tmp.cell_id <- beforeEPI.meta %>%
  filter(Stage %in% Stage.level) %>%
  pull(cell_id)

### corplot
tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20201215.Rds")
hist(tmp.res.cor.1$PCC)
cor.test(rnorm(100),rnorm(100))
test.color <- colorRampPalette(RColorBrewer::brewer.pal(name = "OrRd",n=9))
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))
Stage.level <-  c("Oocyte","zygote","early2cell",
                  "mid2cell","late2cell","4cell",
                  "8cell","16cell","earlyblast",
                  "midblast","lateblast")
paper.color <- c("#FCB03B","#0F6937","#67BD45",
                 "#169192","#B09DBA","#954496",
                 "#4E1550","#8EAF3D","#F05153",
                 "#5358A5","#79A7AA")
scales::show_col(paper.color)
myCorplot <- function(mat=beforeEPI.exp[,tmp.cell_id],
                      mat.meta=beforeEPI.meta,
                      tmp.Lgene=tmp.res.cor.1$Lgene[1],
                      tmp.Rgene=tmp.res.cor.1$Rgene[1],
                      PCC=tmp.res.cor.1$PCC[1],
                      Stage.level=Stage.level,
                      color=test.color(length(Stage.level))){
  
  mat.L <- mat[tmp.Lgene,]
  mat.R <- mat[tmp.Rgene,]
  
  res.df <- data.frame(t(mat.L),
                       t(mat.R),
                       stringsAsFactors = F) %>%
    rownames_to_column("cell_id") %>%
    left_join(mat.meta,by = "cell_id") %>%
    mutate(Stage=factor(Stage,levels = Stage.level))
  
  p<- ggplot(data = res.df,
             aes_string(x=tmp.Lgene,y=tmp.Rgene,color="Stage"))+
    geom_point(size=8)+
    scale_color_manual(values = color)+
    ggtitle(paste0(tmp.Lgene,"_",tmp.Rgene,",cor:",PCC))+
    theme_cowplot(font_size = 15)+
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text( size = 15 ) ,
          axis.text = element_text(size = 15),
          axis.line.x = element_line( size = 1 ),
          axis.line.y = element_line( size = 1 ))
  return(p)
}

tmp.idx <- which(tmp.res.cor.1$Lgene %in% "Fgf4")

plot.list <- list()
for(ii in tmp.idx){
  p <- myCorplot(mat=beforeEPI.exp[,tmp.cell_id],
                 mat.meta=beforeEPI.meta,
                 tmp.Lgene=tmp.res.cor.1$Lgene[ii],
                 tmp.Rgene=tmp.res.cor.1$Rgene[ii],
                 PCC=tmp.res.cor.1$PCC[ii],
                 Stage.level=Stage.level,
                 color=test.color.2(length(Stage.level)))
  plot.list <- c(plot.list,list(p))
}
plot_grid(plotlist = plot.list)





####scatter plot

data.merge.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20210114.rds")
data.merge.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20210114.rds")

Stage.level <-  c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                  "4cell","8cell","16cell","earlyblast","midblast",  
                  "lateblast","E5.25","E5.5","E6.25","E6.5","fibroblast")
tmp.stage <- c("MII oocyte","zygote",
               "early 2-cell","mid 2-cell","late 2-cell",
               "4-cell","8-cell","16-cell",
               "early blastocyst","mid blastocyst","late blastocyst",
               "E5.25","E5.5","E6.25","E6.5","fibroblast")

tmp.cell_id <- data.merge.meta %>%
  filter(Stage %in% Stage.level) %>%
  pull(cell_id)


myScatterPlot <- function(mat = data.merge.exp[,tmp.cell_id],
                          mat.meta = data.merge.meta,
                          tmp.Lgene = tmp.data.plot$Lgene[nrow(tmp.data.plot)],
                          tmp.Rgene = tmp.data.plot$Rgene[nrow(tmp.data.plot)],
                          PCC = tmp.data.plot$PCC[nrow(tmp.data.plot)]){
  mat.L <- mat[tmp.Lgene,]
  mat.R <- mat[tmp.Rgene,]
  
  res.df <- data.frame(t(mat.L),
                       t(mat.R),
                       stringsAsFactors = F) %>%
    rownames_to_column("cell_id") %>%
    left_join(mat.meta,by = "cell_id") %>%
    mutate(Stage=factor(Stage,levels = Stage.level)) %>%
    gather(key = "gene",value = "gene.exp",-cell_id,-Stage) %>%
    mutate(gene=factor(gene,levels = c(tmp.Lgene,tmp.Rgene)))
  
  
  
  p <- ggplot(data = res.df,aes_string(x="Stage",y="gene.exp",colour="gene"))+
    geom_point(alpha=0.1,size=1)+
    geom_smooth(method = "gam",aes_string(group="gene"),se = F)+
    scale_color_manual(values = c("#E41A1C","#377EB8"))+
    ylab("gene.exp")+
    ggtitle(paste0(tmp.Lgene,"_",tmp.Rgene,",cor:",PCC))+
    theme_cowplot(font_size = 18)+
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          legend.text = element_text( size = 18 ) ,
          axis.text = element_text(size = 18),
          axis.line.x = element_line( size = 1 ),
          axis.line.y = element_line( size = 1 ),
          axis.text.x = element_text(hjust = 0.8,angle = 30))
  return(p)
}
#####structure!!
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = "Ctcf",
                   tmp.Rgene = "Gtf3c2",
                   PCC = "Ctcf")+
  ggtitle("Ctcf")+
  scale_x_discrete(label=tmp.stage)

p

tmp.idx <- 681
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = tmp.res.cor.1$PCC[tmp.idx])+
  scale_x_discrete(label=tmp.stage)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


tmp.idx <- 955
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = tmp.res.cor.1$PCC[tmp.idx])+
  scale_x_discrete(label=tmp.stage)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


tmp.idx <- 292
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = tmp.res.cor.1$PCC[tmp.idx])+
  scale_x_discrete(label=tmp.stage)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)



tmp.idx <- 534
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = tmp.res.cor.1$PCC[tmp.idx])+
  scale_x_discrete(label=tmp.stage)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


tmp.idx <- 27
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = tmp.res.cor.1$PCC[tmp.idx])+
  scale_x_discrete(label=tmp.stage)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


tmp.idx <- 19
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = tmp.res.cor.1$PCC[tmp.idx])+
  scale_x_discrete(label=tmp.stage)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


tmp.idx <- 356
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = tmp.res.cor.1$PCC[tmp.idx])+
  scale_x_discrete(label=tmp.stage)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


tmp.idx <- 109
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = tmp.res.cor.1$PCC[tmp.idx])+
  scale_x_discrete(label=tmp.stage)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


tmp.idx <- 1038
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = tmp.res.cor.1$PCC[tmp.idx])+
  scale_x_discrete(label=tmp.stage)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)


tmp.idx <- 1038
p <- myScatterPlot(mat = data.merge.exp[,tmp.cell_id],
                   mat.meta = data.merge.meta,
                   tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                   tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                   PCC = tmp.res.cor.1$PCC[tmp.idx])+
  scale_x_discrete(label=tmp.stage)

p
ggsave(p,filename = myFileName(prefix = paste0("res/fig/eLR_",tmp.res.cor.1$LRpairs[tmp.idx],"_scatter"),suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)




####--------8.5 enrichment by Msigdb datasets-------------
####read eLR

p_eLR <- readRDS(file = "res/R/p_eLR_20201230.rds")
n_eLR <- readRDS(file = "res/R/n_eLR_20210113.rds")

####using mysidb C5 for enrichment
Mm_msigdb <- msigdbr(species = "Mus musculus")

Mm_C5 <-  Mm_msigdb %>%
  dplyr::filter(gs_cat == "C5" & gs_subcat == "GO:BP") %>%
  dplyr::select(gs_name,gs_exact_source,entrez_gene) 

Mm_C5_GO <- Mm_C5[,1:2] %>%
  unique()


####-------p-eLR------
ii <- 1
tmp_eLR <- p_eLR %>%
  dplyr::filter(Cluster == ii) %>%
  pull(LRpairs)
eLR_L <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[1])) 
eLR_R <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[2]))
eLR_gene <- unique(c(eLR_L,eLR_R))


gene.list <- eLR_gene
eg = bitr(gene.list, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)

em <- enricher(eg$ENTREZID, TERM2GENE=Mm_C5[,c(1,3)],pvalueCutoff = -1,qvalueCutoff = -1)
em <- setReadable(em, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
tmp.res.df <- em@result %>%
  dplyr::filter(p.adjust < 0.05 ) %>%
  arrange(p.adjust) %>%
  rownames_to_column("rownames.id")
em.res <- merge(tmp.res.df,Mm_C5_GO,by.x="ID",by.y="gs_name")
em.res <- em.res %>%
  arrange(p.adjust)


p <- myenrichr_GOplot(df = em.res,
                    fill.color = "#C6DBEFFF",
                    show_number = 20,
                    title = paste0("p-eLR cluster ",ii),
                    font.size = 18,
                    p.adjust.cut = 0.05,
                    plot.ylab = NULL,
                    term.pos.adjust = 0)+
  theme(axis.line.y = element_blank())
p

ggsave(p,
       filename = myFileName(prefix = paste0("res/fig/p_eLR_cluster_",ii,"_GO"),suffix = ".jpg"),
       dpi = 250,width = 8,height = 8)


ii <- 2
tmp_eLR <- p_eLR %>%
  dplyr::filter(Cluster == ii) %>%
  pull(LRpairs)
eLR_L <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[1])) 
eLR_R <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[2]))
eLR_gene <- unique(c(eLR_L,eLR_R))


gene.list <- eLR_gene
eg = bitr(gene.list, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)

em <- enricher(eg$ENTREZID, TERM2GENE=Mm_C5[,c(1,3)],pvalueCutoff = -1,qvalueCutoff = -1)
em <- setReadable(em, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
tmp.res.df <- em@result %>%
  dplyr::filter(p.adjust < 0.05 ) %>%
  arrange(p.adjust) %>%
  rownames_to_column("rownames.id")
em.res <- merge(tmp.res.df,Mm_C5_GO,by.x="ID",by.y="gs_name")
em.res <- em.res %>%
  arrange(p.adjust)


p <- myenrichr_GOplot(df = em.res,
                      fill.color = "#C6DBEFFF",
                      show_number = 20,
                      title = paste0("p-eLR cluster ",ii),
                      font.size = 18,
                      p.adjust.cut = 0.05,
                      plot.ylab = NULL,
                      term.pos.adjust = 0)+
  theme(axis.line.y = element_blank())
p

ggsave(p,
       filename = myFileName(prefix = paste0("res/fig/p_eLR_cluster_",ii,"_GO"),suffix = ".jpg"),
       dpi = 250,width = 8,height = 8)



ii <- 3
tmp_eLR <- p_eLR %>%
  dplyr::filter(Cluster == ii) %>%
  pull(LRpairs)
eLR_L <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[1])) 
eLR_R <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[2]))
eLR_gene <- unique(c(eLR_L,eLR_R))


gene.list <- eLR_gene
eg = bitr(gene.list, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)

em <- enricher(eg$ENTREZID, TERM2GENE=Mm_C5[,c(1,3)],pvalueCutoff = -1,qvalueCutoff = -1)
em <- setReadable(em, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
tmp.res.df <- em@result %>%
  dplyr::filter(p.adjust < 0.05 ) %>%
  arrange(p.adjust) %>%
  rownames_to_column("rownames.id")
em.res <- merge(tmp.res.df,Mm_C5_GO,by.x="ID",by.y="gs_name")
em.res <- em.res %>%
  arrange(p.adjust)


p <- myenrichr_GOplot(df = em.res,
                      fill.color = "#C6DBEFFF",
                      show_number = 20,
                      title = paste0("p-eLR cluster ",ii),
                      font.size = 18,
                      p.adjust.cut = 0.05,
                      plot.ylab = NULL,
                      term.pos.adjust = 0)+
  theme(axis.line.y = element_blank())
p

ggsave(p,
       filename = myFileName(prefix = paste0("res/fig/p_eLR_cluster_",ii,"_GO"),suffix = ".jpg"),
       dpi = 250,width = 8,height = 8)



####-----------n-eLR---------------------------
ii <- 1
tmp_eLR <- n_eLR %>%
  dplyr::filter(Cluster == ii) %>%
  pull(LRpairs)
eLR_L <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[1])) 
eLR_R <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[2]))
eLR_gene <- unique(c(eLR_L,eLR_R))


gene.list <- eLR_gene
eg = bitr(gene.list, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)

em <- enricher(eg$ENTREZID, TERM2GENE=Mm_C5[,c(1,3)],pvalueCutoff = -1,qvalueCutoff = -1)
em <- setReadable(em, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
tmp.res.df <- em@result %>%
  dplyr::filter(p.adjust < 0.05 ) %>%
  arrange(p.adjust) %>%
  rownames_to_column("rownames.id")
em.res <- merge(tmp.res.df,Mm_C5_GO,by.x="ID",by.y="gs_name")
em.res <- em.res %>%
  arrange(p.adjust)


p <- myenrichr_GOplot(df = em.res,
                      fill.color = "#C6DBEFFF",
                      show_number = 20,
                      title = paste0("n-eLR cluster ",ii),
                      font.size = 18,
                      p.adjust.cut = 0.05,
                      plot.ylab = NULL,
                      term.pos.adjust = 0)+
  theme(axis.line.y = element_blank())
p

ggsave(p,
       filename = myFileName(prefix = paste0("res/fig/n_eLR_cluster_",ii,"_GO"),suffix = ".jpg"),
       dpi = 250,width = 8,height = 8)




ii <- 2
tmp_eLR <- n_eLR %>%
  dplyr::filter(Cluster == ii) %>%
  pull(LRpairs)
eLR_L <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[1])) 
eLR_R <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[2]))
eLR_gene <- unique(c(eLR_L,eLR_R))


gene.list <- eLR_gene
eg = bitr(gene.list, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)

em <- enricher(eg$ENTREZID, TERM2GENE=Mm_C5[,c(1,3)],pvalueCutoff = -1,qvalueCutoff = -1)
em <- setReadable(em, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
tmp.res.df <- em@result %>%
  dplyr::filter(p.adjust < 0.05 ) %>%
  arrange(p.adjust) %>%
  rownames_to_column("rownames.id")
em.res <- merge(tmp.res.df,Mm_C5_GO,by.x="ID",by.y="gs_name")
em.res <- em.res %>%
  arrange(p.adjust)


p <- myenrichr_GOplot(df = em.res,
                      fill.color = "#C6DBEFFF",
                      show_number = 20,
                      title = paste0("n-eLR cluster ",ii),
                      font.size = 18,
                      p.adjust.cut = 0.05,
                      plot.ylab = NULL,
                      term.pos.adjust = 0)+
  theme(axis.line.y = element_blank())
p

ggsave(p,
       filename = myFileName(prefix = paste0("res/fig/n_eLR_cluster_",ii,"_GO"),suffix = ".jpg"),
       dpi = 250,width = 8,height = 8)


ii <- 3
tmp_eLR <- n_eLR %>%
  dplyr::filter(Cluster == ii) %>%
  pull(LRpairs)
eLR_L <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[1])) 
eLR_R <- unlist(lapply(strsplit(tmp_eLR,split = "_"),FUN = function(x) x[2]))
eLR_gene <- unique(c(eLR_L,eLR_R))


gene.list <- eLR_gene
eg = bitr(gene.list, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)

em <- enricher(eg$ENTREZID, TERM2GENE=Mm_C5[,c(1,3)],pvalueCutoff = -1,qvalueCutoff = -1)
em <- setReadable(em, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
tmp.res.df <- em@result %>%
  dplyr::filter(p.adjust < 0.05 ) %>%
  arrange(p.adjust) %>%
  rownames_to_column("rownames.id")
em.res <- merge(tmp.res.df,Mm_C5_GO,by.x="ID",by.y="gs_name")
em.res <- em.res %>%
  arrange(p.adjust)


p <- myenrichr_GOplot(df = em.res,
                      fill.color = "#C6DBEFFF",
                      show_number = 20,
                      title = paste0("n-eLR cluster ",ii),
                      font.size = 18,
                      p.adjust.cut = 0.05,
                      plot.ylab = NULL,
                      term.pos.adjust = 0)+
  theme(axis.line.y = element_blank())
p

ggsave(p,
       filename = myFileName(prefix = paste0("res/fig/n_eLR_cluster_",ii,"_GO"),suffix = ".jpg"),
       dpi = 250,width = 8,height = 8)

######--------8.6 category eLR----------------

tmp.df.1 <- read.delim(file = "res/txt/p-eLR-top30_anno_20210113.txt.txt")
tmp.df.2 <- read.delim(file = "res/txt/n-eLR-top30_anno_20210113.txt.txt")
tmp.df <- rbind(tmp.df.1,tmp.df.2) %>%
  group_by(Comment) %>%
  summarise(count=n())
ggplot(tmp.df,aes(Comment,count))+
  geom_bar(stat = "identity")+
  xlab(label = "Category")+
  geom_text(aes(label=count),vjust = -0.2,size = 5)+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        axis.text.x = element_text(angle = 15))
ggsave(filename = "res/fig/eLR_category_stat_top30_20210113.jpg",width = 8,height = 8)

#####-------8.7 eLR screen model------



### generate x
x <- seq(from=0,to=1,length.out = 1000)


###pos
simu.df <- data.frame(Stage=rep(x,times=2),
                      gene.exp=c(cospi(x-1.1)+rnorm(1000,sd = 0.1),
                                 cospi(x-1)+rnorm(1000,sd = 0.1)),
                      gene=rep(c("L","R"),each = 1000),
                      stringsAsFactors = F)

### copy and edit in ppt
dev.new()

ggplot(data = simu.df,aes_string(x="Stage",y="gene.exp",colour="gene"))+
  geom_point(alpha=0,size=1)+
  geom_smooth(method = "gam",aes_string(group="gene"),se = F)+
  scale_color_manual(values = c("#E41A1C","#377EB8"))+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())



###neg
simu.df <- data.frame(Stage=rep(x,times=2),
                      gene.exp=c(cospi(x)+rnorm(1000,sd = 0.1),
                                 cospi(x-1)+rnorm(1000,sd = 0.1)),
                      gene=rep(c("L","R"),each = 1000),
                      stringsAsFactors = F)

### copy and edit in ppt
dev.new()
ggplot(data = simu.df,aes_string(x="Stage",y="gene.exp",colour="gene"))+
  geom_point(alpha=0,size=1)+
  geom_smooth(method = "gam",aes_string(group="gene"),se = F)+
  scale_color_manual(values = c("#E41A1C","#377EB8"))+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())


####debug
x <- seq(from=0,to=1,length.out = 1000)

simu.df <- data.frame(Stage=rep(x,times=2),
                      gene.exp=c(cospi(x-1.5)+rnorm(1000,sd = 0.1),
                                 cospi(x-1)+rnorm(1000,sd = 0.1)),
                      gene=rep(c("L","R"),each = 1000),
                      stringsAsFactors = F)

ggplot(data = simu.df,aes_string(x="Stage",y="gene.exp",colour="gene"))+
  geom_point(alpha=0.1,size=1)+
  geom_smooth(method = "gam",aes_string(group="gene"),se = F)+
  scale_color_manual(values = c("#E41A1C","#377EB8"))+
  theme_cowplot(font_size = 18)+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank())


myIScor <- function(x,y,k=20){
  res.1 <- gtools::running(x,y,fun = cor,width = k,align = "left")
  flag <- sum(res.1>0) > floor(1/2*length(res.1))
  res <- ifelse(flag,max(res.1),min(res.1))
  return(res)
}

x <- seq(from=0,to=1,length.out = 1000)
x <- cospi(x-1.5)+rnorm(1000,sd = 0.1)
y <- cospi(x-1)+rnorm(1000,sd = 0.1)

myIScor(x,y,k = 20)

k = 20
res.1 <- gtools::running(x,y,fun = cor,width = k,align = "left")
plot(res.1)



#######--------9. nichenet ---------

# #######-------9.1 load nichenet data--------
# ligand_target_matrix = readRDS("res/R/ligand_target_matrix.rds")
# ligand_target_matrix[1:5,1:5] # target genes in rows, ligands in columns
# ##                 CXCL1        CXCL2        CXCL3        CXCL5
# ## A1BG     3.534343e-04 4.041324e-04 3.729920e-04 3.080640e-04
# ## A1BG-AS1 1.650894e-04 1.509213e-04 1.583594e-04 1.317253e-04
# ## A1CF     5.787175e-04 4.596295e-04 3.895907e-04 3.293275e-04
# ## A2M      6.027058e-04 5.996617e-04 5.164365e-04 4.517236e-04
# ## A2M-AS1  8.898724e-05 8.243341e-05 7.484018e-05 4.912514e-05
# ##                  PPBP
# ## A1BG     2.628388e-04
# ## A1BG-AS1 1.231819e-04
# ## A1CF     3.211944e-04
# ## A2M      4.590521e-04
# ## A2M-AS1  5.120439e-05
# ligand_tf_matrix = readRDS("res/R/ligand_tf_matrix.rds")
# ligand_tf_matrix[1:5,1:5]
# # CXCL1 CXCL2 CXCL3 CXCL5 PPBP
# # A1BG         0     0     0     0    0
# # A1BG-AS1     0     0     0     0    0
# # A1CF         0     0     0     0    0
# # A2M          0     0     0     0    0
# # A2M-AS1      0     0     0     0    0
# lr_network = readRDS("res/R/lr_network.rds")
# head(lr_network)
# ## # A tibble: 6 x 4
# ##   from  to    source         database
# ##   <chr> <chr> <chr>          <chr>
# ## 1 CXCL1 CXCR2 kegg_cytokines kegg
# ## 2 CXCL2 CXCR2 kegg_cytokines kegg
# ## 3 CXCL3 CXCR2 kegg_cytokines kegg
# ## 4 CXCL5 CXCR2 kegg_cytokines kegg
# ## 5 PPBP  CXCR2 kegg_cytokines kegg
# ## 6 CXCL6 CXCR2 kegg_cytokines kegg
# 
# weighted_networks = readRDS("res/R/weighted_networks.rds")
# weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
# weighted_networks_gr = weighted_networks$gr %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))
# 
# head(weighted_networks$lr_sig) # interactions and their weights in the ligand-receptor + signaling network
# ## # A tibble: 6 x 3
# ##   from  to     weight
# ##   <chr> <chr>   <dbl>
# ## 1 A1BG  ABCC6  0.422
# ## 2 A1BG  ACE2   0.101
# ## 3 A1BG  ADAM10 0.0970
# ## 4 A1BG  AGO1   0.0525
# ## 5 A1BG  AKT1   0.0855
# ## 6 A1BG  ANXA7  0.457
# head(weighted_networks$gr) # interactions and their weights in the gene regulatory network
# ## # A tibble: 6 x 3
# ##   from  to     weight
# ##   <chr> <chr>   <dbl>
# ## 1 A1BG  A2M    0.0294
# ## 2 AAAS  GFAP   0.0290
# ## 3 AADAC CYP3A4 0.0422
# ## 4 AADAC IRF8   0.0275
# ## 5 AATF  ATM    0.0330
# ## 6 AATF  ATR    0.0355
# 
# lr_network = lr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
# colnames(ligand_target_matrix) = ligand_target_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
# rownames(ligand_target_matrix) = ligand_target_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
# ligand_target_matrix = ligand_target_matrix %>% .[!is.na(rownames(ligand_target_matrix)), !is.na(colnames(ligand_target_matrix))]
# 
# colnames(ligand_tf_matrix) = ligand_tf_matrix %>% colnames() %>% convert_human_to_mouse_symbols()
# rownames(ligand_tf_matrix) = ligand_tf_matrix %>% rownames() %>% convert_human_to_mouse_symbols()
# ligand_tf_matrix = ligand_tf_matrix %>% .[!is.na(rownames(ligand_tf_matrix)), !is.na(colnames(ligand_tf_matrix))]
# 
# 
# weighted_networks_lr = weighted_networks_lr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
# weighted_networks_gr = weighted_networks_gr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
# weighted_networks <- list(weighted_networks_lr,weighted_networks_gr)
# names(weighted_networks) <- c("lr_sig","gr")
# 
# 
# saveRDS(object = ligand_tf_matrix,file = "res/R/nichenet_ligand_tf_matrix.rds")
# saveRDS(object = lr_network,file = "res/R/nichenet_lr_network_mouse.rds")
# saveRDS(object = ligand_target_matrix,file = "res/R/nichenet_ligand_target_matrix.rds")
# saveRDS(object = weighted_networks_lr,file = "res/R/nichenet_weighted_network_lr_mouse.rds")
# saveRDS(object = weighted_networks,file = "res/R/nichenet_weighted_networks.rds")

# data("sig_network")
# mouse_sig_network <- sig_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
# saveRDS(mouse_sig_network,file = "res/R/nichenet_mouse_sig_network.rds")
# 
# data("gr_network")
# mouse_gr_network <- gr_network %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
# saveRDS(mouse_gr_network,file = "res/R/nichenet_mouse_gr_network.rds")
# 
# weighted_networks <- readRDS(file = "res/R/weighted_networks.rds")
# weighted_networks$lr_sig <- weighted_networks$lr_sig %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
# weighted_networks$gr <- weighted_networks$gr %>% mutate(from = convert_human_to_mouse_symbols(from), to = convert_human_to_mouse_symbols(to)) %>% drop_na()
# saveRDS(weighted_networks,file = "res/R/nichenet_mouse_weighted_network_20210131.rds")


######---------9.2 run nichenet analysis-------------

######---------9.2.1 load data------------
early.scRNAseq.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20210114.rds")
early.scRNAseq.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20210114.rds")

ligand_target_matrix = readRDS("res/R/nichenet_ligand_target_matrix.rds")
lr_network <- readRDS(file = "res/R/nichenet_lr_network_mouse.rds")
mouse_sig_network <- readRDS(file = "res/R/nichenet_mouse_sig_network.rds")
weighted_networks <- readRDS(file = "res/R/nichenet_mouse_weighted_network_20210131.rds")
ligand_tf_matrix = readRDS(file = "res/R/nichenet_ligand_tf_matrix.rds")

#####---------9.2.2 define expressed genes---------
tmp.id <- early.scRNAseq.meta %>%
  filter(Stage == "early2cell") %>%
  pull(cell_id)
expressed_genes_sender <- early.scRNAseq.exp[,tmp.id] %>%
  apply(2,function(x){2^x-1}) %>% 
  apply(1,function(x){log2( mean(x) +1)}) %>% 
  .[. >= 4] %>% 
  names()
expressed_genes_receiver <- expressed_genes_sender



######-------9.2.3 define gene sets, ligands------
zga_gene_zhangyi <- read_delim(file = "database/other_ZGA/Nfya_KD_zhangyi.txt",delim = "\t") %>%
  pull(gene) 
geneset_oi <- zga_gene_zhangyi
background_expressed_genes = expressed_genes_receiver  %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()



####3d genome list 

architecture_protein <- c(paste0("Gtf3c",1:6),"Ctcf","Yy1","Zbtb33","Parp1","Maz","Prdm5","Rad21",
                          "Ncapd2","Ncapd3","Ncapg","Ncapg2","Ncaph","Ncaph2","Smc2","Smc4",
                          paste0("L3mbtl",1:3),"Chd4","Adnp","Cbx5","Smc1a","Stag1","Wapl","Nipbl",
                          "Hnrnpu","Hnrnpk","Safb","Matr3","Pou5f1","Sox2","Klf4","Myc",
                          "Gata1","Myod","Hmgb2","Mecp2")

##architecture_protein <- c("Ctcf","Rad21","Smc1a","Smc2","Smc4")


geneset_oi <- architecture_protein
background_expressed_genes = expressed_genes_receiver  %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()



#### enrichment
names(ligand_target_matrix["Ctcf",]>0)

tmp.target <- rownames(ligand_target_matrix)[ligand_target_matrix[,"Tgfb2"] > 0]
tmp.zga <- zga_gene_zhangyi
tmp.expressed <- early.scRNAseq.exp[,tmp.id] %>%
  apply(2,function(x){2^x-1}) %>% 
  apply(1,function(x){log2( mean(x) +1)}) %>% 
  .[. >= 0] %>% 
  names()

tmp.1 <- intersect(tmp.expressed,tmp.target)
tmp.2 <- setdiff(tmp.expressed,tmp.target)

sum(zga_gene_zhangyi %in% tttt) / length(tttt)
sum(zga_gene_zhangyi %in% background_expressed_genes) / length(background_expressed_genes)


tmp.test <- matrix(c(length(intersect(tmp.1,tmp.zga)),
              length(intersect(tmp.2,tmp.zga)), 
              length(setdiff(tmp.1,tmp.zga)), 
              length(setdiff(tmp.2,tmp.zga))),
              nrow = 2,byrow = T,
              dimnames = list(ZGA = c("ZGA", "Not ZGA"),
                         Exprressed_gene = c("Target", "Not Target")))
tmp.test
fisher.test(tmp.test, alternative = "greater")




########--------9.2.4 Perform NicheNet’s ligand activity analysis on the gene set of interest--------
ligand_activities = predict_ligand_activities(geneset = geneset_oi, 
                                              background_expressed_genes = background_expressed_genes, 
                                              ligand_target_matrix = ligand_target_matrix, 
                                              potential_ligands = potential_ligands)

ligand_activities %>% arrange(-pearson) 
best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand)
best_upstream_ligands
# show histogram of ligand activity scores
p_hist_lig_activity = ggplot(ligand_activities, aes(x=pearson)) + 
  geom_histogram(color="black", fill="darkorange")  + 
  # geom_density(alpha=.1, fill="orange") +
  geom_vline(aes(xintercept=min(ligand_activities %>% top_n(20, pearson) %>% pull(pearson))), color="red", linetype="dashed", size=1) + 
  labs(x="ligand activity (PCC)", y = "# ligands") +
  theme_classic()
p_hist_lig_activity


######------------9.2.5 Infer target genes of top-ranked ligands and visualize in a heatmap ----------

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 250) %>% bind_rows()

nrow(active_ligand_target_links_df)
head(active_ligand_target_links_df)

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0)
nrow(active_ligand_target_links)
## [1] 143
head(active_ligand_target_links)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets = active_ligand_target_links_df$target %>% unique()
vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()

p_ligand_target_network = vis_ligand_target %>% 
  make_heatmap_ggplot("Prioritized ligands",
                      "Architecture Factor", 
                      color = "purple",
                      x_axis = T,
                      legend_position = "top", 
                      x_axis_position = "top",
                      legend_title = "Regulatory potential") + 
  scale_fill_gradient2(low = "whitesmoke",  high = "purple") + 
  theme(axis.text = element_text(size = 12))

p_ligand_target_network
ggsave(filename = "res/fig/test_architecture_2021.jpg",width = 8,height = 8)



#####-----------9.2.6 : Ligand-receptor network inference for top-ranked ligands------------------------


# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

# convert to a matrix
lr_network_top_df = lr_network_top_df %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% 
  dplyr::select(-to) %>% 
  as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

# perform hierarchical clustering to order the ligands and receptors
dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]

dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
p_ligand_receptor_network = vis_ligand_receptor_network %>% t() %>% make_heatmap_ggplot("Prioritized CAF-ligands","Receptors expressed by malignant cells", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network


ligand_pearson_matrix = ligand_activities %>% 
  dplyr::select(pearson) %>% 
  as.matrix() %>% 
  magrittr::set_rownames(ligand_activities$test_ligand)

vis_ligand_pearson = ligand_pearson_matrix[order_ligands, ] %>% as.matrix(ncol = 1) %>% magrittr::set_colnames("Pearson")
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized CAF-ligands","Ligand activity", color = "darkorange",legend_position = "top", x_axis_position = "top", legend_title = "Pearson correlation coefficient\ntarget gene prediction ability)")
p_ligand_pearson

(p_ligand_pearson | p_ligand_target_network) / p_ligand_receptor_network
ggsave(filename = "res/fig/nichenet_res_test_early2cell_20210130.jpg",width = 8,height = 8,dpi = 350)

plot_grid(p_ligand_pearson+ylab("prioritized Ligands"),
          p_ligand_target_network+ylab(NULL)+
            theme(axis.text.y = element_blank(),
                  axis.ticks.y = element_blank()),rel_widths = c(1,8))
ggsave(filename = "res/fig/nichenet_res_test_early2cell_20210130.jpg",width = 16,height = 9,dpi = 350)


######----------9.2.7 get graph------
ligands_all = c("Tgfb2") # this can be a list of multiple ligands if required
targets_all = c("Icam1","Sdc4")


head(lr_network)
weighted_networks = construct_weighted_networks(lr_network, mouse_sig_network, mouse_gr_network,
                                                source_weights_df)

active_signaling_network = get_ligand_signaling_path(ligand_tf_matrix = ligand_tf_matrix, 
                                                     ligands_all = ligands_all, 
                                                     targets_all = targets_all, 
                                                     weighted_networks = weighted_networks)
# For better visualization of edge weigths: normalize edge weights to make them comparable between signaling and gene regulatory interactions
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 5.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 5.75)
graph_min_max = diagrammer_format_signaling_graph(signaling_graph_list = active_signaling_network, 
                                                  ligands_all = ligands_all, 
                                                  targets_all = targets_all, 
                                                  sig_color = "indianred", 
                                                  gr_color = "steelblue")
unique(active_signaling_network_min_max$sig)
# To render the graph: uncomment following line of code
test <- DiagrammeR::render_graph(graph_min_max, layout = "tree",output = "graph")
test
export_svg(test) %>%
  charToRaw %>% rsvg %>% ('res/fig/test_nichenet_path_width.png')
test.tttt <- export_svg(test) %>%
  charToRaw


test <- graph_min_max %>% to_igraph() 

ttt.layout = layout_with_sugiyama(test)

tmp.vertext <- vertex.attributes(test)
tmp.edge <- edge.attributes(test)
tmp.edge
tmp.layout <- ttt.layout$layout
jpeg(filename = "res/fig/test_nichenet_path_width.jpg",width = 8,height = 8,res = 350,units = "in")
plot.igraph(test,
            layout=tmp.layout,
            vertex.size = 25,
            vertex.label.font=3,
            vertex.label.cex=1,
            vertex.frame.color=tmp.vertext$fillcolor,
            vertex.color=tmp.vertext$fillcolor,
            vertex.label.color=tmp.vertext$fontcolor,
            edge.width = 10*tmp.edge$penwidth)
dev.off()




?plot.igraph



#####----------9.2.8 stat CTCF---------------


sum(ligand_target_matrix["Ctcf",]>0)
ncol(ligand_target_matrix)
####load data
data.plot <- readRDS(file = "res/R/data_merge_plot_20210115.rds")

#####Background
data.plot.mean.L <- data.plot %>%
  filter(group=="Lgene") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp),
            gene.exp.sd=sd(gene.exp),
            gene.exp.cv=ifelse(mean(gene.exp)!=0,sd(gene.exp)/mean(gene.exp),0)) %>%
  ungroup()


tmp.ligand <- colnames(ligand_target_matrix)[ligand_target_matrix["Ctcf",]>0]

tmp.eLR <- unique(tmp.res.LR$Lgene)

tmp.res.1 <- data.plot.mean.L %>%
  dplyr::select(-gene.exp.sd,-gene.exp.cv) %>%
  filter(gene.exp.mean > 0.01) %>%
  filter(Gene %in% tmp.eLR) %>%
  mutate( regulate_ctcf = Gene %in% tmp.ligand ) %>%
  filter(regulate_ctcf == TRUE) 

tmp.res.2 <- data.plot.mean.L %>%
  dplyr::select(-gene.exp.sd,-gene.exp.cv) %>%
  filter(gene.exp.mean > 0.01) %>%
  filter(Gene %in% tmp.eLR) %>%
  mutate( regulate_ctcf = Gene %in% tmp.ligand )


tmp.df <- data.frame(table(tmp.res.1$Stage),
                     table(tmp.res.2$Stage),
                     stringsAsFactors = F) %>%
  dplyr::select(-Var1.1) %>%
  mutate(ratio = Freq/Freq.1)

colnames(tmp.df)[1] <- "Stage"
tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst",
                  "E5.25","E5.5","E6.25","E6.5","fibroblast")
ggplot(tmp.df,aes(Stage,ratio))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=round(ratio,2)),vjust = -0.2,size = 5)+
  scale_x_discrete(labels=tmp.levels.2)+
  theme_cowplot(font_size = 18)+
  theme(axis.text.x = element_text(angle = 30,hjust = 0.8))
ggsave(filename = "res/fig/test_ctcf_filter_gene.jpg",width = 8,height = 8,dpi = 350)



tmp.res.1 <- data.plot.mean.L %>%
  dplyr::select(-gene.exp.sd,-gene.exp.cv) %>%
  filter(gene.exp.mean > 0) %>%
  mutate( regulate_ctcf = Gene %in% tmp.ligand ) %>%
  filter(regulate_ctcf == TRUE) 

tmp.res.2 <- data.plot.mean.L %>%
  dplyr::select(-gene.exp.sd,-gene.exp.cv) %>%
  filter(gene.exp.mean > 0) %>%
  mutate( regulate_ctcf = Gene %in% tmp.ligand )


tmp.df <- data.frame(table(tmp.res.1$Stage),
                     table(tmp.res.2$Stage),
                     stringsAsFactors = F) %>%
  dplyr::select(-Var1.1) %>%
  mutate(ratio = Freq/Freq.1)

colnames(tmp.df)[1] <- "Stage"
tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst",
                  "E5.25","E5.5","E6.25","E6.5","fibroblast")
ggplot(tmp.df,aes(Stage,ratio))+
  geom_bar(stat = "identity")+
  geom_text(aes(label=round(ratio,2)),vjust = -0.2,size = 5)+
  scale_x_discrete(labels=tmp.levels.2)+
  theme_cowplot(font_size = 18)+
  theme(axis.text.x = element_text(angle = 30,hjust = 0.8))
ggsave(filename = "res/fig/test_ctcf_filter_gene.jpg",width = 8,height = 8,dpi = 350)


length(unique(tmp.res.LR$Lgene))



