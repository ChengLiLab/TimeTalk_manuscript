####TimeTalk manuscript
####author: wlt
####20220308
####20220411
####This version of script will 
####under go minor revsion in the near future

####Fig1: Show data integration



####-------0. load package---------------

library(RColorBrewer)
library(ggplot2)
library(cowplot)
library(tidyverse)
library(clusterProfiler)
library(future.apply)
library(circlize)
library(ggthemes)
library(ggsci)
library(ggsignif)
library(preprocessCore)
library(factoextra)
library(Seurat)
library(monocle)
library(rtracklayer)
library(GenomicFeatures)
library(patchwork)
library(EnrichedHeatmap)
library(ComplexHeatmap)
#library(CellChat)
library(rstatix)
source(file = "code/myUtils.R")


# ####-------1. Preprocess data-------------------
### The preprocess procedure don't need to run again
### If the readers want to to do reanalysis, just cancel the comments of the 
### following code

### Will include more cells from other datasets in the near future

# ###-----1.1 prepare Oocyte data------
# tmp.1 <- read_delim(file = "data/GSM967487_MII_oocyte_CASTEi_1_RPKM.txt.gz",delim = "\t",skip = 1) %>%
#   dplyr::rename(Oocyte1 = `RPKM/FPKM`,Gene_symbol=`Gene symbol`) %>%
#   dplyr::select(Gene_symbol,Oocyte1)
# 
# tmp.2 <- read_delim(file = "data/GSM967488_MII_oocyte_CASTEi_2_RPKM.txt.gz",delim = "\t",skip = 1) %>%
#   dplyr::rename(Oocyte2 = `RPKM/FPKM`,Gene_symbol=`Gene symbol`) %>%
#   dplyr::select(Gene_symbol,Oocyte2)
# 
# tmp.3 <- read_delim(file = "data/GSM967489_MII_oocyte_CASTEi_3_RPKM.txt.gz",delim = "\t",skip = 1) %>%
#   dplyr::rename(Oocyte3 = `RPKM/FPKM`,Gene_symbol=`Gene symbol`) %>%
#   dplyr::select(Gene_symbol,Oocyte3)
# 
# 
# ####check the symbol coordination
# sum(tmp.1$Gene_symbol==tmp.2$Gene_symbol)
# sum(tmp.2$Gene_symbol==tmp.3$Gene_symbol)
# 
# Oocyte.df <- cbind(tmp.1,Oocyte2=tmp.2$Oocyte2,Oocyte3=tmp.3$Oocyte3) 
# Oocyte.df <- Oocyte.df[!duplicated(Oocyte.df$Gene_symbol),]
# 
# #### fix gene name, eg. Bmp4|Gm15217	
# Oocyte.df$Gene_symbol <- gsub(pattern = "\\|.*",replacement = "",x = Oocyte.df$Gene_symbol)
# 
# #### save results
# saveRDS(Oocyte.df,file = myFileName(prefix = "res/R/Oocyte_GSE38495_raw_rpkm_",suffix = ".rds"))
# 
# 
# ###-----1.2 prepare Deng et.al 2013---------
# 
# ###-----1.2.1 meta data-----
# dir <- "data/test/GSE45719_RAW/"
# samplefiles <- list.files(dir)

# tmp.idx <- which(unlist(lapply(strsplit(samplefiles,split = "_"),length))==3)
# samplefiles[tmp.idx]
# 
# tmp.idx <- which(unlist(lapply(strsplit(samplefiles,split = "_"),length))==5)
# samplefiles[tmp.idx]

# sampleinfor <- gsub(pattern = "_expression.txt","",x = samplefiles)
# sampleinfor <- gsub(pattern = "GSM[0-9]{7}_",replacement = "",sampleinfor)
# GSM <- regmatches(x = samplefiles,m = regexpr(pattern = "GSM[0-9]{7}",text = samplefiles))
# GSM
# 
# sample.df <- data.frame(samplefiles=samplefiles,
#                         sampleinfor=sampleinfor,
#                         GSM=GSM,
#                         embryo_id=gsub(pattern = "-[0-9]*",replacement = "",sampleinfor),
#                         stringsAsFactors = F)
# sample.df$Stage <- gsub(pattern = "_[0-9].*",replacement = "",sample.df$embryo_id)
# length(unique(sample.df$sampleinfor))
# 
# sample.df$embryo_id <- gsub(pattern = "fibroblast_.*_",replacement = "",sample.df$embryo_id)
# sample.df$embryo_id[grep("BxC|CxB",sample.df$embryo_id,perl = T)] <- paste0("fibroblast","_",sample.df$embryo_id[grep("BxC|CxB",sample.df$embryo_id,perl = T)])
# 
# ### check stages
# table(sample.df$Stage)


####-----1.2.2 get gene expression-----
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


####-----1.2.3 merge to get beforeEPI data-------

# #### load before exp data
# data.exp <- readRDS("res/R/deng_science_2013_rpkm.Rds")
# tmp.select <- c("zy1","zy2","zy3","zy4",
#                 "early2cell_0r-1","early2cell_0r-2",
#                 "early2cell_1-1","early2cell_1-2", 
#                 "early2cell_2-1","early2cell_2-2",
#                 "early2cell_3-1","early2cell_3-2",
#                 "mid2cell_0r-1","mid2cell_0r-2",  
#                 "mid2cell_3-1","mid2cell_3-2",
#                 "mid2cell_4-1","mid2cell_4-2",   
#                 "mid2cell_5-1","mid2cell_5-2",
#                 "mid2cell_6-1","mid2cell_6-2",   
#                 "mid2cell_7-1","mid2cell_7-2",
#                 "late2cell_5-1","late2cell_5-2",
#                 "late2cell_6-1","late2cell_6-2",  
#                 "late2cell_7-1","late2cell_7-2",
#                 "late2cell_8-1","late2cell_8-2", 
#                 "late2cell_9-1","late2cell_9-2",
#                 "4cell_1-1","4cell_1-2",
#                 "4cell_1-4","4cell_2-1",
#                 "4cell_2-2","4cell_2-3",
#                 "4cell_2-4","4cell_3-1",
#                 "4cell_3-3","4cell_3-4",
#                 "4cell_4-1","4cell_4-2",
#                 "4cell_4-3","4cell_4-4",
#                 "8cell_1-1","8cell_1-2",            
#                 "8cell_1-4","8cell_1-5",            
#                 "8cell_1-6","8cell_1-7",
#                 "8cell_1-8","8cell_2-1",            
#                 "8cell_2-2","8cell_2-3",            
#                 "8cell_2-4","8cell_2-6",            
#                 "8cell_2-7","8cell_2-8",            
#                 "8cell_5-1","8cell_5-2",            
#                 "8cell_5-3","8cell_5-4",            
#                 "8cell_5-6","8cell_5-7",            
#                 "8cell_5-8","8cell_8-1",            
#                 "8cell_8-2","8cell_8-3",            
#                 "8cell_8-4","8cell_8-6",            
#                 "8cell_8-7","8cell_8-8")
# 
# tmp.select <- c(tmp.select,grep("split",grep("16cell",sampleinfor,perl = T,value = T),value = T,invert = T))
# tmp.select <- c(tmp.select,grep("earlyblast",sampleinfor,perl=T,value = T))
# tmp.select <- c(tmp.select,grep("midblast",sampleinfor,perl=T,value = T))
# 
# tmp.select <- c(tmp.select,grep("lateblast",sampleinfor,perl=T,value = T))
# tmp.select <- c(tmp.select,grep("fibroblast",sampleinfor,perl = T,value = T))
# 
# tmp.sample.df <- data.frame(cell_id=tmp.select,
#                             cell_type=tmp.select,
#                             embryo_id=gsub(pattern = "-.*",
#                                            replacement = "",
#                                            gsub(pattern = ".*_",
#                                                 replacement = "",
#                                                 tmp.select)),
#                             embryo_cell_id = gsub(pattern = ".*-",
#                                                   replacement = "",
#                                                   tmp.select),
#                             Stage=gsub(pattern = "_.*",
#                                        replacement = "",
#                                        tmp.select),
#                             stringsAsFactors = F)
# tmp.sample.df$Stage[1:4] <- "zygote"
# 
# 
# data.exp <- data.exp[,tmp.select]
# beforeEPI.exp <- data.exp %>%
#   rownames_to_column("Gene_symbol")
# 
# 
# ### fix gene id
# ### 20210114
# beforeEPI.exp$Gene_symbol <- gsub(pattern = "\\|.*",replacement = "",beforeEPI.exp$Gene_symbol)
# 
# 
# ####load Oocyte
# Oocyte.exp <- readRDS(file = "res/R/Oocyte_GSE38495_raw_rpkm__20210114.rds")
# Oocyte.meta <- data.frame(cell_id=colnames(Oocyte.exp)[-1],
#                           cell_type=colnames(Oocyte.exp)[-1],
#                           embryo_id=colnames(Oocyte.exp)[-1],
#                           embryo_cell_id=colnames(Oocyte.exp)[-1],
#                           Stage=rep("Oocyte",3),stringsAsFactors = F)
# 
# 
# ####merge data
# beforeEPI.exp <- inner_join(Oocyte.exp,beforeEPI.exp,by="Gene_symbol")
# beforeEPI.exp <- beforeEPI.exp[!duplicated(beforeEPI.exp$Gene_symbol),]
# 
# rownames(beforeEPI.exp) <- NULL
# beforeEPI.exp <- beforeEPI.exp %>%
#   column_to_rownames("Gene_symbol")
# 
# 
# beforeEPI.meta <- rbind(Oocyte.meta,tmp.sample.df)
# 
# 
# 
# ####save to rds
# 
# saveRDS(beforeEPI.exp,file = "res/R/beforeEPI_exp_20210114.rds")
# saveRDS(beforeEPI.meta,file = "res/R/beforeEPI_meta_20210114.rds")

# ####-----------1.3 Add xiewei data-------------------------------
# xiewei.data.exp <- read.delim(file = "data/test/GSE66582_stage_FPKM.txt",stringsAsFactors = F)
# xiewei.data.exp <- xiewei.data.exp %>%
#   column_to_rownames("gene")
# boxplot(xiewei.data.exp)
# colnames(xiewei.data.exp)[4:6] <- c("2cell","4cell","8cell")
# colnames(xiewei.data.exp)
# saveRDS(object = xiewei.data.exp,file = "res/R/xiewei.data.exp_20201012.rds")




####----------2.QC------------------------
####----------2.1 load data-----------

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

### quantile normalize
tmp.names.1 <- colnames(beforeEPI.exp.mean)
tmp.mat.1 <- preprocessCore::normalize.quantiles(as.matrix(beforeEPI.exp.mean))
colnames(tmp.mat.1) <- tmp.names.1

tmp.sample <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                "4cell","8cell","16cell","earlyblast","midblast",  
                "lateblast")
tmp.mat.1 <- tmp.mat.1[,tmp.sample]
####get boxplot
boxplot(tmp.mat.1)



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

beforeEPI.exp.mean$ICM <- beforeEPI.exp.mean %>%
  dplyr::select(c("earlyblast","midblast","lateblast")) %>%
  rowMeans()


####-------revised version-----------
Stage <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell","ICM")
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
         width = 16,height = 16,dpi=350)

myggsave(pp,prefix = "res/fig/figS2_QC_beforeEPI_cor_revised",
         suffix = ".pdf",
         width = 16,height = 16,dpi=350)


#####-----2.2.2 Fig 1B Correlation pheatmap-----

Stage <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell","ICM")
data.plot.1 <- beforeEPI.exp.mean[gene_symbols,]
data.plot.2 <- xiewei.data.exp[gene_symbols,]


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

tmp.stage  <- c("MII oocyte","zygote","early 2-cell","late 2-cell","4-cell","8-cell","ICM")

png(filename = myFileName(prefix = "res/fig/Fig1_QC_pheatmap",suffix = ".png"),
     width = 8,height = 8,res = 350,units = "in")
exprTable_t <- as.data.frame(t(res.mat))
col_dist = dist(exprTable_t)
hclust_1 <- hclust(col_dist)
manual_order = rownames(res.mat)
#dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))))
dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))), agglo.FUN = max)
col_cluster <- as.hclust(dend)
ComplexHeatmap::pheatmap(res.mat,fontsize = 20,
         show_rownames = T,
         display_numbers = T,
         number_color = "white",
         cluster_rows = col_cluster,
         cluster_cols = col_cluster,
         fontsize_number = 22,
         name = "PCC",border_color = "black",
         labels_row = paste0("bulk ",tmp.stage),
         labels_col = paste0("pseudo-bulk ",tmp.stage),
         show_colnames = T,
         heatmap_legend_param = list(title_gp = gpar(fontsize = 20,fontface = "bold"),
                                     labels_gp=gpar(fontsize = 20)),
         angle_col = "315",
         color = rdbu(100))
dev.off()



pdf(file  = myFileName(prefix = "res/fig/Fig1_QC_pheatmap",suffix = ".pdf"),
    width = 8,height = 8)
exprTable_t <- as.data.frame(t(res.mat))
col_dist = dist(exprTable_t)
hclust_1 <- hclust(col_dist)
manual_order = rownames(res.mat)
#dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))))
dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))), agglo.FUN = max)
col_cluster <- as.hclust(dend)

ComplexHeatmap::pheatmap(res.mat,fontsize = 20,
                         show_rownames = T,
                         display_numbers = T,
                         number_color = "white",
                         cluster_rows = col_cluster,
                         cluster_cols = col_cluster,
                         fontsize_number = 22,
                         name = "PCC",border_color="black",
                         labels_row = paste0("bulk ",tmp.stage),
                         labels_col = paste0("pseudo-bulk ",tmp.stage),
                         show_colnames = T,
                         heatmap_legend_param = list(title_gp = gpar(fontsize = 20,fontface = "bold"),
                                                     labels_gp=gpar(fontsize = 20)),
                         angle_col = "315",
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
  dplyr::filter(Stage!="fibroblast") %>%
  pull(cell_id)

data.plot <- data.plot[,tmp.select]
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)

test.res <- pca.res$rotation[,1:2]




#fviz_eig

eig <- factoextra::get_eigenvalue(pca.res) %>%
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

######Next will put into the main figure
######-------------2.3.2  add xiewei's data-----------------

gene.use <- intersect(rownames(xiewei.data.exp),rownames(beforeEPI.exp.mean))
data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  dplyr::filter(Stage!="fibroblast") %>%
  pull(cell_id)
data.plot <- data.plot[gene.use,tmp.select]

tmp.select <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell","ICM")
data.plot <- cbind(xiewei.data.exp[gene.use,tmp.select],data.plot)
data.plot <- median_center(data.plot)
pca.res <- prcomp(t(data.plot),center = T,scale. = F)
eig <- get_eigenvalue(pca.res) %>%
  pull("variance.percent")
ext_labels <- paste0(round(eig, 1), "%")

####add meta data for xieweidata
tmp.select <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell","ICM")
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

data.plot$Stage[c(1,8:10)] <- "MII_oocyte" 

Stage <- c("MII_oocyte","zygote","early2cell","mid2cell","late2cell",
           "4cell","8cell","16cell","ICM","earlyblast","midblast",  
           "lateblast")

Stage.1 <- c("MII oocyte","zygote","early 2-cell",
             "mid 2-cell","late 2-cell","4-cell",
             "8-cell","16-cell","ICM",
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
                               "#4E1550","#8EAF3D","grey",
                               "#F05153","#5358A5","#79A7AA"))+
  #scale_shape_manual(values = c(25,17,16,19,18,15,9))+
  theme_cowplot(font_size = 24)+
  theme(legend.position = "right",aspect.ratio = 1,
        axis.line.x = element_line( size = 2 ),
        axis.line.y = element_line( size = 2 ),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))
p
myggsave(p,prefix = "res/fig/Fig1_QC_pca",suffix = ".jpg",width = 8,height = 8,dpi = 350)
myggsave(p,prefix = "res/fig/Fig1_QC_pca",suffix = ".pdf",width = 8,height = 8,dpi = 350)




#####-------2.4 boxplot---------------

#####-------2.4.1 beforeEPI--------------

#####------2.4.1.1 raw unquantile normlized------

tmp.stage <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                "4cell","8cell","16cell","earlyblast","midblast",  
                "lateblast")
data.plot.1 <- beforeEPI.exp.mean[,tmp.stage]
data.plot.1 <- data.plot.1 %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = "gene.exp",-Gene) 
data.plot.1$Stage <- factor(data.plot.1$Stage,tmp.stage)
tmp.stage <-  c("MII oocyte","zygote","early 2-cell","mid 2-cell","late 2-cell",
                "4-cell","8-cell","16-cell","early blastocyst","mid blastocyst",  
                "late blastocyst")
p <- myBoxplot_advanced(data = data.plot.1,
                        x = "Stage",
                        y = "gene.exp",
                        fontsize = 30,
                        axislinesize = 2,
                        fill = "Stage",color = NULL,
                        mycolor = col.spectral(length(tmp.stage)),
                        ylimts = c(0,20),ybreaks = seq(0,20,5))+
  scale_x_discrete(labels=tmp.stage)+
  ylab("gene.exp.mean")+
  xlab(NULL)+
  theme(legend.position = "none",axis.ticks = element_line(size = 2))
p
myggsave(p = p,
         prefix = "res/fig/figS1_QC_pseudobulk_boxplot_raw",
         suffix = ".png",
         width = 12,height = 9,dpi = 350)

myggsave(p = p,
         prefix = "res/fig/figS1_QC_pseudobulk_boxplot_raw",
         suffix = ".pdf",useDingbats=F,
         width = 12,height = 9,dpi = 350)


####------------2.4.1.2 normalized boxplot------------------

tmp.stage <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
               "4cell","8cell","16cell","earlyblast","midblast",  
               "lateblast")
data.plot.1 <- beforeEPI.exp.mean[,tmp.stage]

tmp.names.1 <- colnames(data.plot.1)
tmp.names.2 <- rownames(data.plot.1)
data.plot.1 <- preprocessCore::normalize.quantiles(as.matrix(data.plot.1))
colnames(data.plot.1) <- tmp.names.1
rownames(data.plot.1) <- tmp.names.2
data.plot.1 <- data.plot.1 %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = "gene.exp",-Gene) 



data.plot.1$Stage <- factor(data.plot.1$Stage,tmp.stage)
tmp.stage <-  c("MII oocyte","zygote","early 2-cell","mid 2-cell","late 2-cell",
                "4-cell","8-cell","16-cell","early blastocyst","mid blastocyst",  
                "late blastocyst","E5.25","E5.5","E6.25","E6.5")
p <- myBoxplot_advanced(data = data.plot.1,
                        x = "Stage",
                        y = "gene.exp",
                        fontsize = 30,
                        axislinesize = 2,
                        fill = "Stage",color = NULL,
                        mycolor = col.spectral(length(tmp.stage)),
                        ylimts = c(0,20),ybreaks = seq(0,20,5))+
  scale_x_discrete(labels=tmp.stage)+
  ylab("gene.exp.mean")+
  xlab(NULL)+
  theme(legend.position = "none",axis.ticks = element_line(size = 2))

myggsave(p,
         prefix = "res/fig/figS1_QC_pseudobulk_boxplot_quantile_norm",
         suffix = ".png",
         width = 12,height = 9,dpi = 350)

myggsave(p,
         prefix = "res/fig/figS1_QC_pseudobulk_boxplot_quantile_norm",
         suffix = ".pdf",
         width = 12,height = 9,dpi = 350)




####-------2.4.2 Prof. xiewei data -----

#####------2.4.2.1 unormlized---------------
tmp.stage <- c("MII_oocyte","zygote","early2cell","late2cell",
               "4cell","8cell","ICM")
data.plot.2 <- xiewei.data.exp[,tmp.stage]
# tmp.names.1 <- colnames(data.plot.2)
# tmp.names.2 <- rownames(data.plot.2)
# data.plot.2 <- preprocessCore::normalize.quantiles(as.matrix(data.plot.2))
# colnames(data.plot.2) <- tmp.names.1
# rownames(data.plot.2) <- tmp.names.2


data.plot.2 <- data.plot.2 %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = gene.exp,-Gene) %>%
  mutate(group = "bulk")


data.plot.2$Stage <- factor(data.plot.2$Stage,levels = tmp.stage)

tmp.stage <- c("MII oocyte","zygote","early 2-cell","late 2-cell",
               "4-cell","8-cell","ICM")
p <- myBoxplot_advanced(data = data.plot.2,
                        x = "Stage",
                        y = "gene.exp",
                        fill = "Stage",
                        color = NULL,fontsize = 30,axislinesize = 2,
                        mycolor = col.spectral(length(tmp.stage)),
                        ylimts = c(0,20),ybreaks = seq(0,20,5))+
  scale_x_discrete(labels=tmp.stage)+
  ylab("log2(RPKM+1)")+
  xlab(NULL)+
  theme(legend.position = "none",axis.ticks = element_line(size = 2))
p
myggsave(p = p,prefix = "res/fig/figS1_QC_bulk_boxplot",
         suffix = ".jpg",width = 8,height = 8,dpi = 350)
myggsave(p = p,prefix = "res/fig/figS1_QC_bulk_boxplot",
         suffix = ".pdf",width = 8,height = 8,dpi = 350)





####normalized data

tmp.stage <- c("MII_oocyte","zygote","early2cell","late2cell",
               "4cell","8cell","ICM")
data.plot.2 <- xiewei.data.exp[,tmp.stage]
tmp.names.1 <- colnames(data.plot.2)
tmp.names.2 <- rownames(data.plot.2)
data.plot.2 <- preprocessCore::normalize.quantiles(as.matrix(data.plot.2))
colnames(data.plot.2) <- tmp.names.1
rownames(data.plot.2) <- tmp.names.2


data.plot.2 <- data.plot.2 %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  gather(key = "Stage",value = gene.exp,-Gene) %>%
  mutate(group = "bulk")


data.plot.2$Stage <- factor(data.plot.2$Stage,levels = tmp.stage)

tmp.stage <- c("MII oocyte","zygote","early 2-cell","late 2-cell",
               "4-cell","8-cell","ICM")
p <- myBoxplot_advanced(data = data.plot.2,
                        x = "Stage",
                        y = "gene.exp",
                        fill = "Stage",
                        color = NULL,fontsize = 30,axislinesize = 2,
                        mycolor = col.spectral(length(tmp.stage)),
                        ylimts = c(0,20),ybreaks = seq(0,20,5))+
  scale_x_discrete(labels=tmp.stage)+
  ylab("log2(RPKM+1)")+
  xlab(NULL)+
  theme(legend.position = "none",axis.ticks = element_line(size = 2))
p
myggsave(p = p,prefix = "res/fig/figS1_QC_bulk_boxplot_normalized",
         suffix = ".jpg",width = 8,height = 8,dpi = 350)
myggsave(p = p,prefix = "res/fig/figS1_QC_bulk_boxplot_normalized",
         suffix = ".pdf",width = 8,height = 8,dpi = 350)




####-------2.5 merge to seurat------------------

####load data
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)

beforeEPI.meta <- beforeEPI.meta %>%
  dplyr::filter(Stage!="fibroblast") 
beforeEPI.exp <- beforeEPI.exp[,beforeEPI.meta$cell_id]


tmp.levels.1 <- c("Oocyte","zygote","early2cell","mid2cell","late2cell",
                  "4cell","8cell","16cell","earlyblast","midblast","lateblast")
tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst")

beforeEPI.meta$Stage <- plyr::mapvalues(beforeEPI.meta$Stage,
                                             from = tmp.levels.1,
                                             to = tmp.levels.2)
beforeEPI.meta$Stage <- factor(beforeEPI.meta$Stage,
                                    levels = tmp.levels.2)


####rename cell types
beforeEPI.meta$cell_type <- beforeEPI.meta$Stage
rownames(beforeEPI.meta) <- beforeEPI.meta$cell_id



#####-------2.5.2 Run seurat steps---------
seu <- myCreateSeurat(data = beforeEPI.exp,
                      meta.data = beforeEPI.meta)

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
seu <- RunUMAP(seu,
               dims = c(1:20),
               umap.method = "umap-learn")
seu <- FindNeighbors(seu,dims = c(1:20))
seu <- FindClusters(seu,resolution = 0.4)
colnames(seu[[]])

p1 <- TSNEPlot(seu,group.by="RNA_snn_res.0.4")
p2 <- TSNEPlot(seu,group.by="Stage")

p1 | p2

p1 | VlnPlot(seu,features = c("Cdx2","Gata3","Sox2","Krt8","Nanog"))
ggsave(filename = myFileName(prefix = "res/fig/Fig1_seurat_showmarker_assign_celltypes",suffix = ".png"),
       width = 12,height = 8,dpi = 350)


### Cdx2, TE
TE.cells <- WhichCells(seu,idents = c("1"))
### ICM. cell
ICM.cells <- WhichCells(seu,idents = c("2"))

TE.cells <- seu[[]][TE.cells,] %>%
  dplyr::filter(Stage %in% c("early blastocyst","mid blastocyst","late blastocyst")) %>%
  rownames()

ICM.cells <- seu[[]][ICM.cells,] %>%
  dplyr::filter(Stage %in% c("early blastocyst","mid blastocyst","late blastocyst")) %>%
  rownames()

### Reassign Cell types
Idents(seu) <- "cell_type"
seu <- SetIdent(seu,cells = TE.cells,value = "TE")
seu <- SetIdent(seu,cells = ICM.cells,value = "ICM")

seu$cell_type <- Idents(seu)



tmp.levels <- c("MII oocyte","zygote",
                "early 2-cell","mid 2-cell","late 2-cell",
                "4-cell","8-cell","16-cell","ICM","TE")
seu$cell_type <- factor(seu$cell_type,levels = tmp.levels)
Idents(seu) <- "cell_type"
#### check assign
# table(seu$cell_type)
# tmp.cells <- WhichCells(seu,expression = Stage %in% c("early blastocyst","mid blastocyst","late blastocyst"))
DimPlot(seu,reduction = "tsne")





saveRDS(seu,file = myFileName(prefix = "res/R/early.scRNAseq.seurat",suffix = ".rds"))



#####-------2.5.3 tsne plot ------------------------
#####UMAP can't see the lineage.....


#####-----2.5.3.1 group.by stage------
seu <- readRDS(file = "res/R/early.scRNAseq.seurat_2022031015.rds")
colnames(seu[[]])
table(seu$cell_type)
tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst")
p <- TSNEPlot(seu,label=F,group="Stage",pt.size=3)+
  ggtitle(NULL)+
  scale_color_manual(values = col.spectral(length(tmp.levels.2)))
  


myFunction <- function(p,tmp.ratio=0.3,tmp.fontsize=6){
  tmp.range.x <- round(range(p$data$tSNE_1))-1
  tmp.range.y <- round(range(p$data$tSNE_2))-1
  tmp.ticks.x <- seq(tmp.range.x[1],tmp.range.x[2],length.out=10)
  tmp.ticks.y <- seq(tmp.range.y[1],tmp.range.y[2],length.out=10)
  tmp.step.x <- tmp.ticks.x[2]-tmp.ticks.x[1]
  tmp.step.y <- tmp.ticks.y[2]-tmp.ticks.y[1]
  list(geom_path(data = data.frame(x=tmp.ticks.x[1:2],
                                   y=tmp.ticks.y[c(1,1)]),
                 size = 1,
                 mapping = aes(x,y),
                 linejoin = "bevel",
                 lineend = "round",
                 arrow = arrow(angle = 45,type = "closed",length = unit(0.1, "inches"))),
       geom_path(data = data.frame(x=tmp.ticks.x[c(1,1)],
                                   y=tmp.ticks.y[1:2]),
                 size = 1,
                 mapping = aes(x,y),
                 linejoin = "bevel",
                 lineend = "round",
                 arrow = arrow(angle = 45,type = "closed",length = unit(0.1, "inches"))),
       annotate("text",x = mean(tmp.ticks.x[1:2]),
                y = tmp.ticks.y[1]-tmp.ratio*tmp.step.y,label="tSNE1",fontface="bold",size=tmp.fontsize),
       annotate("text",x = tmp.ticks.x[1]-tmp.ratio*tmp.step.x,
                y = mean(tmp.ticks.y[1:2]),label="tSNE2",fontface="bold",size=tmp.fontsize,angle=90),
       NoAxes(),
       theme(legend.text = element_text(size = 24)))
}

p + myFunction(p,tmp.ratio = 0.8,tmp.fontsize = 8)  


ggsave(filename = myFileName(prefix = "res/fig/Fig1_tsne_by_stage",suffix = ".jpg"),
       width = 8,height = 6,dpi = 350)

ggsave(filename = myFileName(prefix = "res/fig/Fig1_tsne_by_stage",suffix = ".pdf"),
       width = 8,height = 6,useDingbats=F)

#####-----2.5.3.2 group.by cell_type show markers------

table(seu$cell_type)
tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "ICM","TE")
p <- TSNEPlot(seu,label=F,group="cell_type",pt.size=3)+
  ggtitle(NULL)+
  scale_color_manual(values = col.spectral(length(tmp.levels.2)))



myFunction <- function(p,tmp.ratio=0.3,tmp.fontsize=6){
  tmp.range.x <- round(range(p$data$tSNE_1))-1
  tmp.range.y <- round(range(p$data$tSNE_2))-1
  tmp.ticks.x <- seq(tmp.range.x[1],tmp.range.x[2],length.out=10)
  tmp.ticks.y <- seq(tmp.range.y[1],tmp.range.y[2],length.out=10)
  tmp.step.x <- tmp.ticks.x[2]-tmp.ticks.x[1]
  tmp.step.y <- tmp.ticks.y[2]-tmp.ticks.y[1]
  list(geom_path(data = data.frame(x=tmp.ticks.x[1:2],
                                   y=tmp.ticks.y[c(1,1)]),
                 size = 1,
                 mapping = aes(x,y),
                 linejoin = "bevel",
                 lineend = "round",
                 arrow = arrow(angle = 45,type = "closed",length = unit(0.1, "inches"))),
       geom_path(data = data.frame(x=tmp.ticks.x[c(1,1)],
                                   y=tmp.ticks.y[1:2]),
                 size = 1,
                 mapping = aes(x,y),
                 linejoin = "bevel",
                 lineend = "round",
                 arrow = arrow(angle = 45,type = "closed",length = unit(0.1, "inches"))),
       annotate("text",x = mean(tmp.ticks.x[1:2]),
                y = tmp.ticks.y[1]-tmp.ratio*tmp.step.y,label="tSNE1",fontface="bold",size=tmp.fontsize),
       annotate("text",x = tmp.ticks.x[1]-tmp.ratio*tmp.step.x,
                y = mean(tmp.ticks.y[1:2]),label="tSNE2",fontface="bold",size=tmp.fontsize,angle=90),
       NoAxes(),
       theme(legend.text = element_text(size = 24)))
}

p + myFunction(p,tmp.ratio = 0.8,tmp.fontsize = 8)  


ggsave(filename = myFileName(prefix = "res/fig/Fig1_tsne_by_cell_type",suffix = ".jpg"),
       width = 8,height = 6,dpi = 350)

ggsave(filename = myFileName(prefix = "res/fig/Fig1_tsne_by_cell_type",suffix = ".pdf"),
       width = 8,height = 6,useDingbats=F)





#####show markers



p1 <- FeaturePlot(seu,features = c("Nlrp5","Zscan4c","Gata3","Cdx2","Sox2","Sox7","Nanog"),
                  reduction = "tsne",pt.size = 3) +
  plot_layout(ncol = 3) &
  theme_cowplot(font_size = 22) &
  scale_color_gradientn(colors = rdwhbu(100)) &
  theme(plot.title = element_text(face="bold.italic",hjust = 0.5,size = 28),
        axis.line = element_line(size = 2),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(size = 2)) &
  NoAxes()
#?FeaturePlot
p1

ggsave(p1,filename = myFileName(prefix = "res/fig/FigS2_tsne_feature",suffix = ".png"),
       width = 8,height = 12,dpi = 350)
ggsave(p1,filename = myFileName(prefix = "res/fig/FigS2_tsne_feature",suffix = ".pdf"),
       width = 8,height = 12,dpi = 350)

tmp.markers <- c("Nlrp5","Zscan4f","Zscan4c","Gata3","Cdx2","Sox2","Nanog","Klf4","Gata6")

p1 <- FeaturePlot(seu,features = tmp.markers,
                  reduction = "tsne",pt.size = 3) +
  plot_layout(ncol = 2) &
  theme_cowplot(font_size = 22) &
  scale_color_gradientn(colors = rdwhbu(100)) &
  theme(plot.title = element_text(face="bold.italic",hjust = 0.5,size = 28),
        axis.line = element_line(size = 2),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(size = 2)) &
  NoAxes()
#?FeaturePlot
p1

ggsave(p1,filename = myFileName(prefix = "res/fig/FigS2_tsne_feature",suffix = ".png"),
       width = 6,height = 12,dpi = 350)
ggsave(p1,filename = myFileName(prefix = "res/fig/FigS2_tsne_feature",suffix = ".pdf"),
       width = 6,height = 12,dpi = 350)


p1 <- FeaturePlot(seu,features = tmp.markers,
                  reduction = "tsne",pt.size = 3) +
  plot_layout(ncol = 3) &
  theme_cowplot(font_size = 22) &
  scale_color_gradientn(colors = rdwhbu(100)) &
  theme(plot.title = element_text(face="bold.italic",hjust = 0.5,size = 28),
        axis.line = element_line(size = 2),
        axis.title = element_text(size = 18),
        axis.ticks = element_line(size = 2)) &
  NoAxes()
#?FeaturePlot
p1

ggsave(p1,filename = myFileName(prefix = "res/fig/FigS2_tsne_feature",suffix = ".png"),
       width = 12,height = 9,dpi = 350)

ggsave(p1,filename = myFileName(prefix = "res/fig/FigS2_tsne_feature",suffix = ".jpg"),
       width = 12,height = 9,dpi = 350)

ggsave(p1,filename = myFileName(prefix = "res/fig/FigS2_tsne_feature",suffix = ".pdf"),
       width = 12,height = 9,dpi = 350)

####----------3.Ordering Cell by monocle2------------------

####----------3.1 load Seurat -------------


seu <- readRDS(file = "res/R/early.scRNAseq.seurat_2022031015.rds")

####load data
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)

beforeEPI.meta <- beforeEPI.meta %>%
  dplyr::filter(Stage!="fibroblast") 
beforeEPI.exp <- beforeEPI.exp[,beforeEPI.meta$cell_id]


####Define function

tmp.convert <- function (x, assay = NULL, reduction = NULL,expressionFamily=NULL) {
  
  if (!Seurat:::PackageCheck("monocle", error = FALSE)) {
    stop("Please install monocle from Bioconductor before converting to a CellDataSet object")
  }
  else if (packageVersion(pkg = "monocle") >= package_version(x = "2.99.0")) {
    stop("Seurat can only convert to/from Monocle v2.X objects")
  }
  assay <- assay %||% DefaultAssay(object = x)
  counts <- GetAssayData(object = x, assay = assay, slot = "data")
  cell.metadata <- x[[]]
  feature.metadata <- x[[assay]][[]]
  if (!"gene_short_name" %in% colnames(x = feature.metadata)) {
    feature.metadata$gene_short_name <- rownames(x = feature.metadata)
  }
  pd <- new(Class = "AnnotatedDataFrame", data = cell.metadata)
  fd <- new(Class = "AnnotatedDataFrame", data = feature.metadata)
  cds <- monocle::newCellDataSet(cellData = counts, phenoData = pd, 
                                 featureData = fd, expressionFamily = expressionFamily)
  if ("monocle" %in% names(x = Misc(object = x))) {
    monocle::cellPairwiseDistances(cds = cds) <- Misc(object = x, 
                                                      slot = "monocle")[["cellPairwiseDistances"]]
    monocle::minSpanningTree(cds = cds) <- Misc(object = x, 
                                                slot = "monocle")[["minSpanningTree"]]
    Biobase::experimentData(cds = cds) <- Misc(object = x, 
                                               slot = "monocle")[["experimentData"]]
    Biobase::protocolData(cds = cds) <- Misc(object = x, 
                                             slot = "monocle")[["protocolData"]]
    Biobase::classVersion(cds = cds) <- Misc(object = x, 
                                             slot = "monocle")[["classVersion"]]
    slot(object = cds, name = "lowerDetectionLimit") <- Misc(object = x, 
                                                             slot = "monocle")[["lowerDetectionLimit"]]
    slot(object = cds, name = "dispFitInfo") <- Misc(object = x, 
                                                     slot = "monocle")[["dispFitInfo"]]
    slot(object = cds, name = "auxOrderingData") <- Misc(object = x, 
                                                         slot = "monocle")[["auxOrderingData"]]
    slot(object = cds, name = "auxClusteringData") <- Misc(object = x, 
                                                           slot = "monocle")[["auxClusteringData"]]
  }
  dr.slots <- c("reducedDimS", "reducedDimK", "reducedDimW", 
                "reducedDimA")
  reduction <- reduction %||% Seurat:::DefaultDimReduc(object = x, assay = assay)
  if (!is.null(x = reduction)) {
    if (grepl(pattern = "tsne", x = tolower(x = reduction))) {
      slot(object = cds, name = "dim_reduce_type") <- "tSNE"
      monocle::reducedDimA(cds = cds) <- t(x = Embeddings(object = x[[reduction]]))
    }
    else {
      slot(object = cds, name = "dim_reduce_type") <- reduction
      monocle::reducedDimA(cds = cds) <- Loadings(object = x[[reduction]])
      slot(object = cds, name = "reducedDimS") <- Embeddings(object = x[[reduction]])
    }
    for (ii in dr.slots) {
      if (ii %in% names(x = slot(object = x[[reduction]], 
                                 name = "misc"))) {
        slot(object = cds, name = ii) <- slot(object = x[[reduction]], 
                                              name = "misc")[[ii]]
      }
    }
  }
  return(cds)
}
early.cds <- tmp.convert(seu,reduction = "tsne",
                         expressionFamily = uninormal())

#####------3.2 run monlcle process------------------

tmp.cds <- early.cds
tmp.cds <- estimateSizeFactors(tmp.cds)
expressed_genes <- VariableFeatures(seu)
tmp.cds <- clusterCells(tmp.cds)
tmp.diff <- differentialGeneTest(tmp.cds[expressed_genes,],
                                 fullModelFormulaStr = "~Stage",
                                 cores = 10)
tmp.ordergene <- tmp.diff %>%
  dplyr::filter(qval < 0.01) %>%
  rownames()

tmp.cds <- setOrderingFilter(tmp.cds,tmp.ordergene)
tmp.cds <- reduceDimension(cds = tmp.cds,method="DDRTree",norm_method = "none")
tmp.cds <- orderCells(tmp.cds)

saveRDS(tmp.cds,file = myFileName(prefix = "res/R/early.scRNAseq.moncole",suffix = ".rds"))

#### The following test code will remove in the final version
# tmp.cells <- pData(tmp.cds) %>%
#   dplyr::filter(cell_type %in% c("TE","ICM")) %>%
#   rownames()
# 
# tmp.diff <- FindMarkers(seu,ident.1 = "TE",ident.2 = "ICM") 
# 
# 
# tmp.diff <- FindAllMarkers(seu)
# 
# tmp.diff.1 <- tmp.diff %>%
#   dplyr::filter(!(cluster %in% c("early 2-cell","mid 2-cell","late 2-cell"))) %>%
#   dplyr::filter(avg_logFC > 1.5 & p_val_adj < 0.01 )

# tmp.diff <- FindMarkers(seu,ident.1 = "MII oocyte",ident.2 = "zygote") 
# tmp.diff.2 <- tmp.diff %>%
#   dplyr::filter(abs(avg_logFC) > 1.5 & p_val < 0.01 ) %>%
#   rownames()
# 
# tmp.diff <- FindMarkers(seu,ident.1 = "TE",ident.2 = NULL) 
# tmp.diff.3 <- tmp.diff %>%
#   dplyr::filter(abs(avg_logFC) > 1.5 & p_val_adj < 0.01 ) %>%
#   rownames()
# 
# tmp.diff <- FindMarkers(seu,ident.1 = "ICM",ident.2 = NULL) 
# tmp.diff.4 <- tmp.diff %>%
#   dplyr::filter(abs(avg_logFC) > 1.5 & p_val_adj < 0.01 ) %>%
#   rownames()

# tmp.ordergene <- unique(c(tmp.diff.1,tmp.diff.2,tmp.diff.3,tmp.diff.4))
# 
# tmp.ordergene <- rownames(tmp.diff.1)

# tmp.cells <- pData(tmp.cds) %>%
#   dplyr::filter(!(cell_type %in% c("TE","ICM"))) %>%
#   rownames()
# 
# tmp.diff <- differentialGeneTest(tmp.cds[expressed_genes,tmp.cells],
#                                  fullModelFormulaStr = "~Stage",
#                                  cores = 10)



# tmp.ordergene <- tmp.diff %>%
#   dplyr::filter(qval < 0.01) %>%
#   rownames()




#tmp.ordergene <- VariableFeatures(seu)


#tmp.ordergene <- expressed_genes


#####------3.3 visulization-----------

##### please rember, we just need pseduotime

#####-----3.3.1 show pseudotime----------
tmp.cds <- readRDS(file = "res/R/early.scRNAseq.moncole_2022031116.rds")

str(tmp.cds@featureData@varMetadata)


tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst")

p1 <- plot_cell_trajectory(tmp.cds,color_by = "Stage",cell_size = 3)  +
  scale_color_manual(values = col.spectral(11))+
  theme_cowplot(font_size = 28)+
  NoAxes()
p2 <- plot_cell_trajectory(tmp.cds,color_by = "Pseudotime",cell_size = 3)+
  theme_cowplot(font_size = 28)+
  scale_color_gradientn(colors = col.spectral(100))+
  NoAxes()


p1 / p2 
ggsave(filename = myFileName(prefix = "res/fig/Fig1_monocle_DDRTree",
                             suffix = ".png"),
       width = 9,height = 9,dpi = 350)
ggsave(filename = myFileName(prefix = "res/fig/Fig1_monocle_DDRTree",
                             suffix = ".pdf"),
       width = 9,height = 9,dpi = 350)




######----3.3.2 adjust monocle2 tree-----------

# backbone 
reduced_dim_coords <- reducedDimK(tmp.cds)
ica_space_df <- data.frame(Matrix::t(reduced_dim_coords[c(1,2), ]))
colnames(ica_space_df) <- c("prin_graph_dim_1", "prin_graph_dim_2")
ica_space_df$sample_name <- row.names(ica_space_df)
ica_space_df$sample_state <- row.names(ica_space_df)
dp_mst <- minSpanningTree(tmp.cds)
if (is.null(dp_mst)) {
  stop("You must first call orderCells() before using this function")
}
library(igraph)
edge_list <- as.data.frame(get.edgelist(dp_mst))
colnames(edge_list) <- c("source", "target")
edge_df <- merge(ica_space_df, edge_list, 
                 by.x = "sample_name", 
                 by.y = "source", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "source_prin_graph_dim_1", 
                                   prin_graph_dim_2 = "source_prin_graph_dim_2"))
edge_df <- merge(edge_df, ica_space_df[, c("sample_name", 
                                           "prin_graph_dim_1", "prin_graph_dim_2")], by.x = "target", 
                 by.y = "sample_name", all = TRUE)
edge_df <- plyr::rename(edge_df, c(prin_graph_dim_1 = "target_prin_graph_dim_1", 
                                   prin_graph_dim_2 = "target_prin_graph_dim_2"))
S_matrix <- reducedDimS(tmp.cds)
data_df <- data.frame(t(S_matrix[c(1, 2), ]))
sample_state <- pData(tmp.cds)$State
data_df <- cbind(data_df, sample_state)
colnames(data_df) <- c("data_dim_1", "data_dim_2")
data_df$sample_name <- row.names(data_df)
lib_info_with_pseudo <- pData(tmp.cds)
data_df <- merge(data_df, lib_info_with_pseudo, 
                 by.x = "sample_name", 
                 by.y = "row.names")
return_rotation_mat <- function(theta) {
  theta <- theta/180 * pi
  matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 
         nrow = 2)
}
theta <- 180
tmp <- return_rotation_mat(theta) %*% t(as.matrix(data_df[, 
                                                          c(2, 3)]))
data_df$data_dim_1 <- tmp[1, ]
data_df$data_dim_2 <- tmp[2, ]
tmp <- return_rotation_mat(theta = theta) %*% t(as.matrix(edge_df[, 
                                                                  c("source_prin_graph_dim_1", "source_prin_graph_dim_2")]))
edge_df$source_prin_graph_dim_1 <- tmp[1, ]
edge_df$source_prin_graph_dim_2 <- tmp[2, ]
tmp <- return_rotation_mat(theta) %*% t(as.matrix(edge_df[, 
                                                          c("target_prin_graph_dim_1", "target_prin_graph_dim_2")]))
edge_df$target_prin_graph_dim_1 <- tmp[1, ]
edge_df$target_prin_graph_dim_2 <- tmp[2, ]



# ggplot()+
#   geom_point(data = data_df,aes())
# 
# 
# bg.plot <- data.frame(x = homeostasis.meta$x, y = homeostasis.meta$y, celltype = homeostasis.meta$New.Celltype)


# ggplot(data=homeostasis.meta, aes(x = x, y = y)) + 
#   geom_point(data=bg.plot, aes(x = x, y = y,fill = celltype), shape = 21, size=1, colour="grey", alpha=I(0.3)) +
#   geom_segment(aes_string(x = "source_prin_graph_dim_1", 
#                           y = "source_prin_graph_dim_2", 
#                           xend = "target_prin_graph_dim_1", 
#                           yend = "target_prin_graph_dim_2"), 
#                size = 0.75, 
#                linetype = "solid", na.rm = TRUE, data = edge_df) +
#   geom_point(aes(fill=factor(New.Celltype)), shape=21, colour="black", size=3) +
#   scale_fill_manual(name="Cell type", values = col.values) + 
#   scale_alpha(guide=F)+ theme(text = element_text(size=20)) + 
#   facet_wrap(~New.Celltype, nrow=4) + 
#   theme(legend.position="right", legend.key.height = grid::unit(0.35, "in")) + 
#   theme(legend.key = element_blank()) + 
#   theme(panel.background = element_rect(fill = "white"), axis.line = element_line(colour = "black")) +
#   xlab("Dimension 1") + 
#   ylab("Dimension 2")

tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst")

p <- ggplot()+
  geom_point(aes_string(x = "data_dim_1",
                        y = "data_dim_2",
                        color="Pseudotime"),
             size=3,
             data = data_df)+
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               size = 0.75, 
               linetype = "solid", 
               na.rm = TRUE, 
               data = edge_df)+
  #scale_color_manual(values = col.spectral(length(tmp.levels.2)))+
  scale_color_gradientn(colours = col.spectral(100))+
  theme_cowplot(font_size = 28)+
  theme(legend.position = "none")+
  NoAxes()
p
myggsave(p = p,
         prefix =  "res/fig/fig1_moncole2_for_model_pseudotime",
         suffix = ".pdf",width = 8,height = 8,dpi=350)
myggsave(p = p,
         prefix =  "res/fig/fig1_moncole2_for_model_pseudotime",
         suffix = ".png",width = 8,height = 8,dpi=350)


p1 <- ggplot()+
  geom_point(aes_string(x = "data_dim_1",
                        y = "data_dim_2",
                        color="Pseudotime"),
             size=3,
             data = data_df)+
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               size = 0.75, 
               linetype = "solid", 
               na.rm = TRUE, 
               data = edge_df)+
  #scale_color_manual(values = col.spectral(length(tmp.levels.2)))+
  scale_color_gradientn(colours = warm(100))+
  theme_cowplot(font_size = 28)+
  NoAxes()



p2 <- ggplot()+
  geom_point(aes_string(x = "data_dim_1",
                        y = "data_dim_2",
                        color="Stage"),
             size=3,
             data = data_df)+
  geom_segment(aes_string(x = "source_prin_graph_dim_1", 
                          y = "source_prin_graph_dim_2", 
                          xend = "target_prin_graph_dim_1", 
                          yend = "target_prin_graph_dim_2"), 
               size = 0.75, 
               linetype = "solid", 
               na.rm = TRUE, 
               data = edge_df)+
  scale_color_manual(values = col.spectral(length(tmp.levels.2)))+
  theme_cowplot(font_size = 28)+
  NoAxes()


p2 / p1
ggsave(filename = myFileName(prefix = "res/fig/Fig1_stage_pseudotime",suffix = ".png"),
       width = 8,height = 8,dpi = 350)
ggsave(filename = myFileName(prefix = "res/fig/Fig1_stage_pseudotime",suffix = ".pdf"),
       width = 8,height = 8,dpi = 350,useDingbats = T)



str(tmp.cds)

#####-----3.3.3.1 show marker--------------

tmp.markers <- c("Nlrp5","Zscan4c","Gata3","Cdx2","Sox2","Sox7","Nanog")


pData(tmp.cds)$Nlrp5 <- exprs(tmp.cds)["Nlrp5",]
pData(tmp.cds)$Zscan4c <- exprs(tmp.cds)["Zscan4c",]
pData(tmp.cds)$Gata3 <- exprs(tmp.cds)["Zscan4c",]
pData(tmp.cds)$Cdx2 <- exprs(tmp.cds)["Cdx2",]
pData(tmp.cds)$Sox2 <- exprs(tmp.cds)["Sox2",]
pData(tmp.cds)$Sox7 <- exprs(tmp.cds)["Sox7",]
pData(tmp.cds)$Nanog <- exprs(tmp.cds)["Nanog",]


plot_cell_trajectory(tmp.cds,color_by = "Sox2",cell_size = 3)+
  scale_color_gradientn(colors = rdwhbu(100))+
  NoAxes()


plot.list <- lapply(tmp.markers,FUN = function(ii){
  p <- plot_cell_trajectory(tmp.cds,color_by = ii,cell_size = 3)+
    scale_color_gradientn(colors = rdwhbu(100))+
    theme_cowplot(font_size = 28)+
    NoAxes()+
    theme(legend.position = "top",legend.justification = "center")
  return(p)
})
wrap_plots(plot.list,ncol = 2)
ggsave(filename = myFileName(prefix = "res/fig/Fig1_monocle_DDRTree_cluster_features",
                             suffix = ".png"),
       width = 8,height = 12,dpi = 350)



# tmp.cds.subset <- tmp.cds[tmp.markers,]
# 
# plot_genes_in_pseudotime(tmp.cds.subset,
#                          color_by = "Stage",
#                          trend_formula = "~ sm.ns(Pseudotime,df=15)")
# 
# plot_pseudotime_heatmap(tmp.cds.subset,
#                         show_rownames = T,
#                         hmcols = coldwarm(100),
#                         cluster_rows = F)


#####--------3.3.3.2 plot pheatmap-----------------

tmp.cell.id <- pData(tmp.cds) %>%
  group_by(Stage) %>%
  mutate(tmp.rank = row_number(Pseudotime)) %>%
  arrange(Stage,tmp.rank) %>%
  pull(cell_id)

tmp.mat <- beforeEPI.exp[tmp.markers,tmp.cell.id]
tmp.annotation <- pData(tmp.cds) %>%
  select(cell_type,Stage,Pseudotime)


tmp.cell.id <- pData(tmp.cds) %>%
  mutate(tmp.rank = row_number(Pseudotime)) %>%
  arrange(tmp.rank) %>%
  pull(cell_id)

tmp.mat <- beforeEPI.exp[tmp.markers,tmp.cell.id]

tmp.levels.1 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst")
tmp.Stage.color <- col.spectral(11)
names(tmp.Stage.color) <- tmp.levels.1



tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "ICM","TE")
tmp.cell_type.color <- divergentcolor(10)
names(tmp.cell_type.color) <- tmp.levels.2

ann_colors = list(
  Pseudotime = col.spectral(100),
  cell_type = tmp.cell_type.color,
  Stage = tmp.Stage.color
)

ph <- pheatmap_fixed(mat = tmp.mat,
               annotation_col = tmp.annotation,
               annotation_colors = ann_colors,
               cluster_rows = F,
               name = "Exp",
               fontsize = 18,
               cluster_cols = F,
               show_colnames = F, 
               color = coldwarm(100))

png(filename = myFileName(prefix = "res/fig/Fig1_show_pseudotime_pheatmaap",
                          suffix = ".png"),
    width = 12,height = 8,units = "in",res = 350)
ph
dev.off()

pdf(file  = myFileName(prefix = "res/fig/Fig1_show_pseudotime_pheatmaap",
                          suffix = ".pdf"),
    width = 12,height = 8,useDingbats = F)
ph
dev.off()

####-----------3.3.3.3 get variable gene----------------
tmp.cds <- readRDS(file = "res/R/early.scRNAseq.moncole_2022031116.rds")

####load data
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)

beforeEPI.meta <- beforeEPI.meta %>%
  dplyr::filter(Stage!="fibroblast") 
beforeEPI.exp <- beforeEPI.exp[,beforeEPI.meta$cell_id]

tmp.cell.id <- pData(tmp.cds) %>%
  mutate(tmp.rank = row_number(Pseudotime)) %>%
  arrange(tmp.rank,Stage) %>%
  pull(cell_id)

tmp.gene <- tmp.cds@featureData@data %>%
  dplyr::filter(use_for_ordering == TRUE) %>%
  pull(gene_short_name)


tmp_labels <-  c("Nlrp5","Zscan4f","Zscan4c","Gata3","Cdx2","Sox2","Nanog","Klf4","Gata6")
pt.matrix <- beforeEPI.exp[tmp_labels,tmp.cell.id]
x <- as.numeric(pt.matrix["Zscan4f",])
plot(x)
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#or exprs(cds)
pt.matrix <- as.matrix(pt.matrix)
tmp.colnames <- colnames(pt.matrix)

# smooth the matrix
# smooth and calculate zscore
tmp.mat <- t(apply(pt.matrix,1,function(x){smooth.spline(x)$y}))
x <- tmp.mat["Zscan4f",]
plot(x)
tmp.mat <- t(apply(tmp.mat,1,function(x){(x-mean(x))/sd(x)}))
plot(x)

rownames(tmp.mat) <- tmp_labels
colnames(tmp.mat) <- tmp.colnames


range(tmp.mat)

scale.max <- 3
tmp.mat[tmp.mat > scale.max] <- scale.max
tmp.mat[tmp.mat < -scale.max] <- -scale.max

range(tmp.mat)
head(pData(tmp.cds))
tmp.colannotation <- pData(tmp.cds) %>%
  select(Stage,Pseudotime)
head(tmp.colannotation)

tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst")
tmp.stage.color <- col.spectral(length(tmp.levels.2))
names(tmp.stage.color) <- tmp.levels.2
ann_colors = list(
  Pseudotime = warm(100),
  Stage =  tmp.stage.color)

tmp_labels <-  c("Nlrp5","Zscan4f","Zscan4c","Gata3","Cdx2","Sox2","Nanog","Klf4","Gata6")
tmp_idx <- which(rownames(tmp.mat) %in% tmp_labels)
rownames(tmp.mat)[tmp_idx]

# tmp_anno <- anno_mark(at=tmp_idx,
#                       labels=rownames(tmp.mat)[tmp_idx],
#                       which = "row",
#                       labels_gp = gpar(fontsize=18,
#                                        fontface="italic"))

tmp_labels <- c("Nlrp5","Zscan4f","Zscan4c","Klf4","Gata6","Cdx2","Nanog","Sox2","Gata3")
setdiff(rownames(tmp.mat),tmp_labels)
setdiff(tmp_labels,rownames(tmp.mat))
ph <- pheatmap_fixed(mat = tmp.mat[tmp_labels,],
                     color = coldwarm(100),
                     show_rownames = T,
                     na_col = "black",
                     annotation_colors = ann_colors,
                     show_colnames = F,
                     fontsize = 16,
                     cluster_cols = F,
                     cluster_rows = F,
                     name = "zscore",
                     annotation_col = tmp.colannotation,
                     clustering_method = "ward.D2")


ph

png(filename = myFileName(prefix = "res/fig/fig1_variable_gene_pseudotime",
                          suffix = ".png"),
    width = 6,height = 8,units = "in",res = 300)
ph
dev.off()

pdf(file = myFileName(prefix = "res/fig/fig1_variable_gene_pseudotime",
                      suffix = ".pdf"),
    width = 6,height = 8)
ph
dev.off()




####----------3.3.4 bulk and pseudo-bulk correlation-------------

####get pseudo-bulk
beforeEPI.exp <- seu@assays$RNA@data
beforeEPI.meta <- seu[[]]
tmp.level <- c("MII oocyte","zygote",
               "early 2-cell","mid 2-cell","late 2-cell",
               "4-cell","8-cell","16-cell",
               "ICM","TE")
### test if normalized
# head(beforeEPI.meta)
# boxplot(beforeEPI.exp[,110:100])
beforeEPI.exp.mean <- beforeEPI.exp %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  gather(key = "cell_id",value = "gene.exp",-Gene) %>%
  left_join(beforeEPI.meta,"cell_id") %>%
  group_by(Gene,cell_type) %>%
  summarise(gene.exp.mean=mean(gene.exp)) %>%
  ungroup() %>%
  mutate(cell_type=factor(cell_type,levels = tmp.level)) %>%
  spread(key="cell_type",value = "gene.exp.mean") %>%
  column_to_rownames("Gene")

colnames(beforeEPI.exp.mean) <- c("MII_oocyte","zygote","early2cell","mid2cell","late2cell","4cell","8cell","16cell","ICM","TE")


saveRDS(object = beforeEPI.exp.mean,
        file = myFileName(prefix = "res/R/beforeEPI_exp_mean",suffix = ".rds"))

##### load Prof xiewei data
##### RNA-seq xiewei data
xiewei.data.exp <- readRDS(file = "res/R/xiewei.data.exp_20201012.rds")
xiewei.data.exp <- myNormalize(myRemoveNA(xiewei.data.exp))

colnames(xiewei.data.exp)[3] <- "early2cell"
colnames(xiewei.data.exp)[4] <- "late2cell"
#### overlap with gene symbols
gene_symbols <- intersect(rownames(xiewei.data.exp),
                          rownames(beforeEPI.exp.mean))




####--------3.3.4.1 scatter plot-----------

Stage <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell","ICM")
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
         suffix = ".png",
         width = 16,height = 16,dpi=350)

myggsave(pp,prefix = "res/fig/figS2_QC_beforeEPI_cor_revised",
         suffix = ".pdf",
         width = 16,height = 16,dpi=350)


#####----3.3.4.2 Fig 1B Correlation pheatmap-----

Stage <- c("MII_oocyte","zygote","early2cell","late2cell","4cell","8cell","ICM")
data.plot.1 <- beforeEPI.exp.mean[gene_symbols,]
data.plot.2 <- xiewei.data.exp[gene_symbols,]


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

tmp.stage  <- c("MII oocyte","zygote","early 2-cell","late 2-cell","4-cell","8-cell","ICM")

png(filename = myFileName(prefix = "res/fig/Fig1_QC_pheatmap",suffix = ".png"),
    width = 8,height = 8,res = 350,units = "in")
exprTable_t <- as.data.frame(t(res.mat))
col_dist = dist(exprTable_t)
hclust_1 <- hclust(col_dist)
manual_order = rownames(res.mat)
#dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))))
dend = reorder(as.dendrogram(hclust_1), wts=order(match(manual_order, rownames(exprTable_t))), agglo.FUN = max)
col_cluster <- as.hclust(dend)
ph <- ComplexHeatmap::pheatmap(res.mat,fontsize = 20,
                         show_rownames = T,
                         display_numbers = T,
                         number_color = "white",
                         cluster_rows = col_cluster,
                         cluster_cols = col_cluster,
                         fontsize_number = 22,
                         name = "PCC",border_color = "black",
                         labels_row = paste0("bulk ",tmp.stage),
                         labels_col = paste0("pseudo-bulk ",tmp.stage),
                         show_colnames = T,
                         heatmap_legend_param = list(title_gp = gpar(fontsize = 20,fontface = "bold"),
                                                     labels_gp=gpar(fontsize = 20)),
                         angle_col = "315",
                         color = rdbu(100))
dev.off()



pdf(file  = myFileName(prefix = "res/fig/Fig1_QC_pheatmap",suffix = ".pdf"),
    width = 8,height = 8)
ph
dev.off()










