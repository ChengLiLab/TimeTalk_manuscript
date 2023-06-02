
####TimeTalk manuscript
####Author:wlt
####Fig4
####This figure illustrate the mechanism revealed by TimeTalk
####start at: 20220411
####tidy code at: 20220518
####Final code at: 20220618;
####Remove Epigenetic part
####edit in 20230329
####Fix bugs caused by KEGG update
####Please see https://github.com/YuLab-SMU/clusterProfiler/issues/561

####----keep you alive---------
for(ii in 1:100000){
  cat(ii,sep = "\n")
  Sys.sleep(time = ii)
}
####----0. load package--------
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
library(lmtest)
library(ggrepel)
library(monocle)
#library(netSmooth)
library(clusterProfiler)
### local KEGG db
library(KEGG.db)
library(org.Mm.eg.db)
library(entropy)
library(RTN)
library(RTNduals)
library(snow)
library(ComplexHeatmap)
library(msigdbr)
library(patchwork)
library(universalmotif)
library(GenomicFeatures)
library(rtracklayer)
library(RcisTarget)
library(GSVA)
#library(metaseqR)
#library(CellChat)
source(file = "code/myUtils.R")

####--------1.1 load data------------
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)
beforeEPI.cds <- readRDS(file = "res/R/early.scRNAseq.moncole_2022031116.rds")
beforeEPI.meta <- pData(beforeEPI.cds)
tmp.cell.id <- beforeEPI.meta %>%
  group_by(Stage) %>%
  mutate(tmp.rank = row_number(Pseudotime)) %>%
  arrange(Stage,tmp.rank) %>%
  pull(cell_id)

eLR.df <- readRDS(file = "res/R/putative_eLR_pairs_1094_2022032109.rds") %>%
  filter(abs(PCC) > 0.1) %>%
  #filter(LRtoTF < 0.01 | TFtoLR < 0.01) %>%
  mutate(group = "background") %>%
  mutate(group = ifelse(LRtoTF < 0.01 & TFtoLR > 0.01,"forward",group)) %>%
  mutate(group = ifelse(LRtoTF > 0.01 & TFtoLR < 0.01,"backward",group)) %>%
  mutate(group = ifelse(LRtoTF < 0.01 & TFtoLR < 0.01,"feedback",group)) %>%
  mutate(group = factor(group,levels = c("forward","backward","feedback","background"))) %>%
  mutate(LRtoTF=-log10(LRtoTF)) %>%
  mutate(TFtoLR=-log10(TFtoLR))

tmp.master.regulon <- readRDS(file = "res/R/putative_master_regulon_2022032109.rds")

####-------1.2 test curve plot---------------
#### needn't run,
####20220518
ii <- "forward"

for(ii in c("forward","backward","feedback","background")){
  Lgene.list <- eLR.df %>%
    filter(group == ii) %>%
    pull(Lgene)
  Rgene.list <- eLR.df %>%
    filter(group == ii) %>%
    pull(Rgene)
  tmp.L.exp <- beforeEPI.exp[Lgene.list,tmp.cell.id] 
  tmp.R.exp <- beforeEPI.exp[Rgene.list,tmp.cell.id]
  tmp.TF.exp.score <- colMeans(beforeEPI.exp[tmp.master.regulon$Regulon,tmp.cell.id])
  head(tmp.TF.exp.score)
  tmp.IS <- sqrt(tmp.L.exp*tmp.R.exp)
  rownames(tmp.IS) <- paste0(Lgene.list,"-",Rgene.list)
  #rownames(tmp.TF.exp.score) <- "regulon"
  tmp.IS <- rbind(tmp.IS,regulon=tmp.TF.exp.score)
  
  tmp.meta.df <- beforeEPI.meta %>%
    group_by(Stage) %>%
    mutate(tmp.rank = row_number(Pseudotime)) %>%
    arrange(Stage,tmp.rank) %>%
    dplyr::select(cell_id,Stage,tmp.rank) %>%
    ungroup() %>%
    mutate(tmp.rank = row_number()) 
  
  
  tmp.IS.df <- tmp.IS %>%
    rownames_to_column("LRpairs") %>%
    gather(key = "cell_id",value = "IS",-LRpairs) %>%
    left_join(tmp.meta.df,by = "cell_id") %>%
    mutate(tmp.class=ifelse(LRpairs=="regulon","regulon","LR"))
  
  p1 <- ggplot()+
    geom_smooth(data = tmp.IS.df,
                se=F,
                method = "gam", 
                formula = y ~ s(x, bs = "cs"),
                mapping = aes(x = tmp.rank,y = IS,group=LRpairs,color=tmp.class))+
    xlab("rank")+
    ylab("IS")+
    scale_x_continuous(expand = c(0,0))+
    scale_color_manual(name=NULL,values = divergentcolor(2))+
    ggtitle(ii)+
    theme_cowplot(font_size = 28,rel_small = 9/14)+
    theme(axis.line = element_line(size = 2),
          legend.position = "top",
          legend.justification = "center",
          axis.ticks = element_line(size = 2),
          plot.title = element_text(hjust = 0.5))
  p1
  p2 <- ggplot()+
    geom_smooth(data = tmp.IS.df,
                se=F,
                method = "gam", 
                formula = y ~ s(x, bs = "cs"),
                mapping = aes(x = tmp.rank,y = IS,color=tmp.class))+
    xlab("rank")+
    ylab("IS")+
    scale_x_continuous(expand = c(0,0))+
    scale_color_manual(name=NULL,values = divergentcolor(2))+
    ggtitle(ii)+
    theme_cowplot(font_size = 28,rel_small = 9/14)+
    theme(axis.line = element_line(size = 2),
          
          legend.position = "top",
          legend.justification = "center",
          axis.ticks = element_line(size = 2),
          plot.title = element_text(hjust = 0.5))
  
  p <- p1 | p2
  myggsave(p = p,prefix = paste0("res/fig/fig2_show_",ii,"_correlation"),suffix = ".png",width = 8,height = 6,dpi=350)
  myggsave(p = p,prefix = paste0("res/fig/fig2_show_",ii,"_correlation"),suffix = ".pdf",width = 8,height = 6,dpi=350)
  
}

####---------1.2  show regulon expression----------------------
head(tmp.master.regulon)
tmp.gene <- tmp.master.regulon$Regulon
tmp.exp <- beforeEPI.exp[tmp.gene,tmp.cell.id]

tmp.meta.df <- beforeEPI.meta %>%
  group_by(Stage) %>%
  mutate(tmp.rank = row_number(Pseudotime)) %>%
  arrange(Stage,tmp.rank) %>%
  dplyr::select(cell_id,Stage,Pseudotime,tmp.rank) %>%
  ungroup() %>%
  mutate(tmp.rank = row_number()) %>%
  dplyr::select(cell_id,Stage,Pseudotime) %>%
  column_to_rownames("cell_id")

tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst")
tmp.Stage.color <- col.spectral(length(tmp.levels.2))
names(tmp.Stage.color) <- tmp.levels.2

ann_colors = list(
  Pseudotime = col.spectral(100),
  Stage = tmp.Stage.color
)


p1 <- pheatmap_fixed(tmp.exp[ c("Pou5f1","Zscan4f","Zscan4d","Sox2","Klf5","Cdx2","Tead4"),],
                     scale = "none",color = coolwarm(100),
                     show_colnames  = F,
                     cluster_rows = T,
                     show_rownames = T,
                     name = "Exp",
                     fontsize = 18,
                     annotation_col = tmp.meta.df,
                     annotation_colors = ann_colors,
                     cluster_cols = F)
p2 <- pheatmap_fixed(tmp.exp,
                     scale = "none",color = coolwarm(100),
                     name = "Exp",
                     show_rownames = F,
                     show_colnames = F,
                     cluster_rows = T,
                     cluster_cols = F,
                     fontsize = 18,
                     annotation_colors = ann_colors,
                     annotation_col = tmp.meta.df)
png(filename = myFileName(prefix = "res/fig/heatmap_show_regulon_case",suffix = ".png"),
    width = 8,height = 6,res = 350,units = "in")
p1
dev.off()

png(filename = myFileName(prefix = "res/fig/heatmap_show_regulon_al",suffix = ".png"),
    width = 8,height = 6,res = 350,units = "in")
p2
dev.off()



####---------1.3 filter to get the master regulator of the result-------------

tmp.master.regulon <- readRDS(file = "res/R/putative_master_regulon_2022032109.rds")
head(tmp.master.regulon)
### add mean and sd to data frame
tmp.TF.mean <- rowMeans(beforeEPI.exp[tmp.master.regulon$Regulon,])
tmp.TF.sd <- apply(beforeEPI.exp[tmp.master.regulon$Regulon,],MARGIN = 1,sd)
tmp.master.regulon$mean.exp <- tmp.TF.mean
tmp.master.regulon$sd <- tmp.TF.sd
tmp.df <- tmp.master.regulon %>%
  filter(mean.exp > 1)

####---------1.4 plot the master regulator heatmap -----------

####---------1.4.1 prepare data------------
tmp.min <- -3
tmp.max <- 3
tmp.mat <- beforeEPI.exp[tmp.df$Regulon,tmp.cell.id]
tmp.mat <- pheatmap:::scale_rows(tmp.mat)
tmp.mat[tmp.mat > tmp.max] <- tmp.max
tmp.mat[tmp.mat < tmp.min] <- -tmp.min

####-------1.4.2 draw draft pheatmap----------

####add anno mark
tmp_labels <-  c("Pou5f1","Zscan4f","Sox2","Cdx2","Tead4","Klf6")
tmp_idx <- which(rownames(tmp.mat) %in% tmp_labels)
rownames(tmp.mat)[tmp_idx]
tmp_anno <- anno_mark(at=tmp_idx,
                      labels=rownames(tmp.mat)[tmp_idx],
                      which = "row",
                      labels_gp = gpar(fontsize=18,fontface="italic"))

#### add column information
tmp.meta.df.col <- beforeEPI.meta %>%
  group_by(Stage) %>%
  mutate(tmp.rank = row_number(Pseudotime)) %>%
  arrange(Stage,tmp.rank) %>%
  dplyr::select(cell_id,Stage,Pseudotime,tmp.rank) %>%
  ungroup() %>%
  mutate(tmp.rank = row_number()) %>%
  dplyr::select(cell_id,Stage,Pseudotime) %>%
  column_to_rownames("cell_id")

tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst")
tmp.Stage.color <- col.spectral(length(tmp.levels.2))
names(tmp.Stage.color) <- tmp.levels.2

ann_colors = list(
  Pseudotime = col.spectral(100),
  Stage = tmp.Stage.color
)

ph <- pheatmap_fixed(tmp.mat,
                     scale = "none",
                     show_rownames = F,
                     show_colnames = F,
                     clustering_method =  "ward.D2",
                     cluster_rows = T,
                     cluster_cols = F,
                     name = "zscore",
                     annotation_col = tmp.meta.df.col,
                     annotation_colors = ann_colors,
                     fontsize = 18,
                     color = rdbu(100))
ph
ph <- draw(ph)
ph.row <- row_dend(ph)
ph.row <- as.hclust(ph.row)


num_clusters <- 6
annotation_row <- data.frame(Cluster = factor(cutree(ph.row,num_clusters)))
tmp.cluster.color <- pal_npg()(num_clusters)
names(tmp.cluster.color) <- 1:num_clusters

ann_colors = list(
  Pseudotime = col.spectral(100),
  Stage = tmp.Stage.color,
  Cluster=tmp.cluster.color
)

#####------1.4.3 re annotation levels-----------
ph <- pheatmap_fixed(tmp.mat,
                     scale = "none",
                     show_rownames = F,
                     show_colnames = F,
                     clustering_method =  "ward.D2",
                     cluster_rows = T,
                     cluster_cols = F,
                     annotation_row = annotation_row,
                     name = "zscore",
                     annotation_col = tmp.meta.df.col,
                     annotation_colors = ann_colors,
                     fontsize = 18,
                     color = rdbu(100))
ph <- draw(ph)
tmp.order <- row_order(ph)


annotation_row$Cluster <- plyr::mapvalues(annotation_row$Cluster,
                                          from = c(3,1,6,4,2,5),
                                          to = 1:num_clusters)
annotation_row$Cluster <- factor(annotation_row$Cluster,levels = 1:num_clusters)
annotation_row.df <- annotation_row %>%
  rownames_to_column("gene")
annotation_row.df <- annotation_row.df[tmp.order,] %>%
  group_by(Cluster) %>%
  mutate(tmp.rank = row_number()) %>%
  ungroup() %>%
  arrange(Cluster,tmp.rank) %>%
  dplyr::select(gene,Cluster) %>%
  column_to_rownames("gene")
tmp.mat.plot <- tmp.mat[rownames(annotation_row.df),]

early_embryo_TF <- annotation_row.df %>%
  rownames_to_column("cluster")

early_embryo_TF %>%
  filter(Cluster == 4)


tmp_labels <-  c("Atf2","Nfya","Nanog","Ctcf",
                 "Obox3","Obox5","Obox6",
                 "Yy1","Pou5f1","Zscan4f",
                 "Sox2","Sox21","Cdx2",
                 "Gata4","Gata6","Gata3",
                 "Tead1","Tead4","Klf6")
setdiff(tmp_labels,list.files("database/CISBP"))
setdiff(list.files("database/CISBP"),tmp_labels)


tmp_idx <- which(rownames(tmp.mat.plot) %in% tmp_labels)
rownames(tmp.mat.plot)[tmp_idx]
tmp_anno <- anno_mark(at=tmp_idx,
                      labels=rownames(tmp.mat.plot)[tmp_idx],
                      which = "row",
                      labels_gp = gpar(fontsize=18,
                                       fontface="italic"))

ph <- pheatmap_fixed(tmp.mat.plot,
                     scale = "none",
                     show_rownames = F,
                     show_colnames = F,
                     clustering_method =  "ward.D2",
                     cluster_rows = F,
                     cluster_cols = F,
                     annotation_row = annotation_row.df,
                     name = "zscore",
                     annotation_col = tmp.meta.df.col,
                     annotation_colors = ann_colors,
                     fontsize = 18,
                     color = rdbu(100))+
  rowAnnotation(mark = tmp_anno)



ph


png(filename = myFileName(prefix = "res/fig/fig2_regulon_TF",suffix = ".png"),
    width = 8,height = 6,units = "in",res = 350)
ph
dev.off()

pdf(file = myFileName(prefix = "res/fig/fig2_regulon_TF",suffix = ".pdf"),
    width = 8,height = 6)
ph
dev.off()


####-----------1.4.4 show cluster TF temporal pattern -----------------

tmp.TF.meta <- annotation_row.df %>%
  rownames_to_column("gene")
tmp.TF.meta$Cluster <- paste0("C",tmp.TF.meta$Cluster)
write.table(tmp.TF.meta,file = "res/txt/tTFs_list.txt",quote = F,sep = "\t",row.names = F)



plot.list <- lapply(1:num_clusters,FUN = function(ii){
  tmp.gene <- tmp.TF.meta %>%
    dplyr::filter(Cluster==ii) %>%
    pull(gene)
  
  tmp.gene.exp <- tmp.mat.plot[tmp.gene,]
  
  tmp.meta.df <- beforeEPI.meta %>%
    group_by(Stage) %>%
    mutate(tmp.rank = row_number(Pseudotime)) %>%
    arrange(Stage,tmp.rank) %>%
    dplyr::select(cell_id,Stage,tmp.rank) %>%
    ungroup() %>%
    mutate(tmp.rank = row_number()) 
  
  
  tmp.plot.df <- tmp.gene.exp %>%
    rownames_to_column("gene") %>%
    gather(key = "cell_id",value = "zscore",-gene) %>%
    left_join(tmp.meta.df,by = "cell_id")
  
  p1 <- ggplot()+
    geom_smooth(data = tmp.plot.df,
                se=F,
                method = "gam", 
                formula = y ~ s(x, bs = "cs"),
                mapping = aes(x = tmp.rank,y = zscore,group=gene),
                color = divergentcolor(num_clusters)[ii])+
    xlab("rank")+
    ylab("zscore")+
    scale_x_continuous(expand = c(0,0))+
    ggtitle(paste0("cluster ",ii))+
    theme_cowplot(font_size = 28,rel_small = 9/14)+
    theme(axis.line = element_line(size = 2),
          legend.position = "top",
          legend.justification = "center",
          axis.ticks = element_line(size = 2),
          plot.title = element_text(hjust = 0.5))
  
  return(p1)
})
p <- wrap_plots(plot.list,nrow = 2)
p
myggsave(p,prefix = "res/fig/fig2_TF_cluster",suffix = ".png",width = 12,height = 9,dpi=350)
saveRDS(annotation_row.df,file = myFileName(prefix = "res/R/early_embryo_core_TF_annotation_row",suffix = ".rds"))

#as far as I am concerned, Heatmap is suitable;
#it is better to draw Heatmap

####----------1.4.5 invest the function of the TF------------------

####----------1.4.5.1 GO--------------
tmp.TF.annotation <- readRDS(file = "res/R/early_embryo_core_TF_annotation_row_2022032521.rds")
tmp.df.use <- tmp.TF.annotation %>%
  rownames_to_column("TF")
tmp.gene.list <- lapply(1:6,
                        FUN = function(ii){
                          tmp.gene <- tmp.df.use %>%
                            filter(Cluster == ii) %>%
                            pull(TF)
                          eg = bitr(tmp.gene, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)
                          tmp.gene <- eg$ENTREZID
                          return(tmp.gene)
                        })
names(tmp.gene.list) <- paste0("C",1:6)

Compare <- compareCluster(geneCluster=tmp.gene.list,
                          fun = "enrichGO",
                          OrgDb = org.Mm.eg.db,
                          readable = T,
                          ont = "BP",
                          pvalueCutoff=0.05)

p <- dotplot(Compare,showCategory=20,font.size = 28)+
  scale_color_gradientn(colors  = coolwarm(100))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
tmp.df <- Compare@compareClusterResult

myggsave(p = p,prefix = "res/fig/fig4_TF_GO_test",
         suffix = ".png",
         width = 16,
         height = 12,
         dpi=350)

####---------1.4.5.2 GO with regulon---------------

####load data
tmp.TF.annotation <- readRDS(file = "res/R/early_embryo_core_TF_annotation_row_2022032521.rds")
tmp.df.use <- tmp.TF.annotation %>%
  rownames_to_column("TF")
head(tmp.df.use)

rtna <- readRDS(file = "res/R/early.rtna_20220316.rds")
tmp.list <- RTN::tna.get(rtna,what = "regulons")



tmp.gene.list <- lapply(1:6,FUN = function(ii){
                          tmp.gene <- tmp.df.use %>%
                            filter(Cluster == ii) %>%
                            pull(TF)
                          tmp.res <- lapply(tmp.gene, function(xx){
                            tmp.list[[xx]]
                          })
                          tmp.gene <- unique(unlist(tmp.res))
                          
                          eg = bitr(tmp.gene, 
                                    fromType="SYMBOL",
                                    toType="ENTREZID", 
                                    OrgDb= org.Mm.eg.db)
                          
                          tmp.gene <- eg$ENTREZID
                          return(tmp.gene)
                        })
names(tmp.gene.list) <- paste0("C",1:6)


Compare <- compareCluster(geneCluster=tmp.gene.list,
                          fun = "enrichGO",
                          OrgDb = org.Mm.eg.db,
                          readable = T,
                          ont = "BP",
                          pvalueCutoff=1)



p <- dotplot(Compare,showCategory=20,font.size = 28)+
  scale_color_gradientn(colors  = coolwarm(100))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))

tmp.term <- c("oogenesis","germ cell development","stem cell population maintenance","maintenance of cell number","negative regulation of cell cycle",
              "negative regulation of DNA binding","cell fate commitment","stem cell differentiation",
              "regulation of gene expression, epigenetic","regulation of epithelial cell differentiation","embryonic organ morphogenesis",
              "histone methylation","protein methylation",
              "chromatin silencing","chromatin assembly or disassembly","DNA conformation change",
              "somatic stem cell population maintenance","embryonic organ development","trophectodermal cell differentiation","blastocyst development","blastocyst formation","blastocyst growth")


tmp.df <- Compare@compareClusterResult 

####sapply is good!
tmp.df <- tmp.df %>%
  dplyr::filter(Description %in% tmp.term)  %>%
  mutate(GeneRatio = sapply(GeneRatio,
                            FUN = function(xx) eval(parse(text = xx)),
                            USE.NAMES = F)) %>%
  mutate(Description = factor(Description,levels = rev(tmp.term)))
####head(tmp.df)



p <- ggplot(data = tmp.df,aes(Cluster, Description,size = GeneRatio,color=pvalue))+
  geom_point() +
  scale_color_gradientn(colors = rev(coolwarm(100)),limits = c(NA,0.05))+
  myGridTheme(font.size = 28)+
  ylab(NULL)

p

myggsave(p = p,prefix = "res/fig/fig4_TF_GO_enrichment_term",suffix = ".png",
         width = 14,height = 8,dpi=350)


####------------2. Select eLR by correlation with TF------------------------------------

####------------2.1 load data------------------

LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_2022031117.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)

beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)


beforeEPI.cds <- readRDS(file = "res/R/early.scRNAseq.moncole_2022031116.rds")
beforeEPI.meta <- pData(beforeEPI.cds)
tmp.cell.id <- beforeEPI.meta %>%
  group_by(Stage) %>%
  mutate(tmp.rank = row_number(Pseudotime)) %>%
  arrange(Stage,tmp.rank) %>%
  pull(cell_id)

tmp.TF.meta <- readRDS(file = "res/R/early_embryo_core_TF_annotation_row_2022032521.rds") %>%
  rownames_to_column("gene")




####-------2.2 perform granger causal test -----------

#### Actualy the following code is the TimeTalk version 0.0.1
#### This code will generate 430 eLR

num_clusters <- length(levels(tmp.TF.meta$Cluster))
tmp.res.list <- lapply(1:num_clusters,FUN = function(ii){
  cat(ii,sep = "\n")
  tmp.L.exp <- beforeEPI.exp[Lgene.list,tmp.cell.id] 
  tmp.R.exp <- beforeEPI.exp[Rgene.list,tmp.cell.id]
  tmp.res.PCC <- sapply(1:nrow(tmp.L.exp), function(i) cor(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),
                                                           method = "pearson"))
  tmp.res.SCC <- sapply(1:nrow(tmp.L.exp), function(i) cor(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),
                                                           method = "spearman"))
  tmp.res.rate.L <- sapply(1:nrow(tmp.L.exp), function(i) sum(as.numeric(tmp.L.exp[i,])>0)/length(as.numeric(tmp.L.exp[i,])>0))
  tmp.res.rate.R <- sapply(1:nrow(tmp.R.exp), function(i) sum(as.numeric(tmp.R.exp[i,])>0)/length(as.numeric(tmp.R.exp[i,])>0))
  
  
  tmp.TF.gene <- tmp.TF.meta %>%
    filter(Cluster == ii) %>%
    pull(gene)
  
  tmp.TF.exp.score <- colMeans(beforeEPI.exp[tmp.TF.gene,tmp.cell.id])
  
  tmp.res.granger.LRtoTF <- sapply(1:nrow(tmp.L.exp), function(i) {
    tmp.res.res <-  tryCatch(
      expr = {
        tmp.IS.vec <- as.numeric(tmp.L.exp[i,]) * as.numeric(tmp.R.exp[i,])
        tmp.IS.vec <- sqrt(tmp.IS.vec)
        tmp.res <- grangertest(tmp.IS.vec, tmp.TF.exp.score)
        return(tmp.res$`Pr(>F)`[2])},
      error = function(e) {
        return(1)
      })
    return(tmp.res.res)
  })
  
  tmp.res.granger.TFtoLR <- sapply(1:nrow(tmp.L.exp), function(i) {
    tmp.res.res <-  tryCatch(
      expr = {
        tmp.IS.vec <- as.numeric(tmp.L.exp[i,]) * as.numeric(tmp.R.exp[i,])
        tmp.IS.vec <- sqrt(tmp.IS.vec)
        tmp.res <- grangertest(tmp.TF.exp.score,tmp.IS.vec)
        return(tmp.res$`Pr(>F)`[2])},
      error = function(e) {
        return(1)
      })
    return(tmp.res.res)
  })
  
  
  tmp.res.cor.1 <- data.frame(Lgene=Lgene.list,
                              Rgene=Rgene.list,
                              LRpairs=paste0(Lgene.list,"_",Rgene.list),
                              PCC=tmp.res.PCC,
                              SCC=tmp.res.SCC,
                              detection.rate.L=tmp.res.rate.L,
                              detection.rate.R=tmp.res.rate.R,
                              LRtoTF=tmp.res.granger.LRtoTF,
                              TFtoLR=tmp.res.granger.TFtoLR,
                              TFcluster=paste0("TF_cluster_",ii),
                              stringsAsFactors = F) %>%
    dplyr::filter(detection.rate.L > 0.05 & detection.rate.R > 0.05) %>%
    arrange(-PCC)
  
  
  tmp.res.eLR <- tmp.res.cor.1 %>% 
    filter(abs(PCC) > 0.1) %>%
    filter(abs(PCC) > 0.1) %>%
    #filter(LRtoTF < 0.01 | TFtoLR < 0.01) %>%
    mutate(group = "background") %>%
    mutate(group = ifelse(LRtoTF < 0.01 & TFtoLR > 0.01,"forward",group)) %>%
    mutate(group = ifelse(LRtoTF > 0.01 & TFtoLR < 0.01,"backward",group)) %>%
    mutate(group = ifelse(LRtoTF < 0.01 & TFtoLR < 0.01,"feedback",group)) %>%
    mutate(group = factor(group,levels = c("forward","backward","feedback","background"))) 
  
  return(tmp.res.eLR)
})

tmp.res.eLR.df <- Reduce(rbind,tmp.res.list)
saveRDS(object = tmp.res.eLR.df,
        file = myFileName(prefix = "res/R/early_embryo_eLR_TF_interaction_information",
                          suffix = ".rds"))

#####-------2.3 analysis result, and draw heatmap -----------------

tmp.res.eLR.df <- Reduce(rbind,tmp.res.list)
tmp.res.eLR.annotation <- tmp.res.eLR.df %>%
  mutate(LRpairs = paste0(Lgene,"-",Rgene)) %>%
  dplyr::filter(group!="background") %>%
  dplyr::select(LRpairs,TFcluster,group,Lgene,Rgene) %>%
  spread(key = TFcluster,value = group) %>%
  column_to_rownames("LRpairs")

tmp.res.eLR.annotation[is.na(tmp.res.eLR.annotation)] <- "background"

tmp.L.exp <- beforeEPI.exp[tmp.res.eLR.annotation$Lgene,tmp.cell.id] 
tmp.R.exp <- beforeEPI.exp[tmp.res.eLR.annotation$Rgene,tmp.cell.id]
tmp.IS.mat <- sqrt(tmp.L.exp*tmp.R.exp)
rownames(tmp.IS.mat) <- paste0(tmp.res.eLR.annotation$Lgene,"-",tmp.res.eLR.annotation$Rgene)


annotation_row.df <- tmp.res.eLR.annotation %>%
  dplyr::select(-Lgene,-Rgene)
tmp.IS.mat.plot <- pheatmap:::scale_rows(tmp.IS.mat)
range(tmp.IS.mat.plot)

tmp.meta.df.col <- beforeEPI.meta %>%
  group_by(Stage) %>%
  mutate(tmp.rank = row_number(Pseudotime)) %>%
  arrange(Stage,tmp.rank) %>%
  dplyr::select(cell_id,Stage,Pseudotime,tmp.rank) %>%
  ungroup() %>%
  mutate(tmp.rank = row_number()) %>%
  dplyr::select(cell_id,Stage,Pseudotime) %>%
  column_to_rownames("cell_id")


tmp.res.TF.score <-  lapply(1:6, function(ii){
  tmp.TF.gene <- tmp.TF.meta %>%
    filter(Cluster == ii) %>%
    pull(gene)
  tmp.TF.exp.score <- colMeans(beforeEPI.exp[tmp.TF.gene,tmp.cell.id])
  return(tmp.TF.exp.score)
})
tmp.res.TF.score <- Reduce(cbind,tmp.res.TF.score)
colnames(tmp.res.TF.score) <- paste0("C",1:6)
tmp.meta.df.col <- cbind(tmp.meta.df.col,tmp.res.TF.score)



scale.max <- 3
scale.min <- -3

tmp.IS.mat.plot[tmp.IS.mat.plot < scale.min] <- scale.min
tmp.IS.mat.plot[tmp.IS.mat.plot > scale.max] <- scale.max

tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst")
tmp.Stage.color <- col.spectral(length(tmp.levels.2))
names(tmp.Stage.color) <- tmp.levels.2

tmp.color <- c("brown3","navy","#f2be58","grey")
names(tmp.color) <-  c("forward","backward","feedback","background")
Cluster <- pal_npg()(6)
names(Cluster) <- 1:6
ann_colors = list(
  TF_cluster_1 = tmp.color,
  TF_cluster_2 = tmp.color,
  TF_cluster_3 = tmp.color,
  TF_cluster_4 = tmp.color,
  TF_cluster_5 = tmp.color,
  TF_cluster_6 = tmp.color,
  Cluster = Cluster,
  Pseudotime = col.spectral(100),
  Stage = tmp.Stage.color,
  C1=coldwarm(100),
  C2=coldwarm(100),
  C3=coldwarm(100),
  C4=coldwarm(100),
  C5=coldwarm(100),
  C6=coldwarm(100)
)


ph <- pheatmap_fixed(tmp.IS.mat.plot,
                     show_rownames = F,
                     clustering_method = "ward.D2",
                     show_colnames = F,
                     cluster_cols = F,
                     scale = "none",
                     #annotation_colors = ann_colors,
                     color = rdbu(100),
                     annotation_row = annotation_row.df)
ph <- draw(ph) 
ph.row <- row_dend(ph)
ph.row <- as.hclust(ph.row)
ph.row.order <- row_order(ph)

tmp.number.cluster <- 6
tmp.row.cluster <- data.frame(Cluster=factor(cutree(ph.row,tmp.number.cluster)))
annotation_row.df <- cbind(annotation_row.df,tmp.row.cluster)
annotation_row.df$Cluster <- plyr::mapvalues(annotation_row.df$Cluster,from = c(4,3,2,1,5,6),to = 1:6)
annotation_row.df$Cluster <- factor(annotation_row.df$Cluster,1:6)

tmp.final.annotation <- annotation_row.df[ph.row.order,] %>%
  rownames_to_column("LRpairs") %>%
  group_by(Cluster) %>%
  mutate(tmp.rank = row_number()) %>%
  ungroup() %>%
  arrange(Cluster,tmp.rank) %>%
  dplyr::select(-tmp.rank) %>%
  column_to_rownames("LRpairs")


tmp.mat.plot <- tmp.IS.mat.plot[rownames(tmp.final.annotation),]
tmp_labels <-  c("Fgf4-Fgfr1","Fgf4-Fgfr2","Bmp4-Bmpr2",
                 "Gdf9-Bmpr1b","Ntn1-Dcc","Kitl-Kit",
                 "Rspo2-Lgr5","Rspo2-Lgr6",'Hbegf-Erbb4',
                 "Pdgfc-Pdgfra","Igf2-Igf2r",
                 "Ctsd-Lrp1")
tmp_idx <- which(rownames(tmp.mat.plot) %in% tmp_labels)
setdiff(tmp_labels,rownames(tmp.mat.plot)[tmp_idx])

tmp_anno <- anno_mark(at=tmp_idx,
                      labels=rownames(tmp.mat.plot)[tmp_idx],
                      which = "row",
                      labels_gp = gpar(fontsize=18,
                                       fontface="italic"))


ph <- pheatmap_fixed(tmp.mat.plot,
                     show_rownames = F,
                     clustering_method = "ward.D2",
                     show_colnames = F,
                     cluster_rows = F,
                     cluster_cols = F,
                     scale = "none",
                     annotation_colors = ann_colors,
                     color = rdbu(100),
                     annotation_row = tmp.final.annotation,
                     annotation_col = tmp.meta.df.col,
                     fontsize = 16,name = "zscore")+
  rowAnnotation(mark = tmp_anno)



png(filename = myFileName(prefix = "res/fig/fig2_regulon_TF_eLR",suffix = ".png"),
    width = 12,height = 9,units = "in",res = 350)
ph
dev.off()


pdf(file = myFileName(prefix = "res/fig/fig2_regulon_TF_eLR",suffix = ".pdf"),
    width = 12,height = 9)
ph
dev.off()


saveRDS(tmp.final.annotation,file = myFileName(prefix = "res/R/eLR_row_annotaion_df",suffix = ".rds"))
saveRDS(tmp.meta.df.col,file = myFileName(prefix = "res/R/eLR_col_annoation_df_sc_TFscore",suffix = ".rds"))
####
tmp.check <- tmp.final.annotation %>%
  rownames_to_column("LRpairs")

pdf(file = myFileName(prefix = "res/fig/fig2_regulon_TF_eLR",suffix = ".pdf"),
    width = 12,height = 9,useDingbats = F)
ph
dev.off()



####--------2.4 perform GO-------------



####using mysidb C5 for enrichment
Mm_msigdb <- msigdbr(species = "Mus musculus")
Mm_C5 <-  Mm_msigdb %>%
  dplyr::filter(gs_cat == "C5" & gs_subcat == "GO:BP") %>%
  dplyr::select(gs_name,gs_exact_source,entrez_gene) 
Mm_C5_GO <- Mm_C5[,1:2] %>%
  unique()
saveRDS(Mm_C5_GO,file = "res/R/Mm_C5_GO.rds")


Mm_C5_GO <- readRDS(file = "res/R/Mm_C5_GO.rds")
eLR.df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds") %>%
  rownames_to_column("LRpairs") %>%
  mutate(Lgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[2]})))

for(ii in 1:6){
  
  tmp.Lgene <- eLR.df %>%
    dplyr::filter(Cluster == ii) %>%
    pull(Lgene)
  tmp.Rgene <- eLR.df %>%
    dplyr::filter(Cluster == ii) %>%
    pull(Rgene)
  
  #### ligand
  #test.gene <- unique(c(tmp.Lgene,tmp.Rgene))
  test.gene <- unique(tmp.Lgene)
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
    arrange(p.adjust) %>%
    mutate(Description=gsub(pattern = "GOBP_",replacement = "",Description)) %>%
    mutate(Description=str_replace_all(Description,pattern = "_",replacement = " ")) %>%
    mutate(Description=str_to_sentence(Description)) 
  
  #em.res$Description
  write.table(em.res,
              file = myFileName(prefix = paste0("res/txt/",ii,"_eLR_","ligand","_msigdb_GO"),suffix = ".txt"),
              quote = F,
              sep = "\t",
              row.names = F)
  
  p <- myenrichr_GOplot(df = em.res,
                        fill.color = "#C6DBEFFF",
                        show_number = 20,
                        title = paste0("C",ii," ",length(tmp.Lgene)," eLR ","Ligand ",length(gene.list)," gene"),
                        font.size = 18,
                        p.adjust.cut = 0.05,
                        plot.ylab = NULL,
                        term.pos.adjust = 0)+
    scale_y_continuous(expand = c(0,0))+
    theme(axis.line.y = element_blank(),
          axis.line.x = element_line(size = 1),
          axis.ticks = element_line(size = 1))
  p
  ggsave(filename = myFileName(prefix = paste0("res/fig/","Ligand_Cluster_",ii,"_msigdb_GO"),suffix = ".png"),
         width = 8,height = 8,dpi = 350)
  
  #### receptor
  test.gene <- unique(tmp.Rgene)
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
    arrange(p.adjust) %>%
    mutate(Description=gsub(pattern = "GOBP_",replacement = "",Description)) %>%
    mutate(Description=str_replace_all(Description,pattern = "_",replacement = " ")) %>%
    mutate(Description=str_to_sentence(Description)) 
  
  #em.res$Description
  write.table(em.res,
              file = myFileName(prefix = paste0("res/txt/",ii,"_eLR_","receptor","_msigdb_GO"),suffix = ".txt"),
              quote = F,
              sep = "\t",
              row.names = F)
  
  p <- myenrichr_GOplot(df = em.res,
                        fill.color = "#C6DBEFFF",
                        show_number = 20,
                        title = paste0(ii," ",length(tmp.Lgene)," eLR ","Receptor ",length(gene.list)," gene"),
                        font.size = 18,
                        p.adjust.cut = 0.05,
                        plot.ylab = NULL,
                        term.pos.adjust = 0)+
    scale_y_continuous(expand = c(0,0))+
    theme(axis.line.y = element_blank(),
          axis.line.x = element_line(size = 1),
          axis.ticks = element_line(size = 1))
  p
  ggsave(filename = myFileName(prefix = paste0("res/fig/","Receptor_Cluster_",ii,"_msigdb_GO"),suffix = ".png"),
         width = 8,height = 8,dpi = 350)
  
  
  #### ligand and Receptor
  #test.gene <- unique(c(tmp.Lgene,tmp.Rgene))
  test.gene <- unique(c(tmp.Lgene,tmp.Rgene))
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
    arrange(p.adjust) %>%
    mutate(Description=gsub(pattern = "GOBP_",replacement = "",Description)) %>%
    mutate(Description=str_replace_all(Description,pattern = "_",replacement = " ")) %>%
    mutate(Description=str_to_sentence(Description)) 
  
  #em.res$Description
  write.table(em.res,
              file = myFileName(prefix = paste0("res/txt/",ii,"_eLR_","LR","_msigdb_GO"),suffix = ".txt"),
              quote = F,
              sep = "\t",
              row.names = F)
  
  p <- myenrichr_GOplot(df = em.res,
                        fill.color = "#C6DBEFFF",
                        show_number = 20,
                        title = paste0(ii," ",length(tmp.Lgene)," eLR ","LR ",length(gene.list)," gene"),
                        font.size = 18,
                        p.adjust.cut = 0.05,
                        plot.ylab = NULL,
                        term.pos.adjust = 0)+
    scale_y_continuous(expand = c(0,0))+
    theme(axis.line.y = element_blank(),
          axis.line.x = element_line(size = 1),
          axis.ticks = element_line(size = 1))
  p
  ggsave(filename = myFileName(prefix = paste0("res/fig/","LR_Cluster_",ii,"_msigdb_GO"),suffix = ".png"),
         width = 8,height = 8,dpi = 350)
}


#####-----2.5 eLR KEGG ---------

tmp.gene.list <- lapply(1:6,FUN = function(ii){
  eLR.df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds") %>%
    rownames_to_column("LRpairs") %>%
    mutate(Lgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[1] }))) %>%
    mutate(Rgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[2] })))
  tmp.Lgene <- eLR.df %>%
    dplyr::filter(Cluster == ii) %>%
    pull(Lgene)
  tmp.Rgene <- eLR.df %>%
    dplyr::filter(Cluster == ii) %>%
    pull(Rgene)
  tmp.gene <- unique(c(tmp.Lgene,tmp.Rgene))
  eg = bitr(tmp.gene, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)
  tmp.gene <- eg$ENTREZID
  return(tmp.gene)
})

names(tmp.gene.list) <- paste0("cluster_",1:6)

#tmp.df <- tmp.em@compareClusterResult

# Compare <- compareCluster(geneCluster=tmp.gene.list,
#                           fun="enrichKEGG",
#                           organism = "mmu",
#                           pvalueCutoff=0.05)

Compare <- compareCluster(geneCluster=tmp.gene.list,
                          fun="enrichKEGG",
                          organism = "mmu",
                          pvalueCutoff=0.05,use_internal_data=T)

dotplot(Compare,showCategory=20,font.size = 20)+
  scale_x_discrete(labels=paste0("cluster_",1:6))+
  scale_color_gradientn(colors  = coolwarm(100))+
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5))

# ggsave(filename = "res/fig/fig2_show_eLR_KEGG_important_20230329.pdf",width = 12,height = 12,dpi = 350,useDingbats=F)
# ggsave(filename = "res/fig/fig2_show_eLR_KEGG_important_20230329.png",width = 12,height = 12,dpi = 350)


# tmp.compare <- compareCluster(geneCluster=tmp.gene.list,
#                                          fun = "enrichGO",
#                                          OrgDb = org.Mm.eg.db,
#                                          readable = T,
#                                          ont = "BP",
#                                          pvalueCutoff=0.05)
# 
# dotplot(tmp.compare,showCategory=20,font.size = 20)


tmp.term <- readLines(con  = "res/txt/eLR_KEGG_term_20220531")
Compare <- setReadable(Compare,OrgDb = org.Mm.eg.db,keyType = "ENTREZID")
tmp.df <- Compare@compareClusterResult 
tmp.df$Description <- gsub(pattern = " - Mus musculus \\(.*?\\)",replacement = "",tmp.df$Description,perl = T)
head(tmp.df$Description)
####sapply is good!
tmp.df <- tmp.df %>%
  dplyr::filter(Description %in% tmp.term)  %>%
  mutate(GeneRatio = sapply(GeneRatio,
                            FUN = function(xx) eval(parse(text = xx)),
                            USE.NAMES = F)) %>%
  mutate(Description = factor(Description,levels = rev(tmp.term)))
####head(tmp.df)



p <- ggplot(data = tmp.df,aes(Cluster, Description,size = GeneRatio,color=pvalue))+
  geom_point() +
  scale_x_discrete(labels = seq(1:6))+
  scale_color_gradientn(colors = rev(coolwarm(100)),limits = c(NA,0.05))+
  myGridTheme(font.size = 22)+
  ylab(NULL)

p

myggsave(p,prefix = "res/fig/fig2_show_eLR_KEGG_important",suffix = ".png",
         width = 16,height = 8,dpi = 350)

eLR_row_annotation_df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds")
head(eLR_row_annotation_df)
eLR_row_annotation_df <- eLR_row_annotation_df %>%
  rownames_to_column("eLR") %>%
  mutate(Cluster = paste0("cluster_",Cluster))
head(eLR_row_annotation_df)
write.table(eLR_row_annotation_df,file = "res/txt/eLR_classification_20220907.txt",
            quote = F,sep = "\t",row.names = F)

saveRDS(tmp.df,file = "res/R/eLR_cluster_KEGG_enriched_results_20230329.rds")

###-------2.6 Hippo Signaling pathway--------
### start 20230407

###-------2.6.1 check Hippo signaling pathway related eLR---------

tmp.term <- "Hippo signaling pathway"
tmp.df <- readRDS(file = "res/R/eLR_cluster_KEGG_enriched_results_20230329.rds") %>%
  dplyr::filter(Description %in% tmp.term) 
eLR.df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds") %>%
  rownames_to_column("LRpairs") %>%
  mutate(Lgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[1] }))) %>%
  mutate(Rgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[2] }))) %>%
  mutate(Cluster = paste0("cluster_",Cluster))

ii <- 1 
tmp.df.res <- lapply(1:4, function(ii){
  tmp.gene <- unlist(strsplit(tmp.df$geneID[ii],split = "/"))
  tmp.cluster <- as.character(tmp.df$Cluster)[ii]
  tmp.df.use <- eLR.df %>%
    dplyr::filter(Cluster %in% tmp.cluster) %>%
    dplyr::filter((Lgene %in% tmp.gene) | Rgene %in% tmp.gene) %>%
    dplyr::select(LRpairs,Lgene,Rgene,Cluster)
  return(tmp.df.use)
})
tmp.df.res <- Reduce(f = rbind,tmp.df.res)

write.table(tmp.df.res,file = "res/txt/Hippo_signaling_related_eLR_20230407.txt",sep = "\t",quote = F,row.names = F)



# Cluster       ID             Description  GeneRatio  BgRatio       pvalue     p.adjust       qvalue
# 1 cluster_1 mmu04390 Hippo signaling pathway 0.09821429 157/9070 3.433265e-06 2.980698e-05 1.856263e-05
# 2 cluster_2 mmu04390 Hippo signaling pathway 0.17777778 157/9070 8.439717e-07 2.059291e-05 1.190444e-05
# 3 cluster_4 mmu04390 Hippo signaling pathway 0.11111111 157/9070 9.617365e-10 1.029624e-08 6.431427e-09
# 4 cluster_6 mmu04390 Hippo signaling pathway 0.17647059 157/9070 2.201878e-05 1.456627e-04 9.627646e-05
# geneID Count
# 1                                 Bmp6/Bmp7/Tgfb3/Wnt5a/Tgfb2/Wnt7a/Cdh1/Bmpr1b/Bmpr2/Fzd3/Fzd7    11
# 2                                                Tgfb3/Bmp5/Bmp6/Cdh1/Tgfb2/Tgfbr1/Bmpr2/Bmpr1b     8
# 3 Cdh1/Wnt3a/Areg/Gdf5/Wnt5a/Wnt7a/Fgf1/Tgfbr1/Fzd2/Itgb2/Fzd5/Fzd1/Fzd4/Fzd9/Bmpr1b/Bmpr2/Fzd3    17
# 4                                                               Wnt3a/Wnt7b/Fzd7/Fzd5/Fzd3/Fzd6     6           


tmp.LR <- eLR.df




###------2.6.2 check the core Hippo signaling pathway gene --------

####------2.6.2.1 load 19 core Hippo signaling gene--------------------------
tmp.gene <- c("STK3","TAZ","LATS1","TAOK2","TAOK1",
              "LATS2","WWC1","TEAD4","NF2","YAP1",
              "STK4","TAOK3","FRMD6","TEAD2","TEAD3",
              "SAV1","TEAD1","MOB1B","MOB1A")

bitr(tmp.gene,fromType = "")
cat(tmp.gene,sep = "\n")
####------2.6.2.2 load data -------------------
### visualize gene expression by pseudotime 
###load data
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
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
  dplyr::filter(Stage %in% Stage.level) %>%
  pull(cell_id)
beforeEPI.exp <- beforeEPI.exp[,tmp.cell_id]

####load cds data
beforeEPI.cds <- readRDS(file = "res/R/early.scRNAseq.moncole_2022031116.rds")
beforeEPI.meta <- pData(beforeEPI.cds)

tmp.cell.id <- beforeEPI.meta %>%
  group_by(Stage) %>%
  mutate(tmp.rank = row_number(Pseudotime)) %>%
  arrange(Stage,tmp.rank) %>%
  pull(cell_id)

data.plot <- beforeEPI.meta[tmp.cell.id,] %>%
  mutate(tmp.rank = row_number(),
         tmp.rect = 1) %>%
  dplyr::select(Stage,tmp.rank,tmp.rect)

tmp.levels.2 <- c("MII oocyte","zygote",
                  "early 2-cell","mid 2-cell","late 2-cell",
                  "4-cell","8-cell","16-cell",
                  "early blastocyst","mid blastocyst","late blastocyst")
tmp.Stage.color <- col.spectral(length(tmp.levels.2))
names(tmp.Stage.color) <- tmp.levels.2

p_bottom <- ggplot(data.plot,aes(x = tmp.rank,
                                 y = tmp.rect,fill=Stage))+
  geom_tile()+
  scale_fill_manual(values = tmp.Stage.color)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab(NULL)+
  xlab("rank by pseudotime")+
  theme_cowplot(font_size = 28)+
  NoLegend()+
  theme(axis.ticks = element_line(size = 2),
        axis.text.y = element_blank(),
        #axis.title.x = element_text(size = 22),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_blank())
p_bottom
tmp.Lgene <- "Fgf4"
tmp.Rgene <- "Fgfr2"
tmp.TF <- "Gata6"
pal_aaas()(4)
tmp.color <- c("#ffc000","#00b050","#EE0000FF","#631879FF")
names(tmp.color) <- c("Ligand","Receptor","TF","IS")

tmp.TF.IS.plot <- function(tmp.Lgene,tmp.Rgene,tmp.TF){
  
  tmp.data.plot.1 <- mat.df %>%
    dplyr::filter(gene_symbol == tmp.Lgene) %>%
    left_join(beforeEPI.meta,by = "cell_id") %>%
    dplyr::select(cell_id,gene_symbol,gene.exp,Pseudotime,cell_type,Stage) %>%
    mutate(tmp.rank = row_number(Pseudotime)) %>%
    arrange(Stage,tmp.rank) %>%
    mutate(tmp.rank = row_number()) %>%
    mutate(group = "Ligand")
  
  tmp.data.plot.2 <- mat.df %>%
    dplyr::filter(gene_symbol == tmp.Rgene) %>%
    left_join(beforeEPI.meta,by = "cell_id") %>%
    dplyr::select(cell_id,gene_symbol,gene.exp,Pseudotime,cell_type,Stage) %>%
    mutate(tmp.rank = row_number(Pseudotime)) %>%
    arrange(Stage,tmp.rank) %>%
    mutate(tmp.rank = row_number()) %>%
    mutate(group = "Receptor")
  
  tmp.data.plot.3 <- mat.df %>%
    dplyr::filter(gene_symbol == tmp.TF) %>%
    left_join(beforeEPI.meta,by = "cell_id") %>%
    dplyr::select(cell_id,gene_symbol,gene.exp,Pseudotime,cell_type,Stage) %>%
    mutate(tmp.rank = row_number(Pseudotime)) %>%
    arrange(Stage,tmp.rank) %>%
    mutate(tmp.rank = row_number()) %>%
    mutate(group = "TF")
  
  
  tmp.data.plot.4 <- tmp.data.plot.2 %>%
    mutate(gene.exp = sqrt(tmp.data.plot.1$gene.exp * tmp.data.plot.2$gene.exp)) %>%
    mutate(group = "IS")
  
  tmp.data.plot <- rbind(tmp.data.plot.1,tmp.data.plot.2,
                         tmp.data.plot.3,tmp.data.plot.4)
  
  tmp.data.plot$group <- factor(tmp.data.plot$group,
                                levels = c("Ligand","Receptor","TF","IS"))
  
  
  tmp.x <- tmp.data.plot.3$gene.exp
  tmp.y <- tmp.data.plot.4$gene.exp
  
  tmp.pvalue.TFtoLR <- lmtest::grangertest(tmp.x,tmp.y)
  tmp.pvalue.TFtoLR <- tmp.pvalue.TFtoLR$`Pr(>F)`[2]
  tmp.pvalue.LRtoTF <- lmtest::grangertest(tmp.y,tmp.x)
  tmp.pvalue.LRtoTF <- tmp.pvalue.LRtoTF$`Pr(>F)`[2]
  
  p <- ggplot(tmp.data.plot,aes(tmp.rank,gene.exp,group=group,col=group))+
    geom_point(alpha=0.3)+
    geom_smooth(se = F,method = "gam", formula = y ~ s(x, bs = "cs"),alpha=0.3)+
    xlab("rank by pseudotime")+
    ylab("log2(RPKM+1)")+
    ggtitle(label = paste0(tmp.Lgene,"-",tmp.Rgene,",",tmp.TF,"\n",
                           "p_TFtoLR=",signif(tmp.pvalue.TFtoLR,digits = 4),"\n",
                           "p_LRtoTF=",signif(tmp.pvalue.LRtoTF,digits = 4)))+
    scale_color_manual(values = tmp.color,name=NULL)+
    scale_x_continuous(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme_cowplot(font_size = 28)+
    theme(legend.position = "top",legend.justification = "center",
          axis.line = element_line(size = 2),
          axis.ticks = element_line(size = 2),
          plot.title = element_text(hjust = 0.5))+
    xlab(NULL)+
    theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank())
  
  p_res <- wrap_plots(p,p_bottom,ncol = 1,heights = c(9,1)) &
    theme(plot.margin = margin(t = 0,r = 0,b = 0,l = 0))
  
  return(p_res)
  
}


write.table(tmp.df.res,file = "res/txt/Hippo_signaling_related_eLR_20230407.txt",sep = "\t",quote = F,row.names = F)

Hippo_eLR <- read.delim(file = "res/txt/Hippo_signaling_related_eLR_20230407.txt")

lapply(1:nrow(Hippo_eLR),function(ii){
  tmp.Lgene <- Hippo_eLR$Lgene[ii]
  tmp.Rgene <- Hippo_eLR$Rgene[ii]
  tmp.eLR <- Hippo_eLR$LRpairs[ii]
  tmp.TF <- "Tead2"
  p <- tmp.TF.IS.plot(tmp.Lgene = tmp.Lgene,tmp.Rgene = tmp.Rgene,tmp.TF = tmp.TF)
  ggsave(filename = myFileName(prefix = paste0("res/fig/Hippo_eLR_dynamics_",tmp.eLR,"_",tmp.TF), suffix = ".png"),
         width = 8,height = 8)
  
})

lapply(1:nrow(Hippo_eLR),function(ii){
  tmp.Lgene <- Hippo_eLR$Lgene[ii]
  tmp.Rgene <- Hippo_eLR$Rgene[ii]
  tmp.eLR <- Hippo_eLR$LRpairs[ii]
  tmp.TF <- "Tead4"
  p <- tmp.TF.IS.plot(tmp.Lgene = tmp.Lgene,tmp.Rgene = tmp.Rgene,tmp.TF = tmp.TF)
  ggsave(filename = myFileName(prefix = paste0("res/fig/Hippo_eLR_dynamics_",tmp.eLR,"_",tmp.TF), suffix = ".png"),
         width = 8,height = 8)
  
})






####--------3. investigate network result -----------

####--------3.1 load data----------- 
eLR_row_annotation_df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds")
eLR_embryo_eLR_TF <- readRDS(file = "res/R/early_embryo_eLR_TF_interaction_information_2022060620.rds")
tmp.TF.annotation <- readRDS(file = "res/R/early_embryo_core_TF_annotation_row_2022032521.rds")
tmp.df.use <- tmp.TF.annotation %>%
  rownames_to_column("TF")
nrow(tmp.df.use)
rtna <- readRDS(file = "res/R/early.rtna_20220316.rds")
regulon_list <- RTN::tna.get(rtna,what = "regulons")

tmp.gene.use <- unique(c(eLR_embryo_eLR_TF$Lgene,eLR_embryo_eLR_TF$Rgene))
tmp.TF <- tmp.df.use$TF
tmp.res.length <- unlist(lapply(tmp.TF, function(ii){
  x <- regulon_list[[ii]]
  tmp.res <- length(x[which(x %in% tmp.gene.use )])/length(x)
  return(tmp.res)
}))
names(tmp.res.length) <- tmp.TF
#barplot(tmp.res.length)

tmp.data.plot <- data.frame(ratio = tmp.res.length,stringsAsFactors = F)
tmp.data.plot <- tmp.data.plot %>%
  rownames_to_column("TF")
tmp.data.plot <- merge(tmp.df.use,tmp.data.plot,by="TF")

tmp.data.plot$TF <- factor(tmp.data.plot$TF,
                           levels = rev(rownames(tmp.TF.annotation)))

tmp.data.plot$rect <- 1

tmp.data.plot <- tmp.data.plot %>%
  arrange(TF) %>%
  mutate(tmp.rank = row_number()) %>%
  mutate(group = ifelse(ratio > 0.05,"high","low"))

tmp.data.plot$Cluster <- paste0("C",tmp.data.plot$Cluster)

#####--------3.2 plot ratio------------

p1 <- ggplot(tmp.data.plot,aes(rect,tmp.rank,fill=Cluster))+
  geom_tile()+
  theme_cowplot(font_size = 28)+
  scale_fill_manual(values = monet.sunumbrella.women(6),
                    name = "tTFs",
                    labels = paste0("C",1:6))+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  NoAxes()+
  theme(axis.text.y = element_blank())

p2 <- ggplot(tmp.data.plot,aes(ratio,TF))+
  geom_bar(stat = "identity",fill = "grey")+
  ylab(NULL)+
  xlab("#(eLR gene contained in regulon)/(regulon size)")+
  theme_cowplot(font_size = 28)+
  scale_x_continuous(expand = c(0,0))+
  theme(axis.text.y  = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.line = element_line(size = 1),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_line(size = 1))

p <- wrap_plots(p1,p2,guides = "collect",widths = c(1,19)) & 
  theme(plot.margin = margin(r = 0,  l = 0))
p
ggsave(p,filename = myFileName("res/fig/fig4_regulon_eLR_ratio",
                               suffix = ".png"),
       width = 8,height = 6,dpi = 400)

####---------3.3 plot percent------------

test.plot <- tmp.data.plot %>% 
  group_by(Cluster,group) %>%
  summarise(n = n()) %>%
  mutate(ratio = n/sum(n)) %>%
  ungroup() %>%
  mutate(group = factor(group,levels = c("low","high")))

ggplot(test.plot,aes(Cluster,ratio,fill=group))+
  geom_bar(stat = "identity")+
  theme_cowplot(font_size = 28,rel_large = 9/14,
                rel_small = 10/14)+
  xlab("tTFs")+
  ylab("percent")+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = c("navy","brown3"),
                    name="level",labels=c("low(ratio < 0.05)",
                                           "high(ratio) >= 0.05"))+
  theme(legend.position = "top",
        legend.justification = "center",
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))
ggsave(filename = myFileName(prefix = "res/fig/fig4_regulon_eLR_ratio_statistic",
                             suffix = ".png"),width = 8,height = 6,
       dpi = 400)


####-------3.4 further investigation------------


tmp.list <- RTN::tna.get(rtna,what = "regulons")
tmp.gene.use <- unique(c(eLR_embryo_eLR_TF$Lgene,eLR_embryo_eLR_TF$Rgene))
x <- tmp.list[["Nfya"]]
x[which(x %in% tmp.gene.use )]

tmp.list[["Cdx2"]]

x <- tmp.list[["Tead2"]]
x[which(x %in% tmp.gene.use )]

x <- tmp.list[["Gata3"]]
x
x[which(x %in% tmp.gene.use )]
x <- tmp.list[["Gata4"]]
tmp.gene.use <- unique(c(eLR_embryo_eLR_TF$Lgene,eLR_embryo_eLR_TF$Rgene))
x[which(x %in% tmp.gene.use )]

x <- tmp.list[["Zscan4c"]]
tmp.gene.use <- unique(c(eLR_embryo_eLR_TF$Lgene,eLR_embryo_eLR_TF$Rgene))
x[which(x %in% tmp.gene.use )]

x <- tmp.list[["Sox21"]]
tmp.gene.use <- unique(c(eLR_embryo_eLR_TF$Lgene,eLR_embryo_eLR_TF$Rgene))
x[which(x %in% tmp.gene.use )]

x <- tmp.list[["Obox6"]]
tmp.gene.use <- unique(c(eLR_embryo_eLR_TF$Lgene,eLR_embryo_eLR_TF$Rgene))
x[which(x %in% tmp.gene.use )]

x <- tmp.list[["Nanog"]]
x[which(x %in% tmp.gene.use )]

x <- tmp.list[["Pou5f1"]]
x[which(x %in% tmp.gene.use )]

tTF <- tmp.TF.meta %>%
  pull(gene)

tmp.res.list <- lapply(tTF,FUN = function(ii){
  cat(ii,sep = "\n")
  x <- tmp.list[[ii]]
  tmp.idx <- which(x %in% "Fgfr1")
  if(!purrr::is_empty(tmp.idx)){
    return(ii)
  }
})
unlist(tmp.res.list)


tmp.res.list <- lapply(tTF,FUN = function(ii){
  cat(ii,sep = "\n")
  x <- tmp.list[[ii]]
  return(x[which(x %in% tmp.gene.use )])
})
unlist(tmp.res.list)


trust_TF_network <- read.delim(file = "database/trrust_rawdata.mouse.tsv",
                               stringsAsFactors = F,header = F)



####-------3.5 eLR to TF statistic---------------
####load LR data
ligand_target_matrix = readRDS("res/R/nichenet_ligand_target_matrix.rds")
lr_network <- readRDS(file = "res/R/nichenet_lr_network_mouse.rds")
mouse_sig_network <- readRDS(file = "res/R/nichenet_mouse_sig_network.rds")
weighted_networks <- readRDS(file = "res/R/nichenet_mouse_weighted_network_20210131.rds")
ligand_tf_matrix = readRDS(file = "res/R/nichenet_ligand_tf_matrix.rds")
### head(ligand_tf_matrix)

####load TF data and eLR data
eLR_row_annotation_df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds")
eLR_embryo_eLR_TF <- readRDS(file = "res/R/early_embryo_eLR_TF_interaction_information_2022060620.rds")
tmp.TF.annotation <- readRDS(file = "res/R/early_embryo_core_TF_annotation_row_2022032521.rds")


num_clusters <- 6
tmp.color <- ggsci::pal_npg()(6)


for(jj in 1:6){
  
  cat(paste0("eLR_Cluster_",jj),sep = "\n")
  
  tmp_eLR_df_use <- eLR_row_annotation_df %>%
    dplyr::filter(Cluster == jj)
  tmp.ratio.list <- lapply(1:num_clusters,FUN = function(ii){
    cat(ii,sep = "\n")
    tmp.TF.use <- tmp.TF.annotation %>%
      rownames_to_column("TF") %>%
      filter(Cluster == ii) %>%
      pull(TF)
    tmp.var <- paste0("TF_cluster_",ii)
    
    tmp.eLR.use <- tmp_eLR_df_use %>%
      rownames_to_column("LRpairs") %>%
      dplyr::filter(!!sym(tmp.var) %in% c("forward","feedback")) %>%
      pull(LRpairs)
    tmp.res.hits <- lapply(tmp.eLR.use,FUN = function(xxx){
      tmp.eLR.L <- unlist(lapply(str_split(string = xxx, pattern = "-"), function(x) x[1]))
      tmp.eLR.R <- unlist(lapply(str_split(string = xxx, pattern = "-"), function(x) x[2]))
      tmp.eLR.LR.gene <- unique(c(tmp.eLR.L,tmp.eLR.R))
      tmp.df <- weighted_networks$lr_sig %>%
        dplyr::filter(to %in% tmp.TF.use,
                      from %in% tmp.eLR.LR.gene)
      if(nrow(tmp.df) >= 1){
        tmp.hits <- 1
      }else{
        tmp.hits <- 0
      }
      return(tmp.hits)
    })
    tmp.ratio.nichenet <- sum(unlist(tmp.res.hits)) / nrow(tmp_eLR_df_use)
    tmp.ratio.inference <- length(tmp.eLR.use) / nrow(tmp_eLR_df_use)
    
    tmp.res.df <- data.frame(infer=tmp.ratio.inference,
                             nichenet=tmp.ratio.nichenet,
                             TF_cluster=paste0("C",ii))
    return(tmp.res.df)
  })
  
  tmp.ratio.df <- Reduce(rbind,tmp.ratio.list)
  tmp.data.plot <- tmp.ratio.df %>%
    gather(key = "group",value = "ratio",-TF_cluster)
  
  tmp.color.use <- tmp.color[jj]
  p <- ggplot(tmp.data.plot,aes(TF_cluster,ratio,
                                fill=group))+
    geom_bar(stat = "identity",
             position = position_dodge())+
    xlab("tTFs")+
    ggtitle(paste0("Cluster",jj," eLR to TF ratio"))+
    scale_fill_manual(values = c(tmp.color.use,
                                 alpha(tmp.color.use,
                                       alpha = 0.6)),
                      name=NULL)+
    #scale_x_discrete(expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme_cowplot(font_size = 28)+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "top",
          legend.justification = "center",
          axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 1))
  
  myggsave(p = p,
           prefix = paste0("res/fig/fig4_Cluster_",
                           jj,"_eLR_TF_ratio"),
           suffix = ".png",
           dpi=400,
           width = 8,bg="white",
           height = 8)
}



# ####--------4. analysis motif enrichment result--------
# 
# ####--------4.1 write motif to meme format--------------
# 
# tmp.TF.list <- list.files("database/CISBP/")
# lapply(tmp.TF.list , function(ii){
#   tmp.TF <- ii
#   tmp.path <- list.files(paste0("database/CISBP/",tmp.TF),full.names = T)
#   tmp.data <- read.delim(file = tmp.path,stringsAsFactors = F,header = T)
#   head(tmp.data)
#   tmp.data <- t(tmp.data[,-1])
#   tmp.motif <- create_motif(input = as.matrix(tmp.data),
#                             type = "PPM",
#                             alphabet = "DNA",
#                             name = tmp.TF)
#   ###view_motifs(tmp.motif,relative_entropy = T)
#   tmp.motif@alphabet <- "DNA"
#   dir.create(paste0("res/meme/TimeTalk/",tmp.TF))
#   write_meme(tmp.motif,file = paste0("res/meme/TimeTalk/",tmp.TF,"/",tmp.TF,".meme"))
# })
# 
# ####--------4.2 plot tracks----------------
# 
# tmp.chrom <- read.delim(file = "database/mm9.chrom.sizes",header = F,stringsAsFactors = F)
# head(tmp.chrom)
# colnames(tmp.chrom) <- c("chrom","length")
# tmp.metadata <- data.frame(name="Genome",value="mm9",stringsAsFactors = F)
# ####load txdb
# mm9KG_txdb <- makeTxDbFromGFF(file = "D://Ubuntu/wlt/igenomes/Mus_musculus/UCSC/mm9/Annotation/genes.gtf",
#                               chrominfo = tmp.chrom,
#                               metadata = tmp.metadata)
# 
# #### the following function adapted from signac
# txbd <- mm9KG_txdb
# tx <- myGetGRangesFromsTxDb(txdb = mm9KG_txdb,standard.chromosomes = T,verbose = T)
# saveRDS(tx,file = "res/R/mm9_tx_Granges_20220617.rds")
# 
# tx <- readRDS(file = "res/R/mm9_tx_Granges_20220617.rds")
# tmp.region <- "chr7:152044791-152053648"
# tmp.ylim <- c(0,1)
# 
# tmp.track.plot.function <- function(tx,
#                                     tmp.region,
#                                     tmp.ylim = c(0,1)){
#   
#   p_gene <- myGenePlot(annotation = tx,
#                        region = tmp.region,
#                        arrow_sbreaks = 600,
#                        font_size = 18,label_size = 8)+
#     xlab(NULL)+
#     theme(plot.margin = unit(c(0,0,0,0), "cm"),
#           axis.title.y = element_text(angle = 0,vjust = 0.5,hjust = 0.5),
#           axis.line = element_blank(),axis.text.x = element_blank(),
#           axis.ticks = element_blank())
# 
#   tmp.dir <- "res/atac-seq-process/data/bigwig/"
#   tmp.files <- list.files(path = tmp.dir,pattern = "bw",full.names = T)
#   tmp.levels <- c("2cell_early","2cell","4cell","8cell","ICM")
#   tmp.files <- paste0("res/atac-seq-process/data/bigwig/GSE66581_atac_",tmp.levels,".bw")
#   
#   ###tmp.colors <- shaman.color(length(tmp.files))
#   tmp.colors <- c("#4d311a","#d29156","#7f8a65","#6d82a3","#6e648a")
#   #scales::show_col(tmp.colors)
#   
#   p_ATAC_track_list <- lapply(seq_along(tmp.files),function(ii){
#     cat(ii,sep = "\n")
#     p <- myBigwigTrack(region = as(tmp.region,"GRanges"),
#                        bigwig = tmp.files[ii],
#                        smooth = 100,
#                        lognorm = T,
#                        type = "coverage",
#                        y_label = tmp.levels[ii],
#                        fontsize=18,
#                        track.color=tmp.colors[ii],
#                        tmp.ylimits=tmp.ylim,
#                        max.downsample = 3000,
#                        downsample.rate = 0.3,
#                        tmp.seed=42)
#     
#     return(p)
#   })
#   
#   p_ATAC <- patchwork::wrap_plots(p_ATAC_track_list,ncol = 1) +
#     plot_annotation(title = tmp.region) & 
#     theme(plot.title = element_text(hjust = 0.5,size = 18),
#           plot.margin = unit(c(0,0,0,0), "cm"),
#           axis.title.y = element_text(size = 18,angle = 0,vjust = 0.5,hjust = -1),
#           axis.ticks.x = element_blank(),
#           axis.text.x  = element_blank(),
#           axis.title.x = element_blank()) 
#   
#   
#   tmp.levels <- c("2cell_early","2cell","4cell","8cell","ICM")
#   tmp.TF.name <- "Gata3"
#   
#   tmp.files.list <- paste0("res/meme/TimeTalk/",tmp.TF.name,"/TFbindsites/",
#                            tmp.levels,
#                            "_peaks_",tmp.TF.name,"_bindsite.bed")
#   p_bed_1 <- lapply(seq_along(tmp.files.list), function(ii){
#     cat(ii,sep = "\n")
#     tmp.files <- tmp.files.list[ii]
#     p <- myBedTrack(bed.input = tmp.files,
#                     tmp.region = tmp.region,
#                     track.color = tmp.colors[ii],
#                     tmp.label = paste0(tmp.TF.name,"_",tmp.levels[ii]),
#                     font_size = 18)
#     return(p)
#   })
#   
#   tmp.TF.name <- "Pou5f1"
#   tmp.files.list <- paste0("res/meme/TimeTalk/",tmp.TF.name,"/TFbindsites/",
#                            tmp.levels,
#                            "_peaks_",tmp.TF.name,"_bindsite.bed")
#   
#   p_bed_2 <- lapply(seq_along(tmp.files.list), function(ii){
#     cat(ii,sep = "\n")
#     tmp.files <- tmp.files.list[ii]
#     p <- myBedTrack(bed.input = tmp.files,
#                     tmp.region = tmp.region,
#                     track.color = tmp.colors[ii],
#                     tmp.label = paste0(tmp.TF.name,"_",tmp.levels[ii]),
#                     font_size = 18)
#     return(p)
#   })
#   
#   
#   p_list <- c(p_ATAC_track_list,p_bed_1,p_bed_2)
#   
#   p <- patchwork::wrap_plots(p_list,ncol = 1) +
#     plot_annotation(title = tmp.region) & 
#     theme(plot.title = element_text(hjust = 0.5,size = 18),
#           plot.margin = unit(c(0,0,0,0), "cm"),
#           axis.title.y = element_text(size = 18,angle = 0,vjust = 0.5,hjust = -1),
#           axis.ticks.x = element_blank(),
#           axis.text.x  = element_blank(),
#           axis.title.x = element_blank()) 
#   p / p_gene
#   
# }
# tmp.region <- "chr7:152044791-152053648"
# p1 <- tmp.track.plot.function(tx = tx,tmp.region = tmp.region,tmp.ylim = c(0,1))
# p1
# tmp.region <- "chr7:137303465-137412822"
# p2 <- tmp.track.plot.function(tx = tx,tmp.region = tmp.region,tmp.ylim = c(0,1))
# 
# p2
# 
# p <- p1 | (p2 & theme(axis.text.y = element_blank(),
#                       axis.title.y = element_blank()))
# 
# myggsave(p,prefix = "res/fig/fig4_track_illustrate",
#          suffix = ".png",width = 24,height = 12,dpi = 400)
# 
# myggsave(p,prefix = "res/fig/fig4_track_illustrate",
#          suffix = ".pdf",width = 24,height = 12,dpi = 400)
# 
# 
# 
# ####-----------4.3 read motifs -----------
# 
# ### mm9 gene
# mm9KG_txdb <- makeTxDbFromGFF(file = "D://Ubuntu/wlt/igenomes/Mus_musculus/UCSC/mm9/Annotation/genes.gtf")
# mm9.gene <- genes(mm9KG_txdb)
# mm9.gene.list <- mm9.gene$gene_id
# mm9.tss <- promoters(mm9.gene,upstream = 2500, downstream = 2500)
# mm9.tss
# 
# 
# ### load NicheNet work
# ligand_target_matrix = readRDS("res/R/nichenet_ligand_target_matrix.rds")
# lr_network <- readRDS(file = "res/R/nichenet_lr_network_mouse.rds")
# mouse_sig_network <- readRDS(file = "res/R/nichenet_mouse_sig_network.rds")
# weighted_networks <- readRDS(file = "res/R/nichenet_mouse_weighted_network_20210131.rds")
# ligand_tf_matrix = readRDS(file = "res/R/nichenet_ligand_tf_matrix.rds")
# 
# tmp.gr.network <- weighted_networks$gr
# 
# 
# 
# ### LR gene
# eLR_row_annotation_df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds")
# LRpairs <- rownames(eLR_row_annotation_df)
# Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "-"),FUN = function(x) x[1])) 
# Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "-"),FUN = function(x) x[2])) 
# Lgene <- unique(Lgenelist)
# Rgene <- unique(Rgenelist)
# LRgene <- c(Lgene,Rgene)
# LRgene <- intersect(LRgene,mm9.gene.list)
# 
# 
# tmp.df <- tmp.gr.network %>%
#   dplyr::filter(to %in% c("Fgf4","Fgfr2","Fgfr1"),
#                 from %in% rownames(tmp.TF.annotation))
# 
# ### TFBS
# tmp.TF.name <- "NFYA"
# tmp.dir <- paste0("res/meme/TimeTalk/",tmp.TF.name,"/TFbindsites/")
# #list.files(path = tmp.dir)
# tmp.samples <- c("2cell_early","2cell","4cell","8cell","icm")
# tmp.files.1 <- paste0(tmp.dir,tmp.samples,"_peaks_",tmp.TF.name,"_bindsite.bed")
# tmp.files.2 <- paste0("res/meme/res/bed/",tmp.samples,"_peaks.bed")
# 
# tmp.res.list <- lapply(1:length(tmp.samples),FUN = function(ii){
#   ii <- 2
#   tmp.TF.sites <- import.bed(tmp.files.1[ii])  
#   tmp.ATAC.sites <- import.bed(tmp.files.2[ii])
#   ### TFBS
#   TF_binding_gene <- subsetByOverlaps(mm9.tss,tmp.TF.sites,type = "any")
#   ATAC_binding_gene <- subsetByOverlaps(mm9.tss,tmp.ATAC.sites,type = "any")
#   
#   TF_binding_gene
#   
#   #export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
#   M <- length(TF_binding_gene)
#   N <- length(ATAC_binding_gene)
#   
#   tmp.ATAC.gene <- subsetByOverlaps(mm9.tss[LRgene,],tmp.ATAC.sites,type = "any")
#   tmp.TF_ATAC.gene <- subsetByOverlaps(tmp.ATAC.gene,TF_binding_gene)
#   
#   
#   k <- length(tmp.TF_ATAC.gene)
#   n <- length(tmp.ATAC.gene)
#   
#   x <- M/N
#   y <- k/n
#   p <- phyper(k-1,M,N-M,n,lower.tail = F)
#   #tmp.df <- data.frame(expected = x, obs = y,pvalue = p,sample=tmp.samples[ii],stringsAsFactors = F)
#   tmp.df <- data.frame(sample=tmp.samples[ii],
#                        value=c(x,y),
#                        pvalue=c(p,p),
#                        group=c("expected","observed"),
#                        set = c("all","LRgene"),
#                        stringsAsFactors = F)
#   return(tmp.df)
# })
# 
# 
# tmp.res.df <- Reduce(rbind,tmp.res.list)

#####----------5. plot TRN network-------------

#####----------5.1 load and export all the networks----------

rtna <- readRDS(file = "res/R/early.rtna_20220316.rds")
eLR_row_annotation_df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds")
eLR_embryo_eLR_TF <- readRDS(file = "res/R/early_embryo_eLR_TF_interaction_information_2022060620.rds")
tmp.TF.annotation <- readRDS(file = "res/R/early_embryo_core_TF_annotation_row_2022032521.rds")
tmp.df.use <- tmp.TF.annotation %>%
  rownames_to_column("TF")

eLR_row_annotation_df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds")
LRpairs <- rownames(eLR_row_annotation_df)
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "-"),FUN = function(x) x[1]))
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "-"),FUN = function(x) x[2]))
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)



length(tmp.df.use$TF)

tmp.TFs <-   c("Atf2","Nfya","Nanog","Ctcf",
               "Obox3","Obox5","Obox6",
               "Yy1","Pou5f1","Zscan4f",
               "Sox2","Sox21","Cdx2",
               "Gata4","Gata6","Gata3",
               "Tead1","Tead4","Klf6")


tmp.ttt <- tna.get(rtna,what = "regulons.and.mode")

tTFs <- tmp.df.use$TF

tmp.res.list <- lapply(tTFs,FUN = function(ii){
  cat(ii,sep = "\n")
  x <- tmp.ttt[[ii]]
  tmp.res <- data.frame(source = ii,
                        target = names(x),
                        value = x,
                        stringsAsFactors = F)
  return(tmp.res)
})


tmp.res.df <- Reduce(rbind,tmp.res.list,)

tmp.res.df <- tmp.res.df %>%
  filter(target %in% c(LRgene,tTFs))

tmp.res.df <- merge(tmp.res.df,
                    tmp.TF.meta,
                    by.x="source",
                    by.y="gene")
tmp.res.df$Cluster <- paste0("C",tmp.res.df$Cluster)
tmp.edge <- tmp.res.df
tmp.node <- data.frame(gene=c(Lgene,Rgene),
                       group=c(rep("L",
                                   times = length(Lgene)),
                               rep("R",times = length(Rgene))),
                       Cluster=c(rep("L",times = length(Lgene)),
                               rep("R",times = length(Rgene))),
                       stringsAsFactors = F)

tmp.node.1 <- tmp.TF.meta %>%
  mutate(group = "TF",
         Cluster = paste0("C",Cluster)) %>%
  dplyr::select(gene,group,Cluster)

tmp.node <- rbind(tmp.node,tmp.node.1)



write.table(tmp.edge,
            file = myFileName(prefix = "res/cytoscape/fig5_TRN_eLR_TF_edge",suffix = ".txt"),
            sep = "\t",
            row.names = F,col.names = T,quote = F)
write.table(tmp.node,
            file = myFileName(prefix = "res/cytoscape/fig5_TRN_eLR_TF_node",suffix = ".txt"),
            sep = "\t",
            row.names = F,col.names = T,quote = F)

# dir.create("res/biotapesy")
# write.table(tmp.edge,
#             file = myFileName(prefix = "res/biotapesy/fig5_TRN_eLR_TF_edge",suffix = ".csv"),
#             sep = ",",
#             row.names = F,col.names = F,quote = F)
# write.table(tmp.node,
#             file = myFileName(prefix = "res/biotapesy/fig5_TRN_eLR_TF_node",suffix = ".csv"),
#             sep = ",",
#             row.names = F,col.names = F,quote = F)

####---------5.2 select sub networks------------

weighted_networks <- readRDS(file = "res/R/nichenet_mouse_weighted_network_20210131.rds")


tmp.TFs <-   c("Atf2","Nfya","Nanog","Ctcf",
               "Obox3","Obox5","Obox6",
               "Yy1","Pou5f1","Zscan4f",
               "Sox2","Sox21","Cdx2",
               "Gata4","Gata6","Gata3",
               "Tead1","Tead4","Klf6")

tmp.sig <- weighted_networks$lr_sig %>%
  filter(from %in% c("Fgf4","Fgfr1","Fgfr2"),
         to %in% tmp.TFs)


tmp.edge <- tmp.res.df %>%
  dplyr::filter(source %in% tmp.TFs,
                target %in% c(tmp.TFs,
                              c("Fgf4","Fgfr1","Fgfr2","Bmp4","Bmpr2")))

tmp.idx <- grep(pattern = "Fgf|Bmp",
                x = tmp.edge$target,
                value = T)


tmp.node.TF <- tmp.TF.meta %>%
  filter(gene %in% tmp.TFs)


tmp.edge.1 <- eLR_embryo_eLR_TF %>%
  dplyr::filter(LRpairs %in% c("Fgf4_Fgfr1",
                               "Fgf4_Fgfr2",
                               "Bmp4_Bmpr2")) %>%
  dplyr::filter(group != "background") %>%
  dplyr::select(Lgene,Rgene,TFcluster,group) %>%
  mutate(TFcluster = gsub(pattern = "TF_cluster_",
                          replacement = "C",
                          TFcluster))




write.table(tmp.edge,
            file = myFileName(prefix = "res/cytoscape/fig5_TRN_eLR_TF_edge_subnetworks",
                              suffix = ".txt"),
            sep = "\t",
            row.names = F,col.names = T,quote = F)



