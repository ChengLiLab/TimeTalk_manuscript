####This scritpt was to address the R1.4
####To test the sensitivity of TimeTalk
####Author: wlt
####Date: 20230502


####--------0. load packages------------
library(Seurat)
#library(monocle)
library(monocle3)
library(DropletUtils)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(cowplot)
library(RColorBrewer)
## it's better to install ComplexHeatmap for github
library(ComplexHeatmap)
library(ggthemes)
library(CellChat)
library(future.apply)
library(GEOquery)
library(RTN)
library(RTNduals)
library(snow)
#library(Rmagic)
library(netSmooth)
library(circlize)
library(future.apply)
library(cellAlign)
library(lmtest)
library(AUCell)
library(eulerr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(TimeTalk)
library(tidyverse)

source(file = "code/myUtils.R")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
options(future.globals.maxSize= 10*1024^3)

###--------1.Test TimeResults-------------
####-------1.1 load data---------------------
tmp.mra.res <- readRDS(file = "res/R/B_blastoid_RTN_mra_result.rds")
LRpairs.df <- read.delim(file = "database/Ligand-Receptor-Pairs/Mouse/Mouse-2020-Shao-LR-pairs.txt",stringsAsFactors = F)

####load sc object
seu <- readRDS(file = "res/R/B_blastoid_seurat_inter_2_data_2022042218.rds")
cds <- readRDS(file = "res/R/B_blastoid_monocle3_cds_2022042219.rds")
####check cell types
seu$CellType <- Idents(seu)


#DimPlot(seu)

TimeTalk.result <- RunTimeTalk(tmp.cds=cds,
                               tmp.seu=seu,
                               tmp.orig.ident = "blastocyst",
                               tmp.ident.1="EPI",
                               tmp.ident.2 = "PE",
                               LRpairs.df = LRpairs.df,
                               tmp.mra.res = tmp.mra.res,
                               tmp.winsz = 0.1,
                               tmp.lags = 1,
                               numPts = 200,
                               tmp.SCC.cutoff = 0.2,
                               tmp.granger.cutoff = 1e-2)

tmp.res <- TimeTalk.result %>%
  filter(category == "PASS") %>%
  pull(LR) %>%
  unique()

####----------1.2 Robust test---------------------------

####----------1.2.1 Test window size-----------------

####----------1.2.1.1 Run on different parameters------------
test_winsz <- lapply(seq(0.1,0.5,by=0.1),
                            FUN = function(x){
                              TimeTalk.result <- RunTimeTalk(tmp.cds=cds,
                                                             tmp.seu=seu,
                                                             tmp.orig.ident = "blastocyst",
                                                             tmp.ident.1="EPI",
                                                             tmp.ident.2 = "PE",
                                                             LRpairs.df = LRpairs.df,
                                                             tmp.mra.res = tmp.mra.res,
                                                             tmp.winsz = x,
                                                             tmp.lags = 1,
                                                             numPts = 200,
                                                             tmp.cores = 3,
                                                             tmp.SCC.cutoff = 0.2,
                                                             tmp.granger.cutoff = 1e-2)
                              
                              tmp.res <- TimeTalk.result 
                              return(tmp.res)
                            })

saveRDS(object = test_winsz,file = "res/R/test_winsz_20230502.rds")

###--------1.2.1.2 calculate overlap ratio-----------



####load previsous results
tmp.TimeTalk.res <- readRDS(file = "res/R/TimeTalk_res_B_blastoid_EPI_PE_2022043013.rds") 
tmp.plot.1 <- tmp.TimeTalk.res %>%
  filter(orig.ident == "blastocyst") %>%
  group_by(LR,orig.ident) %>%
  summarise(LRtoTF.ens = min(LRtoTF),
            TFtoLR.ens = min(TFtoLR),
            PCC.ens = max(PCC),
            SCC.ens = max(SCC)) %>%
  ungroup() %>%
  #filter(abs(SCC.ens) > 0.8 ) %>%
  filter(abs(SCC.ens) > 0.8 ) %>%
  #filter((LRtoTF.ens < 1e-20) | (TFtoLR.ens > 1e-20)) %>%
  arrange(-SCC.ens) %>%
  mutate(Rank = row_number())

####tmp_gold standard
tmp_gold_LR <- tmp.plot.1$LR


test_winsz <- readRDS(file = "res/R/test_winsz_20230502.rds")
test_winsz_processed <- lapply(test_winsz, function(tmp_df){
  tmp_df_res <- tmp_df %>%
    dplyr::filter(category == "PASS") %>%
    group_by(LR) %>%
    summarise(LRtoTF.ens = min(LRtoTF),
              TFtoLR.ens = min(TFtoLR),
              PCC.ens = max(PCC),
              SCC.ens = max(SCC)) %>%
    ungroup() %>%
    #filter(abs(SCC.ens) > 0.8 ) %>%
    filter(abs(SCC.ens) > 0.8 ) %>%
    #filter((LRtoTF.ens < 1e-20) | (TFtoLR.ens > 1e-20)) %>%
    arrange(-SCC.ens) %>%
    mutate(Rank = row_number())
  return(tmp_df_res)
})

### winsz = 0.1 is the default results
identical(test_winsz_processed[[1]]$LR,tmp_gold_LR)

tmp.data.plot <- lapply(test_winsz, function(tmp_df){
  tmp.res <- tmp_df %>%
    dplyr::filter(category == "PASS") %>%
    group_by(LR) %>%
    summarise(LRtoTF.ens = min(LRtoTF),
              TFtoLR.ens = min(TFtoLR),
              PCC.ens = max(PCC),
              SCC.ens = max(SCC)) %>%
    ungroup() %>%
    #filter(abs(SCC.ens) > 0.8 ) %>%
    filter(abs(SCC.ens) > 0.8 ) %>%
    #filter((LRtoTF.ens < 1e-20) | (TFtoLR.ens > 1e-20)) %>%
    arrange(-SCC.ens) %>%
    mutate(Rank = row_number()) %>%
    pull(LR)
  tmp.res.ratio <- length(intersect(tmp.res,tmp_gold_LR)) / length(tmp_gold_LR)
  return(tmp.res.ratio)
})

tmp.data.plot <- data.frame(mywinsz = seq(0.1,0.5,by=0.1),
                            ratio = unlist(tmp.data.plot),
                            stringsAsFactors = F)

ggplot(tmp.data.plot,
       aes(mywinsz,ratio,
           fill = as.character(mywinsz),
           label = sprintf("%.1f%%", signif(ratio*100,3))))+
  geom_bar(stat = "identity",color="black",width = 0.08)+
  geom_text(vjust = -0.5,size = 6)+
  scale_fill_manual(values = blues(5))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
  xlab("winsz")+
  ylab("Percent of overlaped interactions\n(compared to winsz = 0.1)")+
  theme_cowplot(font_size = 22)+
  theme(legend.position = "none",
        axis.ticks = element_line(size = 1),
        axis.line = element_line(size = 1))
ggsave(filename = myFileName(prefix = "res/fig/FigS_revision_robust_winsz",
                             suffix = ".jpg"),
       width = 6,height = 6,dpi = 350)
ggsave(filename = myFileName(prefix = "res/fig/FigS_revision_robust_winsz",
                             suffix = ".pdf"),
       width = 6,height = 6)

####----------1.2.2 Test the remaining parameters-------

####----------1.2.2.1 Perform analysis-----------
test_lag <- lapply(seq(1,5,by=1),
                     FUN = function(x){
                       TimeTalk.result <- RunTimeTalk(tmp.cds=cds,
                                                      tmp.seu=seu,
                                                      tmp.orig.ident = "blastocyst",
                                                      tmp.ident.1="EPI",
                                                      tmp.ident.2 = "PE",
                                                      LRpairs.df = LRpairs.df,
                                                      tmp.mra.res = tmp.mra.res,
                                                      tmp.winsz = 0.1,
                                                      tmp.lags = x,
                                                      numPts = 200,
                                                      tmp.cores = 3,
                                                      tmp.SCC.cutoff = 0.2,
                                                      tmp.granger.cutoff = 1e-2)
                       
                       tmp.res <- TimeTalk.result 
                       return(tmp.res)
                     })

saveRDS(object = test_lag,
        file = "res/R/test_lag_20230502.rds")

###Change the Point Size
test_Pt <- lapply(seq(100,500,by=100),
                   FUN = function(x){
                     TimeTalk.result <- RunTimeTalk(tmp.cds=cds,
                                                    tmp.seu=seu,
                                                    tmp.orig.ident = "blastocyst",
                                                    tmp.ident.1="EPI",
                                                    tmp.ident.2 = "PE",
                                                    LRpairs.df = LRpairs.df,
                                                    tmp.mra.res = tmp.mra.res,
                                                    tmp.winsz = 0.1,
                                                    tmp.lags = 1,
                                                    numPts = x,
                                                    tmp.cores = 3,
                                                    tmp.SCC.cutoff = 0.2,
                                                    tmp.granger.cutoff = 1e-2)
                     
                     tmp.res <- TimeTalk.result 
                     return(tmp.res)
                   })

saveRDS(object = test_Pt,
        file = "res/R/test_Pt_20230502.rds")

###--------1.2.2.2 plot resulst---------------

###--------1.2.2.2.1 Plot lags-------------
test_lag <- readRDS(file = "res/R/test_lag_20230502.rds")

tmp.data.plot <- lapply(test_lag, function(tmp_df){
  tmp.res <- tmp_df %>%
    dplyr::filter(category == "PASS") %>%
    group_by(LR) %>%
    summarise(LRtoTF.ens = min(LRtoTF),
              TFtoLR.ens = min(TFtoLR),
              PCC.ens = max(PCC),
              SCC.ens = max(SCC)) %>%
    ungroup() %>%
    #filter(abs(SCC.ens) > 0.8 ) %>%
    filter(abs(SCC.ens) > 0.8 ) %>%
    #filter((LRtoTF.ens < 1e-20) | (TFtoLR.ens > 1e-20)) %>%
    arrange(-SCC.ens) %>%
    mutate(Rank = row_number()) %>%
    pull(LR)
  tmp.res.ratio <- length(intersect(tmp.res,tmp_gold_LR)) / length(tmp_gold_LR)
  return(tmp.res.ratio)
})

head(tmp.data.plot)

tmp.data.plot <- data.frame(para_list = seq(1,5,by=1),
                            ratio = unlist(tmp.data.plot),
                            stringsAsFactors = F)
ggplot(tmp.data.plot,
       aes(para_list,ratio,
           fill = as.character(para_list),
           label = sprintf("%.1f%%", signif(ratio*100,3))))+
  geom_bar(stat = "identity",color="black",width = 0.8)+
  geom_text(vjust = -0.5,size = 6)+
  scale_fill_manual(values = blues(5))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
  xlab("lags")+
  ylab("Percent of overlaped interactions\n(compared to lags = 1)")+
  theme_cowplot(font_size = 22)+
  theme(legend.position = "none",
        axis.ticks = element_line(size = 1),
        axis.line = element_line(size = 1))
ggsave(filename = myFileName(prefix = "res/fig/FigS_revision_robust_lags",
                             suffix = ".jpg"),
       width = 6,height = 6,dpi = 350)
ggsave(filename = myFileName(prefix = "res/fig/FigS_revision_robust_lags",
                             suffix = ".pdf"),
       width = 6,height = 6)



###--------1.2.2.2.2 Plot WinSz-------------
test_Pt <- readRDS(file = "res/R/test_Pt_20230502.rds")

tmp.data.plot <- lapply(test_Pt, function(tmp_df){
  tmp.res <- tmp_df %>%
    dplyr::filter(category == "PASS") %>%
    group_by(LR) %>%
    summarise(LRtoTF.ens = min(LRtoTF),
              TFtoLR.ens = min(TFtoLR),
              PCC.ens = max(PCC),
              SCC.ens = max(SCC)) %>%
    ungroup() %>%
    #filter(abs(SCC.ens) > 0.8 ) %>%
    filter(abs(SCC.ens) > 0.8 ) %>%
    #filter((LRtoTF.ens < 1e-20) | (TFtoLR.ens > 1e-20)) %>%
    arrange(-SCC.ens) %>%
    mutate(Rank = row_number()) %>%
    pull(LR)
  tmp.res.ratio <- length(intersect(tmp.res,tmp_gold_LR)) / length(tmp_gold_LR)
  return(tmp.res.ratio)
})

tmp.data.plot

tmp.data.plot <- data.frame(para_list = seq(100,500,by=100),
                            ratio = unlist(tmp.data.plot),
                            stringsAsFactors = F)
ggplot(tmp.data.plot,
       aes(para_list,ratio,
           fill = as.character(para_list),
           label = sprintf("%.1f%%", signif(ratio*100,3))))+
  geom_bar(stat = "identity",color="black",width = 80)+
  geom_text(vjust = -0.5,size = 6)+
  scale_fill_manual(values = blues(5))+
  scale_y_continuous(expand = c(0,0),limits = c(0,1.1))+
  xlab("numPts")+
  ylab("Percent of overlaped interactions\n(compared to numPts = 100)")+
  theme_cowplot(font_size = 22)+
  theme(legend.position = "none",
        axis.ticks = element_line(size = 1),
        axis.line = element_line(size = 1))

ggsave(filename = myFileName(prefix = "res/fig/FigS_revision_robust_numPts",
                             suffix = ".jpg"),
       width = 6,height = 6,dpi = 350)
ggsave(filename = myFileName(prefix = "res/fig/FigS_revision_robust_numPts",
                             suffix = ".pdf"),
       width = 6,height = 6)


####--------1.3 Tune Multiple Parameters--------------


x <- seq(0.1,0.5,by = 0.1)
y <- seq(100,500,by = 100)
# x <- seq(0.1,0.2,by = 0.1)
# y <- seq(100,200,by = 100)

outer(0:5, 0:6, Vectorize(function(x,y)sum(x,y)))



### We need to vectorize the funciton to perform elementwise evaluation
### For example, this is right
### outer(0:5, 0:6, Vectorize(function(x,y)sum(x,y)))
### But this is wrong
### outer(0:5, 0:6, sum)

.temp_overlap_ratio <- function(x,y,tmp.lags=1){
  TimeTalk.result <- RunTimeTalk(tmp.cds=cds,
                                 tmp.seu=seu,
                                 tmp.orig.ident = "blastocyst",
                                 tmp.ident.1="EPI",
                                 tmp.ident.2 = "PE",
                                 LRpairs.df = LRpairs.df,
                                 tmp.mra.res = tmp.mra.res,
                                 tmp.winsz = x,
                                 tmp.lags = tmp.lags,
                                 numPts = y,
                                 tmp.cores = 5,
                                 tmp.SCC.cutoff = 0.2,
                                 tmp.granger.cutoff = 1e-2)
  
  tmp.res <- TimeTalk.result %>%
    dplyr::filter(category == "PASS") %>%
    group_by(LR) %>%
    summarise(LRtoTF.ens = min(LRtoTF),
              TFtoLR.ens = min(TFtoLR),
              PCC.ens = max(PCC),
              SCC.ens = max(SCC)) %>%
    ungroup() %>%
    #filter(abs(SCC.ens) > 0.8 ) %>%
    filter(abs(SCC.ens) > 0.8 ) %>%
    #filter((LRtoTF.ens < 1e-20) | (TFtoLR.ens > 1e-20)) %>%
    arrange(-SCC.ens) %>%
    mutate(Rank = row_number()) %>%
    pull(LR)
  tmp.res.ratio <- length(intersect(tmp.res,tmp_gold_LR)) / length(tmp_gold_LR)
  return(tmp.res.ratio)
}
ttt.ratio <- .temp_overlap_ratio(x = 0.01,y = 200,tmp.lags = 1)




Vectorize(function(x,y) .temp_overlap_ratio(x,y))

x <- seq(0.1,0.5,by = 0.1)
y <- seq(100,500,by = 100)

tmp.res.list <- lapply(1:5,FUN = function(ii){
  tmp.res <- outer(x,y,Vectorize(function(x,y) 
    .temp_overlap_ratio(x,y,tmp.lags = ii)))
  rownames(tmp.res) <- x
  colnames(tmp.res) <- y
  return(tmp.res)
})
saveRDS(tmp.res.list,file = "res/R/FigS_revision_vary_tags_20230503.rds")
head(tmp.res.list[[1]])


# ht <- pheatmap_fixed(tmp.res*100,
#                      color = viridisLite::viridis(n = 100),
#                      cluster_cols = F,
#                      cluster_rows = F,
#                      number_color = "white",
#                      display_numbers = signif(tmp.res*100,digits = 4),
#                      number_format = "%.2f",
#                      fontsize = 16,
#                      fontsize_number = 20,
#                      border_color = "black",
#                      row_title = "winsz",
#                      row_title_gp = gpar(fontsize = 18,fontface = "bold"),
#                      column_title = "lag=1\nnumPts",
#                      column_title_side = "top",
#                      column_title_gp = gpar(fontsize = 18,fontface = "bold"),
#                      name = "Percent %")


tmp.res.list[[1]]



ht.list <- lapply(1:5,FUN = function(ii){
  ht <- pheatmap_fixed(tmp.res.list[[ii]]*100,
                       color = viridisLite::viridis(n = 100),
                       cluster_cols = F,
                       cluster_rows = F,
                       number_color = "white",
                       display_numbers = signif(tmp.res.list[[ii]]*100,digits = 4),
                       number_format = "%.2f",
                       fontsize = 16,
                       fontsize_number = 20,
                       border_color = "black",
                       row_title = "winsz",
                       row_title_gp = gpar(fontsize = 18,fontface = "bold"),
                       column_title = paste0("lag=",ii,"\nnumPts"),
                       column_title_side = "top",
                       column_title_gp = gpar(fontsize = 18,fontface = "bold"),
                       name = "Percent %")
  return(ht)
})

draw(Reduce("+",ht.list))

pdf(file = myFileName(prefix = "res/fig/FigS_revision_robust_lags_vary_heatmap",
                      suffix = ".pdf"),
    width = 24,
    height = 6)
draw(Reduce("+",ht.list))
dev.off()

jpeg(filename = myFileName(prefix = "res/fig/FigS_revision_robust_lags_vary_heatmap",
                      suffix = ".jpg"),
    width = 24,
    height = 6,
    units = "in",
    res = 350)
draw(Reduce("+",ht.list))
dev.off()


jpeg(filename = myFileName(prefix = "res/fig/FigS_revision_robust_lags_1_heatmap",
                           suffix = ".jpg"),
     width = 6,
     height = 6,
     units = "in",
     res = 350)
draw(ht)
dev.off()


####---------2. check trajectory------------

####---------2.1 visualization----------------

cds <- readRDS(file = "res/R/B_blastoid_monocle3_cds_2022042219.rds")

p1 <- plot_cells(cds, 
                 color_cells_by = "pseudotime", 
                 cell_size = 1.5,
                 group_label_size = 9,
                 label_principal_points = F,
                 label_groups_by_cluster = F, 
                 label_leaves = F, 
                 label_roots = T,
                 label_branch_points = F)+
  scale_color_gradientn(name = "pseudotime",colours = warm(100))+
  theme_cowplot(font_size = 28) +
  NoAxes()+
  theme(legend.justification = "center",
        aspect.ratio = 1,
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
p1

tmp.colors <- c("#F2BE58",divergentcolor(7)[1:4],
                         grey.colors(3,start = 0.5,rev = T))
set.seed(42)
tmp.cell.id <- sample(rownames(pData(cds)))
p2 <- plot_cells(cds[,tmp.cell.id], 
                 color_cells_by = "CellType", 
                 cell_size = 1.5,
                 group_label_size = 9,
                 label_principal_points = F,
                 label_groups_by_cluster = FALSE, 
                 label_leaves = FALSE, 
                 label_roots = F,
                 label_branch_points = F)+
  scale_color_manual(values = tmp.colors)+
  theme_cowplot(font_size = 28) +
  NoAxes()+
  theme(legend.position = "none",
        legend.justification = "center",
        aspect.ratio = 1,
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))

p2 | p1
ggsave(filename = myFileName(prefix = "res/fig/fig5_moncle3_trajectory_show_root_cell",
                             suffix = ".png"),
       width = 12,height = 6,dpi = 350)

devtools::install_local("/mnt/d/softwarepackage/CytoTRACE_0.3.3.tar.gz")


###------2.2 try ctyotrace-------------------
###------2.2.1 load package---------

library(reticulate)
use_condaenv("cytotraceR")
library(TimeTalk)
library(CytoTRACE)
library(Seurat)
library(tidyverse)
library(monocle3)
source("code/myUtils.R")
library(cowplot)
# Test data
# datasets <- list(marrow_10x_expr, marrow_plate_expr)
# results <- iCytoTRACE(datasets)
# ### The resulsts is a list 
# 
# head(marrow_10x_expr[1:3,1:3])
# head(marrow_plate_expr[1:3,1:3])
###------2.2.2 prepare datasets-----------

seu <- readRDS(file = "res/R/B_blastoid_seurat_inter_2_data_2022042218.rds")
cds <- readRDS(file = "res/R/B_blastoid_monocle3_cds_2022042219.rds")

###add pseudotime to metadata of seurat
tmp.df <- data.frame(pseudotime = pseudotime(cds,
                                             reduction_method = "UMAP"),
                     stringsAsFactors = F)
identical(rownames(seu[[]]),rownames(tmp.df))
seu <- AddMetaData(seu,metadata = tmp.df)


seu.list <- SplitObject(seu,split.by = "orig.ident")
mat_blastocyst <- GetAssayData(seu.list$blastocyst,slot = "counts")
mat_blastoid <- GetAssayData(seu.list$EPS_blastoid,slot = "counts")

mat_blastocyst <- as.data.frame(mat_blastocyst)
mat_blastoid <- as.data.frame(mat_blastoid)

datasets <- list(mat_blastocyst,mat_blastoid)

###-----2.2.3 run cytotrace-----------
system.time(results <- iCytoTRACE(datasets = datasets))
# user  system elapsed 
# 56.226  20.955  50.082 
system.time(results.blastocyst <- CytoTRACE(mat = mat_blastocyst))

tmp.df <- data.frame(CytoTRACE = results$CytoTRACE,
                     CytoTRACErank = results$CytoTRACErank)
head(tmp.df)
seu <- AddMetaData(seu,metadata = tmp.df)
saveRDS(seu,file = "res/R/B_blastoid_seurat_add_cytotrace_20230508.rds")

seu.blastocyst <- seu.list$blastocyst
tmp.df <- data.frame(CytoTRACE = results.blastocyst$CytoTRACE,
                     CytoTRACErank = results.blastocyst$CytoTRACErank)
seu.blastocyst <- AddMetaData(seu.blastocyst,metadata = tmp.df)

FeaturePlot(seu.blastocyst,features = "CytoTRACE")

###-----2.3. visualization of results-------------


###-----2.3.1 load data ----------------
seu <- readRDS(file = "res/R/B_blastoid_seurat_add_cytotrace_20230508.rds")

FeaturePlot(seu,features = "CytoTRACE")
seu$pseudotime <- MinMaxScale(seu$pseudotime)

FeaturePlot(seu,features = c("pseudotime","CytoTRACE")) &
  scale_color_gradientn(name = "pseudotime",colours = warm(100))
ggsave(filename = myFileName(prefix = "res/fig/figS_revsion_cytotrace_validation_of_trajectory",
                             suffix = ".jpg"),
       width = 12,height = 6)

head(seu[[]])

tmp.cells.1 <- seu[[]] %>%
  rownames_to_column("cell_id") %>%
  dplyr::filter(orig.ident == "blastocyst",CellType == "EPI") %>%
  dplyr::select(cell_id,pseudotime,CytoTRACE) %>%
  dplyr::mutate(pseudotime_rank = row_number(pseudotime),
                CytoTRACE_rank = row_number(CytoTRACE))
tmp.cells.2 <- seu[[]] %>%
  rownames_to_column("cell_id") %>%
  dplyr::filter(orig.ident == "blastocyst",CellType == "PE") %>%
  dplyr::select(cell_id,pseudotime,CytoTRACE) %>%
  dplyr::mutate(pseudotime_rank = row_number(pseudotime),
                CytoTRACE_rank = row_number(CytoTRACE))

p1 <- ggplot(tmp.cells.1,aes(pseudotime_rank,CytoTRACE_rank))+
  ggtitle("The rank of EPI cells")+
  geom_line()+
  theme_cowplot(font_size = 22)
p2 <- ggplot(tmp.cells.2,aes(pseudotime_rank,CytoTRACE_rank))+
  ggtitle("The rank of PE cells")+
  geom_line()+
  theme_cowplot(font_size = 22)

p1 | p2
ggsave(filename = myFileName(prefix = paste0("res/fig/Figs_revision_cytotrace_rank_of_cells"),
                             suffix = ".jpg"),
       width = 12,
       height = 6,
       dpi = 450)


####-----2.3.2 show cases of the exist causal relationship---------
#### You should believe your hope and capacity

tmp.df <- seu[[]]

tmp.mra.res <- readRDS(file = "res/R/B_blastoid_RTN_mra_result.rds")
LRpairs.df <- read.delim(file = "database/Ligand-Receptor-Pairs/Mouse/Mouse-2020-Shao-LR-pairs.txt",stringsAsFactors = F)
seu <- readRDS(file = "res/R/B_blastoid_seurat_inter_2_data_2022042218.rds")
cds <- readRDS(file = "res/R/B_blastoid_monocle3_cds_2022042219.rds")


# x.df <- tmp.df %>%
#   dplyr::filter(orig.ident == "blastocyst" & CellType == "EPI") %>%
#   dplyr::select(pseudotime,CytoTRACE)
# 
# y.df <- tmp.df %>%
#   dplyr::filter(orig.ident == "blastocyst" & CellType == "PE") %>%
#   dplyr::select(pseudotime,CytoTRACE)
# 
# ggplot(x.df,aes(pseudotime,CytoTRACE))+
#   geom_point()
# 
# ggplot(y.df,aes(pseudotime,CytoTRACE))+
#   geom_point()
# 
# cor(y.df$pseudotime,y.df$CytoTRACE)
# TimeTalk.result <- RunTimeTalk(tmp.cds=cds,
#                                tmp.seu=seu,
#                                tmp.orig.ident = "blastocyst",
#                                tmp.ident.1 = "EPI",
#                                tmp.ident.2 = "PE",
#                                LRpairs.df = LRpairs.df,
#                                tmp.mra.res = tmp.mra.res,
#                                tmp.winsz = 0.1,
#                                tmp.lags = 1,
#                                numPts = 200,
#                                tmp.cores = 3,
#                                tmp.SCC.cutoff = 0.2,
#                                tmp.granger.cutoff = 1e-2)
###------2.3.2.1 let me run the data ----------
tmp.cds=cds
tmp.seu=seu
tmp.orig.ident = "blastocyst"
tmp.ident.1 = "EPI"
tmp.ident.2 = "PE"
LRpairs.df = LRpairs.df
tmp.mra.res = tmp.mra.res
tmp.winsz = 0.1
tmp.lags = 1
numPts = 200
tmp.cores = 3
tmp.SCC.cutoff = 0.2
tmp.granger.cutoff = 1e-2


cat(paste0("prepare ligand and receptor genelist"), sep = "\n")
LRpairs <- LRpairs.df$lr_pair
Lgenelist <- LRpairs.df$ligand_gene_symbol
Rgenelist <- LRpairs.df$receptor_gene_symbol

### using data in RNA assay
tmp.data <- GetAssayData(object = seu, slot = "data", assay = "RNA")
gene_symbols <- rownames(tmp.data)
l.remove <- setdiff(Lgenelist, gene_symbols)
r.remove <- setdiff(Rgenelist, gene_symbols)
index.remove <- c(which(Lgenelist %in% l.remove), which(Rgenelist %in% r.remove))
LRpairs <- LRpairs[-index.remove]
Lgenelist <- Lgenelist[-index.remove]
Rgenelist <- Rgenelist[-index.remove]


cat(paste0(tmp.ident.1, "-", tmp.ident.2, " start:"), sep = "\n")
tmp.df <- data.frame(
  pseudotime = pseudotime(cds, reduction_method = "UMAP"),
  stringsAsFactors = F
)
tmp.df <- cbind(seu[[]], tmp.df)
### Check CellTypes
if (!c("CellType") %in% colnames(tmp.df)) {
  stop("Please add CellType annotation in seurat object!")
}

if (!c("orig.ident") %in% colnames(tmp.df)) {
  stop("Please add CellType annotation in seurat object!")
}

cat(paste0("scale pseudotime"), sep = "\n")
tmp.cell.meta.1 <- tmp.df %>%
  rownames_to_column("cell_id") %>%
  dplyr::filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.1) %>%
  arrange(pseudotime)
x <- tmp.cell.meta.1$pseudotime
tmp.cell.meta.1$pseudotime <- MinMaxScale(x)

tmp.cell.meta.2 <- tmp.df %>%
  rownames_to_column("cell_id") %>%
  dplyr::filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.2) %>%
  arrange(pseudotime)
x <- tmp.cell.meta.2$pseudotime
tmp.cell.meta.2$pseudotime <- MinMaxScale(x)

tmp.mat.1 <- tmp.data[, tmp.cell.meta.1$cell_id]
tmp.mat.2 <- tmp.data[, tmp.cell.meta.2$cell_id]
head(tmp.mat.1)


tmp.mat.pseudotime.1 <- tmp.cell.meta.1$pseudotime
tmp.mat.pseudotime.2 <- tmp.cell.meta.2$pseudotime

names(tmp.mat.pseudotime.1) <- tmp.cell.meta.1$cell_id
names(tmp.mat.pseudotime.2) <- tmp.cell.meta.2$cell_id

inter.tmp.mat.1 <- cellAlign::interWeights(
  expDataBatch = tmp.mat.1,
  trajCond = tmp.mat.pseudotime.1,
  winSz = tmp.winsz,
  numPts = numPts
)
inter.tmp.mat.2 <- cellAlign::interWeights(
  expDataBatch = tmp.mat.2,
  trajCond = tmp.mat.pseudotime.2,
  winSz = tmp.winsz,
  numPts = numPts
)

inter.tmp.mat.1 <- cellAlign::scaleInterpolate(inter.tmp.mat.1)
inter.tmp.mat.2 <- cellAlign::scaleInterpolate(inter.tmp.mat.2)
time <- inter.tmp.mat.1$traj

inter.tmp.mat.1 <- myRemoveNA(inter.tmp.mat.1$scaledData)
inter.tmp.mat.2 <- myRemoveNA(inter.tmp.mat.2$scaledData)

head(inter.tmp.mat.1)


###-------2.3.2.2 Prepare data to plot-----------

tmp.L.gene <- "Fgf4"
tmp.R.gene <- "Fgfr2"
tmp.LR <- "Fgf4-Fgfr2"
tmp.TF.gene.use <- "Gata6"

x.TF <- "Nanog"
y.TF <- "Gata6"


x <- inter.tmp.mat.1[tmp.L.gene,]
y <- inter.tmp.mat.2[tmp.R.gene,]

tmp.TF.level <- inter.tmp.mat.2[tmp.TF.gene.use, ]
tmp.IS <- sqrt(x * y)
tmp.res.LRtoTF.pvalue <- tryCatch(
  expr = {
    tmp.res <- grangertest(tmp.IS, tmp.TF.level, order = tmp.lags)
    tmp.res$`Pr(>F)`[2]
  },
  error = function(e) {
    1
  }
)

tmp.res.TFtoLR.pvalue <- tryCatch(
  expr = {
    tmp.res <- grangertest(tmp.TF.level, tmp.IS, order = tmp.lags)
    tmp.res$`Pr(>F)`[2]
  },
  error = function(e) {
    1
  }
)

### IS and TF
tmp.PCC <- cor(tmp.IS, tmp.TF.level, method = "pearson")
tmp.PCC <- ifelse(is.na(tmp.PCC), 0, tmp.PCC)
tmp.SCC <- cor(tmp.IS, tmp.TF.level, method = "spearman")
tmp.SCC <- ifelse(is.na(tmp.SCC), 0, tmp.SCC)

###--------2.3.2.3 L data to plot ----------
LR.color <- c("#ffc000","#00b050","#EE0000FF","#631879FF")
names(LR.color) <- c("L","R","TF","IS")

x.TF
x.TF.level <- inter.tmp.mat.1[x.TF,]
tmp.xxx <- cor(x.TF.level,
               x, 
               method = "spearman")
tmp.xxx <- ifelse(is.na(tmp.xxx), 0,tmp.xxx)
tmp.SCC.TF.to.L <- tmp.xxx



tmp.res.LtoTF.pvalue <- tryCatch(
  expr = {
    tmp.res <- grangertest(x, x.TF.level, 
                           order = tmp.lags)
    tmp.res$`Pr(>F)`[2]
  },
  error = function(e) {
    1
  }
)

tmp.res.TFtoL.pvalue <- tryCatch(
  expr = {
    tmp.res <- grangertest(x.TF.level,x, order = tmp.lags)
    tmp.res$`Pr(>F)`[2]
  },
  error = function(e) {
    1
  }
)
tmp.res.LtoTF.pvalue
tmp.res.TFtoL.pvalue
tmp.lags


tmp.data.plot <- data.frame(Index = 1:length(x),
                            value = c(x,x.TF.level),
                            group = c(base::rep("L",times = length(x)),
                                      base::rep("TF",times = length(x.TF.level))),
                            stringsAsFactors = F)

tmp.data.plot$group <- factor(tmp.data.plot$group,
                              levels = c("L","TF"))

ggplot(tmp.data.plot,aes(Index,
                         value,
                         color = group))+
  geom_point()+
  xlab("T")+
  ggtitle(paste0("TF-L in sender cell","\n",
                 "SCC=",sprintf("%.3f",
                                tmp.SCC.TF.to.L),"\n",
                 "p_L_to_TF=",sprintf("%.3f",tmp.res.LtoTF.pvalue),"\n",
                 "p_TF_to_L=",sprintf("%.3f",tmp.res.TFtoL.pvalue)))+
  scale_color_manual(values = LR.color,
                     name = NULL,
                     labels = c(tmp.L.gene,x.TF))+
  guides(colour = guide_legend(override.aes = list(size=8)))+
  theme_cowplot(font_size = 22)+
  theme(axis.line  = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.position = "bottom",
        legend.justification = "center",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 28))

ggsave(filename = myFileName(prefix = paste0("res/fig/figS_revision_R1_4_",tmp.L.gene,"_",x.TF,"_"),
                             suffix = ".jpg"),
       width = 6,height = 6,dpi = 350)

ggsave(filename = myFileName(prefix = paste0("res/fig/figS_revision_R1_4_",tmp.L.gene,"_",x.TF,"_"),
                             suffix = ".pdf"),
       width = 6,height = 6)


####
###--------2.3.2.4 R data to plot ----------
LR.color <- c("#ffc000","#00b050","#EE0000FF","#631879FF")
names(LR.color) <- c("L","R","TF","LR")

y.TF
y.TF.level <- inter.tmp.mat.2[y.TF,]
tmp.xxx <- cor(y.TF.level,
               y, 
               method = "spearman")
tmp.xxx <- ifelse(is.na(tmp.xxx), 0,tmp.xxx)
tmp.SCC.TF.to.R <- tmp.xxx


tmp.res.RtoTF.pvalue <- tryCatch(
  expr = {
    tmp.res <- grangertest(y, y.TF.level, 
                           order = tmp.lags)
    tmp.res$`Pr(>F)`[2]
  },
  error = function(e) {
    1
  }
)

tmp.res.TFtoR.pvalue <- tryCatch(
  expr = {
    tmp.res <- grangertest(y.TF.level,y, order = tmp.lags)
    tmp.res$`Pr(>F)`[2]
  },
  error = function(e) {
    1
  }
)
tmp.res.RtoTF.pvalue
tmp.res.TFtoR.pvalue
tmp.lags


tmp.data.plot <- data.frame(Index = 1:length(y),
                            value = c(y,y.TF.level),
                            group = c(base::rep("R",times = length(y)),
                                      base::rep("TF",times = length(y.TF.level))),
                            stringsAsFactors = F)

tmp.data.plot$group <- factor(tmp.data.plot$group,
                              levels = c("R","TF"))

ggplot(tmp.data.plot,aes(Index,
                         value,
                         color = group))+
  geom_point()+
  xlab("T")+
  ggtitle(paste0("R-TF in receiver cell","\n",
                 "SCC=",sprintf("%.3f",
                                tmp.SCC.TF.to.R),"\n",
                 "p_R_to_TF=",signif(tmp.res.RtoTF.pvalue,digits = 3),"\n",
                 "p_TF_to_R=",signif(tmp.res.TFtoR.pvalue,digits = 3)))+
  scale_color_manual(values = LR.color,
                     name = NULL,
                     labels = c(tmp.R.gene,y.TF))+
  guides(colour = guide_legend(override.aes = list(size=8)))+
  theme_cowplot(font_size = 22)+
  theme(axis.line  = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.position = "bottom",
        legend.justification = "center",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 28))

ggsave(filename = myFileName(prefix = paste0("res/fig/figS_revision_R1_4_",tmp.R.gene,"_",y.TF,"_"),
                             suffix = ".jpg"),
       width = 6,height = 6,dpi = 350)

ggsave(filename = myFileName(prefix = paste0("res/fig/figS_revision_R1_4_",tmp.R.gene,"_",y.TF,"_"),
                             suffix = ".pdf"),
       width = 6,height = 6)



####--------2.3.2.5 LR data to plot---------------

tmp.data.plot <- data.frame(Index = 1:length(tmp.IS),
                            value = c(tmp.IS,y.TF.level),
                            group = c(base::rep("LR",times = length(tmp.IS)),
                                      base::rep("TF",times = length(y.TF.level))),
                            stringsAsFactors = F)

tmp.data.plot$group <- factor(tmp.data.plot$group,
                              levels = c("LR","TF"))

ggplot(tmp.data.plot,aes(Index,
                         value,
                         color = group))+
  geom_point()+
  xlab("T")+
  ggtitle(paste0("LR-TF in receiver cell","\n",
                 "SCC=",sprintf("%.3f",
                                tmp.SCC),"\n",
                 "p_LR_to_TF=",signif(tmp.res.LRtoTF.pvalue,digits = 3),"\n",
                 "p_TF_to_LR=",signif(tmp.res.TFtoLR.pvalue,digits = 3)))+
  scale_color_manual(values = LR.color,
                     name = NULL,
                     labels = c(tmp.LR,y.TF))+
  guides(colour = guide_legend(override.aes = list(size=8)))+
  theme_cowplot(font_size = 22)+
  theme(axis.line  = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        legend.position = "bottom",
        legend.justification = "center",
        plot.title = element_text(hjust = 0.5),
        legend.text = element_text(size = 28))

ggsave(filename = myFileName(prefix = paste0("res/fig/figS_revision_R1_4_",tmp.LR,"_",y.TF,"_"),
                             suffix = ".jpg"),
       width = 6,height = 6,dpi = 350)

ggsave(filename = myFileName(prefix = paste0("res/fig/figS_revision_R1_4_",tmp.LR,"_",y.TF,"_"),
                             suffix = ".pdf"),
       width = 6,height = 6)


length(lag(tmp.IS,n = 1))

tttt <- lag(tmp.IS,n = 2)
tttt <- myRemoveNA(tttt)
cor(tttt,y.TF.level,method = "spearman")



###-------3. simulation -------------
library(lmtest)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggsignif)
library(TimeTalk)
library(ggbeeswarm)
source("code/myUtils.R")

###------3.1 test and EDA---------------

### generate white noise  

N <- 100
x <- as.numeric(arima.sim(model = list(order = c(0,0,0)),n = N))
x <- MinMaxScale(x)
y <- as.numeric(arima.sim(model = list(order = c(0,0,0)),n = N))
y <- MinMaxScale(y)

### simulate the time points
?sort
tmp_indx_xx <- sort(sample(1:N,size = 50))
xx <- x[tmp_indx_xx]
xx_trajCond <- tmp_indx_xx/N
tmp_indx_yy <- sort(sample(1:N,size = 50))
yy <- y[tmp_indx_yy]
yy_trajCond <- tmp_indx_yy/N
plot(xx_trajCond,xx,type = "l")
plot(1:N/N,x,type = "l")

### generate white noise  
x <- as.numeric(arima.sim(model = list(order = c(0,0,0)),n = 200))
x <- MinMaxScale(x)
y <- as.numeric(arima.sim(model = list(order = c(0,0,0)),n = 200))
y <- MinMaxScale(y)

.temp_sample_recover <- function(x,y,
                                 numPts=length(x),
                                 tmp.winsz=0.1,
                                 tmp.sample.size = 50){
  ### simulate the time points
  N <- numPts
  tmp_indx_xx <- sort(sample(1:N,size = tmp.sample.size))
  xx <- x[tmp_indx_xx]
  xx_trajCond <- tmp_indx_xx/N
  tmp_indx_yy <- sort(sample(1:N,size = tmp.sample.size))
  yy <- y[tmp_indx_yy]
  yy_trajCond <- tmp_indx_yy/N
  
  ### perform cell aling imputation
  
  xx <- t(as.data.frame(xx))
  yy <- t(as.data.frame(yy))
  xx <- rbind(xx,xx)
  yy <- rbind(yy,yy)
  
  inter.tmp.mat.1 <- cellAlign::interWeights(
    expDataBatch = xx,
    trajCond = xx_trajCond,
    winSz = tmp.winsz,
    numPts = numPts
  )
  inter.tmp.mat.2 <- cellAlign::interWeights(
    expDataBatch = xx,
    trajCond = yy_trajCond,
    winSz = tmp.winsz,
    numPts = numPts) 
  
  inter.tmp.mat.1 <- cellAlign::scaleInterpolate(inter.tmp.mat.1)
  inter.tmp.mat.2 <- cellAlign::scaleInterpolate(inter.tmp.mat.2)
  time <- inter.tmp.mat.1$traj
  
  inter.tmp.mat.1 <- myRemoveNA(inter.tmp.mat.1$scaledData)
  inter.tmp.mat.2 <- myRemoveNA(inter.tmp.mat.2$scaledData)
  
  x_impute <- inter.tmp.mat.1[1,]
  y_impute <- inter.tmp.mat.2[1,]
  
  res.list <- list(x_impute = x_impute,
                   y_impute = y_impute)
  return(res.list)
}

### generate white noise  
N <- 200
x <- as.numeric(arima.sim(model = list(order = c(0,0,0)),n = N))
x <- MinMaxScale(x)
y <- as.numeric(arima.sim(model = list(order = c(0,0,0)),n = N))
y <- MinMaxScale(y)
var(x)/mean(x)

### try the scale the point
N <- 200
x <- as.numeric(arima.sim(model = list(ma = 0.5),n = N))
error <- as.numeric(arima.sim(model = list(order = c(0,0,0)),n = 200))
y <- 0
for(ii in 2:N){
  y[ii] = x[ii-1] + error[ii]
}
x <- MinMaxScale(x)
y <- MinMaxScale(y)

### generate random walk
simulate_ar <- function(n, phi, sigma = .1) {
  y <- rep(0, n)
  
  for (t in seq(2, n)) {
    y[t] <- phi*y[t-1] + rnorm(1, 0, sigma)
  }
  
  y
}
x <- simulate_ar(n = N,phi = 1,sigma = 1)
y <- simulate_ar(n = N,phi = 1,sigma = 1)

res.list <- .temp_sample_recover(x = x,
                                 y = y,
                                 numPts = length(x),
                                 tmp.winsz = 0.1,
                                 tmp.sample.size = 50)
var(res.list$x_impute)/mean(res.list$x_impute)
var(res.list$y_impute)/mean(res.list$y_impute)

plot.ts(cbind(x,res.list$x_impute,y,res.list$y_impute))

res.p.1 <- grangertest(x,y)
res.p.1 <- res.p.1$`Pr(>F)`[2]
res.p.2 <- grangertest(res.list$x_impute,
                       res.list$y_impute)
res.p.2 <- res.p.2$`Pr(>F)`[2]

res.p.3 <- grangertest(y,x)
res.p.3 <- res.p.3$`Pr(>F)`[2]
res.p.4 <- grangertest(y_impute,x_impute)
res.p.4 <- res.p.4$`Pr(>F)`[2]

#plot.ts(cbind(x,x_impute,y,y_impute))
res.df <- data.frame(model_xy = res.p.1,
                     cellalign_xy = res.p.2,
                     model_yx = res.p.3,
                     cellalign_yx = res.p.4,
                     stringsAsFactors = F)
head(res.df)


####------3.2 white noise model---------------------

####------3.2.1 run the model windowsize 0.1---------
res.list.df <- lapply(1:10,FUN = function(ii){
  cat(paste0("simulation round:",ii),sep = "\n")
  res.list <- replicate(n = 1000,expr = {
    ### generate white noise  
    x <- as.numeric(arima.sim(model = list(order = c(0,0,0)),n = 200))
    x <- MinMaxScale(x)
    y <- as.numeric(arima.sim(model = list(order = c(0,0,0)),n = 200))
    y <- MinMaxScale(y)
    
    ### simulate the time points
    tmp.res.list <- .temp_sample_recover(x = x,
                                       y = y,
                                       numPts = length(x),
                                       tmp.winsz = 0.1,
                                       tmp.sample.size = 50)
    x_impute <- tmp.res.list$x_impute
    y_impute <- tmp.res.list$y_impute
    
    plot.ts(cbind(x,x_impute,y,y_impute))
    res.p.1 <- grangertest(x,y)
    res.p.1 <- res.p.1$`Pr(>F)`[2]
    res.p.2 <- grangertest(x_impute,y_impute)
    res.p.2 <- res.p.2$`Pr(>F)`[2]
    
    res.p.3 <- grangertest(y,x)
    res.p.3 <- res.p.3$`Pr(>F)`[2]
    res.p.4 <- grangertest(y_impute,x_impute)
    res.p.4 <- res.p.4$`Pr(>F)`[2]
    
    
    res.df <- data.frame(model_xy = res.p.1,
                         cellalign_xy = res.p.2,
                         model_yx = res.p.3,
                         cellalign_yx = res.p.4,
                         stringsAsFactors = F)
    head(res.df)
    return(res.df)
  },simplify = F)
  
  tmp.data.plot <- Reduce(rbind,res.list)
  
  tmp.data.plot <- tmp.data.plot %>%
    mutate(result_state_xy = "FAIL",
           result_state_yx = "FAIL") %>%
    mutate(result_state_xy = ifelse(model_xy > 0.05 & cellalign_xy > 0.05,"PASS",result_state_xy)) %>%
    mutate(result_state_yx = ifelse(model_yx > 0.05 & cellalign_yx > 0.05,"PASS",result_state_yx))
  tmp.data.plot$simu_round <- ii
  return(tmp.data.plot)
})

head(tmp.data.plot)
####---------3.2.2 visualization--------------
tmp.data.plot.1 <- Reduce(rbind,res.list.df) %>%
  dplyr::group_by(simu_round,result_state_xy) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(ratio = n/sum(n)) %>%
  dplyr::rename(result_state = result_state_xy) %>%
  dplyr::mutate(result_state = paste0(result_state,"_xy"))
head(tmp.data.plot.1)


tmp.data.plot.2 <- Reduce(rbind,res.list.df) %>%
  dplyr::group_by(simu_round,result_state_yx) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(ratio = n/sum(n)) %>%
  dplyr::rename(result_state = result_state_yx) %>%
  dplyr::mutate(result_state = paste0(result_state,"_yx"))

tmp.data.plot <- rbind(tmp.data.plot.1,
                       tmp.data.plot.2) %>%
  mutate(result_state = factor(result_state,
                               levels = c("FAIL_xy","PASS_xy",
                                          "FAIL_yx","PASS_yx"))) 


x <- tmp.data.plot %>%
  filter(result_state == "FAIL_xy") %>%
  pull(ratio)
y <- tmp.data.plot %>%
  filter(result_state == "PASS_xy") %>%
  pull(ratio)
tmp.pvalue_xy <- wilcox.test(x,y,alternative = "less")$p.value


x <- tmp.data.plot %>%
  filter(result_state == "FAIL_yx") %>%
  pull(ratio)
y <- tmp.data.plot %>%
  filter(result_state == "PASS_yx") %>%
  pull(ratio)
tmp.pvalue_yx <- wilcox.test(x,y,alternative = "less")$p.value

tmp.data.plot


ggplot(tmp.data.plot,
       aes(result_state,ratio,color = result_state))+
  geom_quasirandom()+
  geom_signif(comparisons = list(c("FAIL_xy","PASS_xy"),
                                 c("FAIL_yx","PASS_yx")),
              annotations = c(paste0("p=",signif(tmp.pvalue_xy,4)),
                              paste0("p=",signif(tmp.pvalue_yx,4))),
              y_position = c(0.9),
              size = 1,color = "black",
              textsize = 8,
              vjust = -0.5)+
  xlab(NULL)+
  ggtitle("winsz=0.1")+
  scale_color_manual(values  = c("#33658A","#F6AE2D","#2F4858","#F26419"))+
  scale_y_continuous(breaks = seq(0,1,by=0.2),
                     limits = c(0,1))+
  theme_cowplot(font_size = 22)+
  theme(legend.position = "none")

ggsave(filename = myFileName(prefix = "res/fig/FigS_revision_model_without_grangercausal_simu_test",
                             suffix = ".png"),
       width = 6,height = 6,dpi = 400)


####------3.2.3 run the model windowsize 0.01---------


res.list.df <- lapply(1:10,FUN = function(ii){
  cat(paste0("simulation round:",ii),sep = "\n")
  res.list <- replicate(n = 1000,expr = {
    ### generate white noise  
    x <- as.numeric(arima.sim(model = list(order = c(0,0,0)),n = 200))
    x <- MinMaxScale(x)
    y <- as.numeric(arima.sim(model = list(order = c(0,0,0)),n = 200))
    y <- MinMaxScale(y)
    
    ### simulate the time points
    N <- 200
    tmp_indx_xx <- sample(1:N,size = 50)
    xx <- x[tmp_indx_xx]
    xx_trajCond <- tmp_indx_xx/N
    tmp_indx_yy <- sample(1:N,size = 50)
    yy <- y[tmp_indx_yy]
    yy_trajCond <- tmp_indx_yy/N
    
    ### perform cell aling imputation
    tmp.winsz = 0.01
    numPts = 200
    
    xx <- t(as.data.frame(xx))
    yy <- t(as.data.frame(yy))
    xx <- rbind(xx,xx)
    yy <- rbind(yy,yy)
    
    inter.tmp.mat.1 <- cellAlign::interWeights(
      expDataBatch = xx,
      trajCond = xx_trajCond,
      winSz = tmp.winsz,
      numPts = numPts
    )
    inter.tmp.mat.2 <- cellAlign::interWeights(
      expDataBatch = xx,
      trajCond = yy_trajCond,
      winSz = tmp.winsz,
      numPts = numPts) 
    
    inter.tmp.mat.1 <- cellAlign::scaleInterpolate(inter.tmp.mat.1)
    inter.tmp.mat.2 <- cellAlign::scaleInterpolate(inter.tmp.mat.2)
    time <- inter.tmp.mat.1$traj
    
    inter.tmp.mat.1 <- myRemoveNA(inter.tmp.mat.1$scaledData)
    inter.tmp.mat.2 <- myRemoveNA(inter.tmp.mat.2$scaledData)
    
    x_impute <- inter.tmp.mat.1[1,]
    y_impute <- inter.tmp.mat.2[1,]
    
    res.p.1 <- grangertest(x,y)
    res.p.1 <- res.p.1$`Pr(>F)`[2]
    res.p.2 <- grangertest(x_impute,y_impute)
    res.p.2 <- res.p.2$`Pr(>F)`[2]
    
    res.p.3 <- grangertest(y,x)
    res.p.3 <- res.p.3$`Pr(>F)`[2]
    res.p.4 <- grangertest(y_impute,x_impute)
    res.p.4 <- res.p.4$`Pr(>F)`[2]
    
    
    res.df <- data.frame(model_xy = res.p.1,
                         cellalign_xy = res.p.2,
                         model_yx = res.p.3,
                         cellalign_yx = res.p.4,
                         stringsAsFactors = F)
    #head(res.df)
    return(res.df)
  },simplify = F)
  
  tmp.data.plot <- Reduce(rbind,res.list)
  
  tmp.data.plot <- tmp.data.plot %>%
    mutate(result_state_xy = "FAIL",
           result_state_yx = "FAIL") %>%
    mutate(result_state_xy = ifelse(model_xy > 0.05 & cellalign_xy > 0.05,"PASS",result_state_xy)) %>%
    mutate(result_state_yx = ifelse(model_yx > 0.05 & cellalign_yx > 0.05,"PASS",result_state_yx))
  tmp.data.plot$simu_round <- ii
  return(tmp.data.plot)
})

head(tmp.data.plot)
####---------3.2.2 visualization--------------
tmp.data.plot.1 <- Reduce(rbind,res.list.df) %>%
  dplyr::group_by(simu_round,result_state_xy) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(ratio = n/sum(n)) %>%
  dplyr::rename(result_state = result_state_xy) %>%
  dplyr::mutate(result_state = paste0(result_state,"_xy"))
head(tmp.data.plot.1)


tmp.data.plot.2 <- Reduce(rbind,res.list.df) %>%
  dplyr::group_by(simu_round,result_state_yx) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(ratio = n/sum(n)) %>%
  dplyr::rename(result_state = result_state_yx) %>%
  dplyr::mutate(result_state = paste0(result_state,"_yx"))

tmp.data.plot <- rbind(tmp.data.plot.1,
                       tmp.data.plot.2) %>%
  mutate(result_state = factor(result_state,
                               levels = c("FAIL_xy","PASS_xy",
                                          "FAIL_yx","PASS_yx"))) 


x <- tmp.data.plot %>%
  filter(result_state == "FAIL_xy") %>%
  pull(ratio)
y <- tmp.data.plot %>%
  filter(result_state == "PASS_xy") %>%
  pull(ratio)
tmp.pvalue_xy <- wilcox.test(x,y,alternative = "less")$p.value


x <- tmp.data.plot %>%
  filter(result_state == "FAIL_yx") %>%
  pull(ratio)
y <- tmp.data.plot %>%
  filter(result_state == "PASS_yx") %>%
  pull(ratio)
tmp.pvalue_yx <- wilcox.test(x,y,alternative = "less")$p.value

tmp.data.plot


ggplot(tmp.data.plot,
       aes(result_state,ratio,color = result_state))+
  geom_quasirandom()+
  geom_signif(comparisons = list(c("FAIL_xy","PASS_xy"),
                                 c("FAIL_yx","PASS_yx")),
              annotations = c(paste0("p=",signif(tmp.pvalue_xy,4)),
                              paste0("p=",signif(tmp.pvalue_yx,4))),
              y_position = c(0.9),
              size = 1,color = "black",
              textsize = 8,
              vjust = -0.5)+
  xlab(NULL)+
  ggtitle("winsz=0.1")+
  scale_color_manual(values  = c("#33658A","#F6AE2D","#2F4858","#F26419"))+
  scale_y_continuous(breaks = seq(0,1,by=0.2),
                     limits = c(0,1))+
  theme_cowplot(font_size = 22)+
  theme(legend.position = "none")

ggsave(filename = myFileName(prefix = "res/fig/FigS_revision_model_without_grangercausal_simu_test",
                             suffix = ".png"),
       width = 6,height = 6,dpi = 400)


####-------3.3 try another model ---------------------------------

####-------3.3.1 run the simulation model-----------------

res.list.df <- lapply(1:10,FUN = function(ii){
  cat(paste0("simulation round:",ii),sep = "\n")
  res.list <- replicate(n = 1000,expr = {
    
    ### try the scale the point
    N <- 200
    x <- as.numeric(arima.sim(model = list(ma = 0.5),n = N))
    error <- as.numeric(arima.sim(model = list(order = c(0,0,0)),n = 200))
    y <- 0
    for(ii in 2:N){
      y[ii] = x[ii-1] + error[ii]
    }
    x <- MinMaxScale(x)
    y <- MinMaxScale(y)
    
    
    ### simulate the time points
    tmp_indx_xx <- sample(1:N,size = 50)
    xx <- x[tmp_indx_xx]
    xx_trajCond <- tmp_indx_xx/N
    tmp_indx_yy <- sample(1:N,size = 50)
    yy <- y[tmp_indx_yy]
    yy_trajCond <- tmp_indx_yy/N
    
    ### perform cell aling imputation
    tmp.winsz = 0.1
    numPts = 200
    
    xx <- t(as.data.frame(xx))
    yy <- t(as.data.frame(yy))
    xx <- rbind(xx,xx)
    yy <- rbind(yy,yy)
    
    
    inter.tmp.mat.1 <- cellAlign::interWeights(
      expDataBatch = xx,
      trajCond = xx_trajCond,
      winSz = tmp.winsz,
      numPts = numPts
    )
    inter.tmp.mat.2 <- cellAlign::interWeights(
      expDataBatch = xx,
      trajCond = yy_trajCond,
      winSz = tmp.winsz,
      numPts = numPts) 
    
    inter.tmp.mat.1 <- cellAlign::scaleInterpolate(inter.tmp.mat.1)
    inter.tmp.mat.2 <- cellAlign::scaleInterpolate(inter.tmp.mat.2)
    time <- inter.tmp.mat.1$traj
    
    inter.tmp.mat.1 <- myRemoveNA(inter.tmp.mat.1$scaledData)
    inter.tmp.mat.2 <- myRemoveNA(inter.tmp.mat.2$scaledData)
    
    x_impute <- inter.tmp.mat.1[1,]
    y_impute <- inter.tmp.mat.2[1,]
    
    #plot.ts(cbind(x,x_impute,y,y_impute))
    #try the granger test
    res.p.1 <- grangertest(x,y)
    res.p.1 <- res.p.1$`Pr(>F)`[2]
    res.p.2 <- grangertest(x_impute,y_impute)
    res.p.2 <- res.p.2$`Pr(>F)`[2]
    
    res.p.3 <- grangertest(y,x)
    res.p.3 <- res.p.3$`Pr(>F)`[2]
    res.p.4 <- grangertest(y_impute,x_impute)
    res.p.4 <- res.p.4$`Pr(>F)`[2]
    
    
    res.df <- data.frame(model_xy = res.p.1,
                         cellalign_xy = res.p.2,
                         model_yx = res.p.3,
                         cellalign_yx = res.p.4,
                         stringsAsFactors = F)
    #head(res.df)
    return(res.df)
  },simplify = F)
  
  tmp.data.plot <- Reduce(rbind,res.list)
  
  tmp.data.plot <- tmp.data.plot %>%
    mutate(result_state_xy = "FAIL",
           result_state_yx = "FAIL") %>%
    mutate(result_state_xy = ifelse(model_xy < 0.05 & cellalign_xy < 0.05,"PASS",result_state_xy)) %>%
    mutate(result_state_yx = ifelse(model_yx < 0.05 & cellalign_yx < 0.05,"PASS",result_state_yx))
  tmp.data.plot$simu_round <- ii
  return(tmp.data.plot)
})

#####-----------3.3.2 visualization results----------------
res.df <- Reduce(rbind,res.list.df)
head(res.df)
head(res.list.df[[1]])

tmp.data.plot.1 <- Reduce(rbind,res.list.df) %>%
  dplyr::group_by(simu_round,result_state_xy) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(ratio = n/sum(n)) %>%
  dplyr::rename(result_state = result_state_xy) %>%
  dplyr::mutate(result_state = paste0(result_state,"_xy"))
head(tmp.data.plot.1)


tmp.data.plot.2 <- Reduce(rbind,res.list.df) %>%
  dplyr::group_by(simu_round,result_state_yx) %>%
  dplyr::summarise(n = n()) %>%
  dplyr::mutate(ratio = n/sum(n)) %>%
  dplyr::rename(result_state = result_state_yx) %>%
  dplyr::mutate(result_state = paste0(result_state,"_yx"))

head(tmp.data.plot.1)
head(tmp.data.plot.2)

tmp.data.plot <- rbind(tmp.data.plot.1,
                       tmp.data.plot.2) %>%
  mutate(result_state = factor(result_state,
                               levels = c("FAIL_xy","PASS_xy",
                                          "FAIL_yx","PASS_yx"))) 

head(tmp.data.plot)


x <- tmp.data.plot %>%
  filter(result_state == "FAIL_xy") %>%
  pull(ratio)
y <- tmp.data.plot %>%
  filter(result_state == "PASS_xy") %>%
  pull(ratio)
tmp.pvalue_xy <- wilcox.test(x,y,alternative = "less")$p.value


x <- tmp.data.plot %>%
  filter(result_state == "FAIL_yx") %>%
  pull(ratio)
y <- tmp.data.plot %>%
  filter(result_state == "PASS_yx") %>%
  pull(ratio)
tmp.pvalue_yx <- wilcox.test(x,y,alternative = "less")$p.value

ggplot(tmp.data.plot,
       aes(result_state,ratio,color = result_state))+
  geom_quasirandom()+
  geom_signif(comparisons = list(c("FAIL_xy","PASS_xy"),
                                 c("FAIL_yx","PASS_yx")),
              annotations = c(paste0("p=",signif(tmp.pvalue_xy,4)),
                              paste0("p=",signif(tmp.pvalue_yx,4))),
              y_position = c(0.9),
              size = 1,color = "black",
              textsize = 8,
              vjust = -0.5)+
  scale_color_manual(values  = c("#33658A","#F6AE2D","#2F4858","#F26419"))+
  scale_y_continuous(breaks = seq(0,1,by=0.2),
                     limits = c(0,1))+
  theme_cowplot(font_size = 22)+
  theme(legend.position = "none")+
  xlab(NULL)

ggsave(filename = myFileName(prefix = "res/fig/FigS_revision_model_with_grangercausal_simu_test",
                             suffix = ".png"),
       width = 6,height = 6,dpi = 400)
