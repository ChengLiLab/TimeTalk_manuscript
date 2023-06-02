####TimeTalk manuscript
####Author:wlt
####Fig5
####This figure will extend TimeTalk to blastoid datasets
####This aims to extend TimeTalk to show its advance in solving problem compare to cellchat and cellcall;
####start at: 20220411
####break through2: 20220425
####Rewrite the code started in 20230201


for(ii in 1:100000){
  cat(ii,sep = "\n")
  Sys.sleep(time = ii)
}

####reset
rm(list = ls())
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
source(file = "code/myUtils.R")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")


####-------1.1 load B-blastoid dataset ------------

####EPS blastoid
tmp.dir <- "data/blastoid/GSE135701/EPS_blastoid/"
tmp.mat <- Read10X(data.dir = tmp.dir)
#br.out <- DropletUtils::barcodeRanks(tmp.mat)

# # Making a plot test the barcode rank
# plot(br.out$rank, 
#      br.out$total, 
#      log="xy", 
#      xlab="Rank", 
#      ylab="Total")
# o <- order(br.out$rank)
# lines(br.out$rank[o], br.out$fitted[o], col="red")
# abline(h=metadata(br.out)$knee, col="dodgerblue", lty=2)
# abline(h=metadata(br.out)$inflection, col="forestgreen", lty=2)
# legend("bottomleft", lty=2, col=c("dodgerblue", "forestgreen"), 
#        legend=c("knee", "inflection"))
# range(br.out$total[!is.na(br.out$fitted)])

tttt <- DropletUtils::emptyDrops(tmp.mat,BPPARAM=SnowParam(workers = 10))
tmp.cell.id <- as.data.frame(tttt) %>%
  dplyr::filter(FDR < 0.01) %>%
  rownames()
tmp.mat.1 <- tmp.mat[,tmp.cell.id]
seu <- CreateSeuratObject(tmp.mat.1,
                          min.features = 0,
                          min.cells = 0,
                          project = "EPS_blastoid")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
VlnPlot(seu,features = c("nCount_RNA","nFeature_RNA","percent.mt"))
seu <- subset(seu, subset = nFeature_RNA > 2000 & percent.mt < 20)
VlnPlot(seu,features = c("nCount_RNA","nFeature_RNA","percent.mt"))
seu.1 <- seu

####blastocyst
tmp.dir <- "data/blastoid/GSE135701/Blastocyst/"
tmp.mat <- Read10X(data.dir = tmp.dir)
tttt <- DropletUtils::emptyDrops(tmp.mat,BPPARAM=SnowParam(workers = 10))
tmp.cell.id <- as.data.frame(tttt) %>%
  dplyr::filter(FDR < 0.01) %>%
  rownames()
tmp.mat.1 <- tmp.mat[,tmp.cell.id]
seu <- CreateSeuratObject(tmp.mat.1,
                          min.features = 0,
                          min.cells = 0,
                          project = "blastocyst")
#colnames(seu@meta.data)
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")
VlnPlot(seu,features = c("nCount_RNA",
                         "nFeature_RNA",
                         "percent.mt"))
seu <- subset(seu, subset = nFeature_RNA > 2000 & percent.mt < 20)
VlnPlot(seu,features = c("nCount_RNA","nFeature_RNA","percent.mt"))
seu.2 <- seu

seu <- seu.2
tmp.cell.id <- readLines("res/txt/B_blastoid_blastocyst_to_remove.txt")
tmp.cell.id <- WhichCells(seu,cells = tmp.cell.id,invert = T)
seu.2 <- subset(seu.2,cells=tmp.cell.id)


###-----1.2 get additional blastocyst dataset----------


####-----1.2.1 get data----------

### download matrix 
tmp.id <- "GSE84892"
getGEOSuppFiles(GEO = tmp.id,baseDir = "data",makeDirectory = F)

### download metadata
test <- GEOquery::getGEO(tmp.id)
tmp.meta <- test$GSE84892_series_matrix.txt.gz@phenoData@data
write.table(tmp.meta,
            file = paste0("data/",tmp.id,"_metadata.txt"),
            sep = "\t",
            quote = F,
            row.names = F,
            col.names = T)

####----1.2.2 convert to seurat----------

# ### I hvae try it, but for simlicity, this part was deprecated
# tmp.mat <- read_delim(file = "data/GSE84892_counts.txt.gz",delim = "\t") %>%
#   column_to_rownames("X1")
# tmp.meta <- read_delim(file = "data/GSE84892_metadata.txt",delim = "\t")
# tmp.meta <- tmp.meta %>%
#   dplyr::select(title,`lineage:ch1`,`Stage:ch1`) %>%
#   dplyr::rename(cell_id=title,
#                 lineage=`lineage:ch1`,
#                 stage=`Stage:ch1`) %>%
#   dplyr::mutate(tmp.rownames=cell_id) %>%
#   column_to_rownames("tmp.rownames")
# ### filter by meta data
# tmp.mat <- tmp.mat[,tmp.meta$cell_id]
# seu <- CreateSeuratObject(tmp.mat,
#                           min.features = 0,
#                           min.cells = 0,
#                           meta.data = tmp.meta,
#                           project = "blastocyst_2")
# seu[["percent.mt"]] <- PercentageFeatureSet(seu,pattern = "^mt-")
# seu.3 <- seu


####------1.3 batch effect correct by CCA pipeline -------------

####------1.3.1 orginial pipeline----------
tmp.list <- list(seu.1,seu.2)
names(tmp.list) <- c("EPS_blastoid","blastocyst")
saveRDS(object = tmp.list,file = myFileName(prefix = "res/R/B_blastoid_seurat_list",suffix = ".rds"))
tmp.list <- lapply(names(tmp.list),FUN = function(x){
  tmp.list[[x]] <- NormalizeData(tmp.list[[x]], verbose = FALSE)
  tmp.list[[x]] <- FindVariableFeatures(tmp.list[[x]], 
                                        selection.method = "vst", 
                                        nfeatures = 2000,
                                        verbose = FALSE)
  return(tmp.list[[x]])
})

#?FindIntegrationAnchors
tmp.list <- FindIntegrationAnchors(object.list = tmp.list, 
                                   dims = 1:30,
                                   anchor.features = 2000,
                                   k.anchor = 5,
                                   k.filter = 30)
#?IntegrateData
seu.integrated <- IntegrateData(anchorset = tmp.list, dims = 1:30)
###RunSeurat pipeline
DefaultAssay(seu.integrated) <- "integrated"
seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
seu.integrated <- RunPCA(seu.integrated, npcs = 30, verbose = FALSE)
seu.integrated <- RunUMAP(seu.integrated, 
                          reduction = "pca",
                          dims = 1:30,
                          umap.method = "umap-learn",
                          seed.use = 42L)
seu.integrated <- FindNeighbors(seu.integrated,dims = 1:30)
seu.integrated <- FindClusters(seu.integrated,
                               resolution = 0.3)
DefaultAssay(seu.integrated) <- "RNA"


# ####-------1.3.2 test approach-----------------
# 
# tmp.list <- readRDS(file = "res/R/B_blastoid_seurat_list_2022042116.rds")
# tmp.list <- c(tmp.list,list(seu.3))
# names(tmp.list) <-  c("EPS_blastoid","blastocyst_1","blastocyst_2")
# tmp.list <- lapply(names(tmp.list),FUN = function(x){
#   tmp.list[[x]] <- NormalizeData(tmp.list[[x]], verbose = FALSE)
#   tmp.list[[x]] <- FindVariableFeatures(tmp.list[[x]], 
#                                         selection.method = "vst", 
#                                         nfeatures = 3000,
#                                         verbose = FALSE)
#   return(tmp.list[[x]])
# })
# 
# #?FindIntegrationAnchors
# tmp.list <- FindIntegrationAnchors(object.list = tmp.list, 
#                                    dims = 1:20,
#                                    anchor.features = 3000,
#                                    k.anchor = 5,
#                                    k.filter = 30)
# #?IntegrateData
# seu.integrated <- IntegrateData(anchorset = tmp.list, dims = 1:20)
# ###RunSeurat pipeline
# DefaultAssay(seu.integrated) <- "integrated"
# seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
# seu.integrated <- RunPCA(seu.integrated, npcs = 30, verbose = FALSE)
# seu.integrated <- RunUMAP(seu.integrated, 
#                           reduction = "pca",
#                           dims = 1:30,
#                           umap.method = "umap-learn",
#                           seed.use = 42L)
# seu.integrated <- FindNeighbors(seu.integrated,dims = 1:30)
# seu.integrated <- FindClusters(seu.integrated,
#                                resolution = 0.3)
# DefaultAssay(seu.integrated) <- "RNA"
# 
# saveRDS(seu.integrated,file = myFileName("res/R/B_blastoid_seurat_inter_3_data",suffix = ".rds"))
####-------1.4 assign clusters-------------

####-------1.4.1 plot markers----------

TE.markers <- c("Cdx2","Krt8","Krt18","Ascl2","Tacstd2")
ICM.markers <- c("Pou5f1","Nanog","Sox2","Pramel5","Crxos")
PE.markers <- c("Gata4","Gata6","Sox17","Pdgfra","Col4a1")
p1 <- UMAPPlot(seu.integrated,label=T,pt.size=1.5)
p1
tmp.markers <- c(TE.markers,
                 "Gapdh",ICM.markers,
                 "Ppia",PE.markers)
# tmp.seu <- seu.integrated
# tmp.seu.list <- SplitObject(tmp.seu,split.by = "orig.ident")
# seu.blastoid <- tmp.seu.list$EPS_blastoid
# seu.blastocyst <- tmp.seu.list$blastocyst
# 
# 
# FeaturePlot(seu.integrated,features = c("Cdx2","Krt8","Tead4"),split.by = "orig.ident")


# tmp.cell.id <- WhichCells(seu.blastocyst,
#                           invert = F,
#                           expression = UMAP_1 > 10 & UMAP_2 > 10)
# cat(tmp.cell.id,
#     file = "res/txt/B_blastoid_blastocyst_to_remove.txt",
#     sep = "\n")
# 
# tmp.cell.id <- readLines("res/txt/B_blastoid_blastocyst_to_remove.txt")



# p2 <- FeaturePlot(seu.integrated,features = c(TE.markers,
#                                    "Gapdh",ICM.markers,
#                                    "Ppia",PE.markers),combine = F) 

p2 <- lapply(tmp.markers, function(x){
  res.p <- FeaturePlot(seu.integrated,features = x)+
    scale_color_gradientn(colours = rdwhbu(256))
  return(res.p)
})

p <- wrap_plots(c(list(p1),p2),nrow = 3) & 
  theme_cowplot(font_size = 18) &
  NoAxes() &
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5,
                                  face = "bold.italic"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
ggsave(plot = p,
       filename = myFileName(prefix = "res/fig/Fig5_EPS_blastoid_markers_raw",suffix = ".png"),
       dpi = 350,width = 16,height = 6)

UMAPPlot(seu.integrated,split.by="orig.ident")


####-----1.4.2 assign cluster-------
levels(seu.integrated)

head(seu.integrated[[]])
seu.integrated <- RenameIdents(seu.integrated,
                               "0"="inter1",
                               "1"="PE",
                               "2"="inter3",
                               "3"="EPI",
                               "4"="inter2",
                               "5"="TE",
                               "6"="ICM")
levels(seu.integrated) <- c("ICM","EPI","PE","TE","inter1","inter2","inter3")

tmp.colors <- c(divergentcolor(7)[1:4],
                grey.colors(3,start = 0.5,rev = T))
scales::show_col(tmp.colors)

p1 <- UMAPPlot(seu.integrated) + 
  scale_color_manual(values = tmp.colors)+
  ggtitle("Cell Type")
p1
p <- wrap_plots(c(list(p1),p2),nrow = 3) & 
  theme_cowplot(font_size = 18) &
  NoAxes() &
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5,
                                  face = "bold.italic"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
p
ggsave(plot = p,
       filename = myFileName(prefix = "res/fig/Fig5_EPS_blastoid_markers",
                                      suffix = ".png"),
       dpi = 350,width = 18,height = 8)


seu.integrated$CellType <- Idents(seu.integrated)
saveRDS(seu.integrated,
        file = myFileName(prefix = "res/R/EPS_blastoid",
                          suffix = ".rds"))


####-------------1.5 show markers-------------------

####In this section, we need to detail maake it clear the cell types we have here

###--------1.5.1 tSNE plot -----------------
####load data
seu.integrated <- readRDS(file = "res/R/EPS_blastoid_2022042002.rds")
tmp.seu <- seu.integrated
DefaultAssay(tmp.seu) <- "RNA"
tmp.seu <- RunTSNE(tmp.seu,
                   dims = 1:30,
                   reduction = "pca",
                   seed.use = 42)

### tSNE plot by cell type
tmp.colors <- c(divergentcolor(7)[1:4],
                grey.colors(3,start = 0.5,rev = T))

DimPlot(tmp.seu,
        reduction = "tsne",
        pt.size = 1.5,
        group.by = "CellType",
        split.by = "orig.ident")+
  ggtitle(NULL)+
  scale_color_manual(values = tmp.colors)+
  theme_cowplot(font_size = 28) +
  NoAxes() +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = "pt"),
        aspect.ratio = 1,
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))

ggsave(filename = myFileName(prefix = "res/fig/fig5_B_blastoids_tSNE",suffix = ".png"),
       width = 8,height = 4,dpi = 400)
ggsave(filename = myFileName(prefix = "res/fig/fig5_B_blastoids_tSNE",suffix = ".pdf"),
       width = 8,height = 4,dpi = 400)

### tSNE plot by orig.ident
scales::show_col(blues(3))

DimPlot(tmp.seu,
        reduction = "tsne",
        pt.size = 1.5,
        group.by = "orig.ident")+
  ggtitle(label = NULL) +
  scale_color_manual(values = c("#2369ad","#b9d9e9"))+
  theme_cowplot(font_size = 28) +
  NoAxes() +
  theme(legend.position = "top",
        legend.justification = "center",
        aspect.ratio = 1,
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
ggsave(filename = myFileName(prefix = "res/fig/fig5_B_blastoids_tSNE_by_source",suffix = ".png"),
       width = 6,height = 4,dpi = 400)
ggsave(filename = myFileName(prefix = "res/fig/fig5_B_blastoids_tSNE_by_source",suffix = ".pdf"),
       width = 6,height = 4,dpi = 400)


###--------1.5.2 UMAP plot -----------------

### UMAP plot by cell type
tmp.colors <- c(divergentcolor(7)[1:4],
                grey.colors(3,start = 0.5,rev = T))

DimPlot(tmp.seu,
        reduction = "umap",
        pt.size = 1.5,
        group.by = "CellType",
        split.by = "orig.ident")+
  ggtitle(NULL)+
  scale_color_manual(values = tmp.colors)+
  theme_cowplot(font_size = 28) +
  NoAxes() +
  theme(plot.margin = margin(0.5,0.5,0.5,0.5,unit = "pt"),
        aspect.ratio = 1,
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))

ggsave(filename = myFileName(prefix = "res/fig/fig5_b_blastoids_UMAP",suffix = ".png"),
       width = 8,height = 4,dpi = 400)
ggsave(filename = myFileName(prefix = "res/fig/fig5_B_blastoids_UMAP",suffix = ".pdf"),
       width = 8,height = 4,dpi = 400)

### UMAP plot by orig.ident
scales::show_col(blues(3))

DimPlot(tmp.seu,
        reduction = "umap",
        pt.size = 1.5,
        group.by = "orig.ident")+
  ggtitle(label = NULL) +
  scale_color_manual(values = c("#2369ad","#b9d9e9"))+
  theme_cowplot(font_size = 28) +
  NoAxes() +
  theme(legend.position = "top",
        legend.justification = "center",
        aspect.ratio = 1,
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
ggsave(filename = myFileName(prefix = "res/fig/fig5_B_blastoids_UMAP_by_source",suffix = ".png"),
       width = 6,height = 4,dpi = 400)
ggsave(filename = myFileName(prefix = "res/fig/fig5_B_blastoids_UMAP_by_source",suffix = ".pdf"),
       width = 6,height = 4,dpi = 400)


####----------1.5.3 Featureplot -------------------

tmp.colors <- c(divergentcolor(7)[1:4],
                grey.colors(3,start = 0.5,rev = T))
scales::show_col(tmp.colors)

####UMAP plot

p1 <- UMAPPlot(tmp.seu) + 
  scale_color_manual(values = tmp.colors)+
  ggtitle("Cell Type")
p1

####Featureplot
TE.markers <- c("Cdx2","Krt8","Krt18","Ascl2","Tacstd2")
ICM.markers <- c("Pou5f1","Nanog","Sox2","Pramel5","Crxos")
PE.markers <- c("Gata4","Gata6","Sox17","Pdgfra","Col4a1")


p2 <- lapply(tmp.markers, function(x){
  res.p <- FeaturePlot(tmp.seu,features = x)+
    scale_color_gradientn(colours = rdwhbu(256))
  return(res.p)
})


p <- wrap_plots(c(list(p1),p2),nrow = 3) & 
  theme_cowplot(font_size = 18) &
  NoAxes() &
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5,
                                  face = "bold.italic"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
ggsave(plot = p,
       filename = myFileName(prefix = "res/fig/Fig5_EPS_blastoid_markers",
                             suffix = ".png"),
       dpi = 350,width = 18,height = 8)


####--------1.5.4 Detail investigation the location of ICM and mourla---------

seu <- readRDS(file = "res/R/EPS_blastoid_2022042002.rds")
table(seu$orig.ident)
FeaturePlot(seu,
            features = c("Zscan4c","Zscan4d","Zscan4f"),
            order = T,
            keep.scale = "all",
            split.by = "orig.ident") & 
  NoAxes() &
  theme(legend.position = "right",
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
ggsave(filename = myFileName(prefix = "res/fig/Fig5_check_the_2C_like_feature_expression",suffix = ".png"),
       width = 6,height = 6,dpi = 350)


features_2C <- list(c("Zscan4c","Zscan4d","Zscan4f"))
seu <- AddModuleScore(
  object = seu,
  features = features_2C,
  name = 'features_2C_score'
)
p <- FeaturePlot(seu,pt.size = 1.5,features = "features_2C_score1",order = T) + 
  scale_color_gradientn(colours = rdwhbu(256)) +
  theme_cowplot(font_size = 28)+
  NoAxes() +
  ggtitle("2C score")+
  theme(legend.position = "right",
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
p
myggsave(p,prefix = "res/fig/fig5_B_balstoid_2C_like_features",
         suffix = ".png",width = 6,height = 6,dpi=400)
myggsave(p,prefix = "res/fig/fig5_B_balstoid_2C_like_features",
         suffix = ".pdf",width = 6,height = 6,dpi=400)


tmp.cell.id <- WhichCells(object = seu,idents = "ICM",expression = features_2C_score1>0)
seu <- SetIdent(object = seu,cells = tmp.cell.id,value = "2C_like")
levels(seu) <- c("2C_like","ICM","EPI","PE","TE","inter1","inter2","inter3")
tmp.colors <- c("#F2BE58",divergentcolor(7)[1:4],
                grey.colors(3,start = 0.5,rev = T))
p <- DimPlot(seu,pt.size = 1.5,reduction = "umap",shuffle = T)+
  theme_cowplot(font_size = 28)+
  scale_color_manual(values = tmp.colors)+
  NoAxes() +
  theme(legend.position = "right",
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
myggsave(p,prefix = "res/fig/fig5_B_balstoid_fine_cluster",
         suffix = ".png",width = 6,height = 6,dpi=400)
myggsave(p,prefix = "res/fig/fig5_B_balstoid_fine_cluster",
         suffix = ".pdf",width = 6,height = 6,dpi=400)

saveRDS(object = seu,
        file = myFileName(prefix = "res/R/res/R/B_blastoid_seurat_inter_2_data",suffix = ".rds"))



FeaturePlot(seu,
            features = c("T"),
            order = T,
            keep.scale = "all",
            split.by = "orig.ident") & 
  NoAxes() &
  theme(legend.position = "right",
        aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
ggsave(filename = myFileName(prefix = "res/fig/Fig5_check_the_intermediate_feature_expression",
                             suffix = ".jpg"),
       width = 6,height = 3,dpi = 350)




####----------2. TimeTalk ------------

# ####----------2.1 calculate pseudotime------------
# 
# ####----------2.1.1 load data----------
# seu.integrated <- readRDS(file = "res/R/EPS_blastoid_2022041922.rds")
# tmp.seu <- seu.integrated
# DefaultAssay(tmp.seu) <- "RNA"
# tmp.seu <- RunTSNE(tmp.seu,
#                    dims = 1:30,
#                    reduction = "pca",
#                    seed.use = 42)
# 
# tmp.seu.list <- SplitObject(tmp.seu,split.by = "orig.ident")
# seu.blastoid <- tmp.seu.list$EPS_blastoid
# seu.blastocyst <- tmp.seu.list$blastocyst
# 
# ####-------2.1.2 test monocle2 pipeline-----------
# 
# ####Define function
# tmp.convert <- function (x, assay = NULL, reduction = NULL,
#                          expressionFamily=NULL) {
#   
#   if (!Seurat:::PackageCheck("monocle", error = FALSE)) {
#     stop("Please install monocle from Bioconductor before converting to a CellDataSet object")
#   }
#   else if (packageVersion(pkg = "monocle") >= package_version(x = "2.99.0")) {
#     stop("Seurat can only convert to/from Monocle v2.X objects")
#   }
#   assay <- assay %||% DefaultAssay(object = x)
#   counts <- GetAssayData(object = x, assay = assay, slot = "data")
#   cell.metadata <- x[[]]
#   feature.metadata <- x[[assay]][[]]
#   if (!"gene_short_name" %in% colnames(x = feature.metadata)) {
#     feature.metadata$gene_short_name <- rownames(x = feature.metadata)
#   }
#   pd <- new(Class = "AnnotatedDataFrame", data = cell.metadata)
#   fd <- new(Class = "AnnotatedDataFrame", data = feature.metadata)
#   cds <- monocle::newCellDataSet(cellData = counts, phenoData = pd, 
#                                  featureData = fd, expressionFamily = expressionFamily)
#   if ("monocle" %in% names(x = Misc(object = x))) {
#     monocle::cellPairwiseDistances(cds = cds) <- Misc(object = x, 
#                                                       slot = "monocle")[["cellPairwiseDistances"]]
#     monocle::minSpanningTree(cds = cds) <- Misc(object = x, 
#                                                 slot = "monocle")[["minSpanningTree"]]
#     Biobase::experimentData(cds = cds) <- Misc(object = x, 
#                                                slot = "monocle")[["experimentData"]]
#     Biobase::protocolData(cds = cds) <- Misc(object = x, 
#                                              slot = "monocle")[["protocolData"]]
#     Biobase::classVersion(cds = cds) <- Misc(object = x, 
#                                              slot = "monocle")[["classVersion"]]
#     slot(object = cds, name = "lowerDetectionLimit") <- Misc(object = x, 
#                                                              slot = "monocle")[["lowerDetectionLimit"]]
#     slot(object = cds, name = "dispFitInfo") <- Misc(object = x, 
#                                                      slot = "monocle")[["dispFitInfo"]]
#     slot(object = cds, name = "auxOrderingData") <- Misc(object = x, 
#                                                          slot = "monocle")[["auxOrderingData"]]
#     slot(object = cds, name = "auxClusteringData") <- Misc(object = x, 
#                                                            slot = "monocle")[["auxClusteringData"]]
#   }
#   dr.slots <- c("reducedDimS", "reducedDimK", "reducedDimW", 
#                 "reducedDimA")
#   reduction <- reduction %||% Seurat:::DefaultDimReduc(object = x, assay = assay)
#   if (!is.null(x = reduction)) {
#     if (grepl(pattern = "tsne", x = tolower(x = reduction))) {
#       slot(object = cds, name = "dim_reduce_type") <- "tSNE"
#       monocle::reducedDimA(cds = cds) <- t(x = Embeddings(object = x[[reduction]]))
#     }
#     else {
#       slot(object = cds, name = "dim_reduce_type") <- reduction
#       monocle::reducedDimA(cds = cds) <- Loadings(object = x[[reduction]])
#       slot(object = cds, name = "reducedDimS") <- Embeddings(object = x[[reduction]])
#     }
#     for (ii in dr.slots) {
#       if (ii %in% names(x = slot(object = x[[reduction]], 
#                                  name = "misc"))) {
#         slot(object = cds, name = ii) <- slot(object = x[[reduction]], 
#                                               name = "misc")[[ii]]
#       }
#     }
#   }
#   return(cds)
# }
# 
# 
# ####---------2.1.2.1 blastocyst--------------
# 
# tmp.levels <- levels(tmp.seu)
# tmp.seu <- tmp.seu[,WhichCells(tmp.seu,
#                                idents = tmp.levels)]
# tmp.cds <- tmp.convert(tmp.seu,reduction = "tsne",
#                        expressionFamily = negbinomial.size())
# tmp.cds <- estimateSizeFactors(tmp.cds)
# tmp.cds <- clusterCells(tmp.cds)
# tmp.diff <- differentialGeneTest(tmp.cds,
#                                  fullModelFormulaStr = "~CellType",
#                                  cores = 10)
# tmp.ordergene <- tmp.diff %>%
#   dplyr::filter(qval < 0.05) %>%
#   rownames()
# tmp.cds <- setOrderingFilter(tmp.cds,tmp.ordergene)
# tmp.cds <- reduceDimension(cds = tmp.cds,
#                            method="DDRTree",
#                            norm_method = "none")
# tmp.cds <- orderCells(tmp.cds)
# 
# p1 <- plot_cell_trajectory(tmp.cds, 
#                            color_by = "CellType")+
#   ggtitle("blastocyst")+
#   theme(axis.line = element_line(size = 1),
#         plot.title = element_text(face = "bold",
#                                   hjust = 0.5,size = 18))
# p1
# ggsave(filename = myFileName(prefix = "res/fig/fig5_test_timeTalk_1",suffix = ".png"),
#        width = 8,height = 8,dpi = 350)
# 
# ####---------2.1.2.1 blastoid--------------
# 
# tmp.seu <- seu.blastoid
# tmp.levels <- levels(tmp.seu)
# tmp.seu <- tmp.seu[,WhichCells(tmp.seu,
#                                idents = tmp.levels)]
# tmp.cds <- tmp.convert(tmp.seu,reduction = "tsne",
#                        expressionFamily = negbinomial.size())
# tmp.cds <- estimateSizeFactors(tmp.cds)
# expressed_genes <- VariableFeatures(tmp.seu)
# tmp.cds <- clusterCells(tmp.cds)
# tmp.diff <- differentialGeneTest(tmp.cds[expressed_genes,],
#                                  fullModelFormulaStr = "~CellType",
#                                  cores = 10)
# tmp.ordergene <- tmp.diff %>%
#   dplyr::filter(qval < 0.05) %>%
#   rownames()
# tmp.cds <- setOrderingFilter(tmp.cds,tmp.ordergene)
# tmp.cds <- reduceDimension(cds = tmp.cds,
#                            method="DDRTree",
#                            norm_method = "none")
# tmp.cds <- orderCells(tmp.cds)
# 
# p1 <- plot_cell_trajectory(tmp.cds, 
#                            color_by = "CellType")+
#   ggtitle("blastoid")+
#   theme(axis.line = element_line(size = 1),
#         plot.title = element_text(face = "bold",
#                                   hjust = 0.5,size = 18))
# p1
# ggsave(filename = myFileName(prefix = "res/fig/fig5_test_timeTalk_2",suffix = ".png"),
#        width = 8,height = 8,dpi = 350)

####------2.2 use monocle3, confirmed--------------
####20220420
####------2.2.1 load data and preprocess------

seu.integrated <- readRDS(file = "res/R/B_blastoid_seurat_inter_2_data_2022042218.rds")
tmp.seu <- seu.integrated
tmp.seu$CellType <- Idents(tmp.seu)
table(tmp.seu$CellType)


mat <- GetAssayData(tmp.seu,assay = "RNA",slot = "counts")
tmp.seu.pdata <- tmp.seu[[]]
tmp.seu.fdata <- data.frame(gene_short_name = rownames(mat))
rownames(tmp.seu.fdata) <- rownames(mat)
cds <- new_cell_data_set(expression_data = mat,
                         cell_metadata = tmp.seu.pdata,
                         gene_metadata = tmp.seu.fdata)

####----2.2.2 run monocle3 pipeline ----------
cds <- preprocess_cds(cds = cds,num_dim = 50)
cds <- reduce_dimension(cds,reduction_method= "UMAP")
p1 <- plot_cells(cds,group_label_size = 6,
           reduction_method = "UMAP",
           color_cells_by = "orig.ident")+
  scale_color_manual(values = c("#2369ad","#b9d9e9"))+
  theme_cowplot(font_size = 28) +
  NoAxes() +
  theme(legend.position = "none",
        legend.justification = "center",
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

p2 <- plot_cells(cds,group_label_size = 6,
                 reduction_method = "UMAP",
                 color_cells_by = "CellType")+
  scale_color_manual(values = tmp.colors)+
  theme_cowplot(font_size = 28) +
  NoAxes() +
  theme(legend.position = "none",
        legend.justification = "center",
        aspect.ratio = 1,
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
p2
p1 | p2
ggsave(filename = myFileName(prefix = "res/fig/fig5_moncle3_test",
                             suffix = ".png"),
       width = 8,height = 4,dpi = 350)

#### load integrated umap
tmp.cell.id <- rownames(cds@int_colData$reducedDims$UMAP)
DefaultAssay(tmp.seu) <- "integrated"
tmp.umap.embed <- Embeddings(tmp.seu,reduction = "umap")[tmp.cell.id,]
cds@int_colData$reducedDims$UMAP <- tmp.umap.embed

p1 <- plot_cells(cds,
                 cell_size = 1.5,
                 group_label_size = 6,
                 reduction_method = "UMAP",
                 label_cell_groups = F,
                 color_cells_by = "orig.ident")+
  theme_cowplot(font_size = 28) +
  scale_color_manual(values = c("#2369ad","#b9d9e9"))+
  guides(color=guide_legend(nrow = 2))+
  NoAxes() +
  theme(legend.title  = element_blank(),
        legend.position = "top",
        legend.justification = "center",
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

tmp.cell.id <- sample(rownames(pData(cds)))
p2 <- plot_cells(cds[,tmp.cell.id],
                 cell_size = 1.5,
                 group_label_size = 6,
                 reduction_method = "UMAP",
                 color_cells_by = "CellType")+
  scale_color_manual(values = tmp.colors)+
  theme_cowplot(font_size = 28) +
  NoAxes() +
  theme(legend.position = "none",
        legend.justification = "center",
        aspect.ratio = 1,
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
p2
p1 | p2
ggsave(filename = myFileName(prefix = "res/fig/fig5_moncle3_test",
                             suffix = ".png"),
       width = 6,height = 4,dpi = 350)

####clustering
cds <- cluster_cells(cds = cds)

p1 <- plot_cells(cds = cds, 
                 cell_size = 1.5,
                 group_label_size = 9,
                 reduction_method = "UMAP")+
  scale_color_manual(values = divergentcolor(6))+
  ggtitle("by clusterID") +
  theme_cowplot(font_size = 28) +
  NoAxes() +
  theme(legend.position = "none",
        legend.justification = "center",
        aspect.ratio = 1,
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))

p2 <- plot_cells(cds = cds, 
                 cell_size = 1.5,
                 group_label_size = 9,
                 color_cells_by = "partition",
                 reduction_method = "UMAP")+
  scale_color_manual(values = divergentcolor(6))+
  ggtitle("by partionID") + 
  theme_cowplot(font_size = 28) +
  NoAxes() +
  theme(legend.position = "none",
        legend.justification = "center",
        aspect.ratio = 1,
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
p1 | p2
ggsave(filename = myFileName(prefix = "res/fig/fig5_moncle3_by_cluster",
                             suffix = ".png"),
       width = 6,height = 4,dpi = 350)


## learn trajectory
####orignial 
cds <- learn_graph(cds,
                   close_loop = F,
                   learn_graph_control = list(minimal_branch_len = 5,
                                              orthogonal_proj_tip = T))

# cds <- learn_graph(cds,
#                    close_loop = F,
#                    learn_graph_control = list(minimal_branch_len = 5,
#                                               orthogonal_proj_tip = F))


# plot_cells(cds, 
#            color_cells_by = "seurat_clusters", 
#            cell_size = 1.5,
#            group_label_size = 9,
#            label_principal_points = T,
#            label_groups_by_cluster = FALSE, 
#            label_leaves = FALSE, 
#            label_branch_points = T)


head(pData(cds))

tmp.colors <- c("#F2BE58",divergentcolor(7)[1:4],
                              grey.colors(3,start = 0.5,rev = T))

### shuffle to avoid overlay 
### 42 is amazing number
set.seed(42)
tmp.cell.id <- sample(rownames(pData(cds)))
p <- plot_cells(cds[,tmp.cell.id], 
               color_cells_by = "CellType",
               cell_size = 1.5,
               group_label_size = 9,
               label_principal_points = T,
               label_groups_by_cluster = FALSE, 
               label_leaves = FALSE, 
               label_branch_points = T)+
  scale_color_manual(values = tmp.colors)+
  theme_cowplot(font_size = 28) +
  theme(legend.position = "none",
        legend.justification = "center",
        aspect.ratio = 1,
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))
p



####recognize root cells

# cds <- order_cells(cds)
# p + 
#   geom_vline(xintercept = seq(1,2,0.25)) + 
#   geom_hline(yintercept = seq(7,8,0.25))
# 
# ggsave(filename = myFileName(prefix = "res/fig/fig5_moncle3_trajectory_show_root_cell",
#                              suffix = ".png"),
#        width = 6,height = 6,dpi = 350)
# 
# embed <- data.frame(Embeddings(tmp.seu, reduction = "umap"))
# root.cell <- WhichCells(tmp.seu, expression = UMAP_1 > 1.75 & UMAP_1 < 2 & UMAP_2 > 7 & UMAP_2 < 7.25)

root.cell <- WhichCells(tmp.seu,idents = "2C_like", expression = UMAP_1 > 0 & UMAP_2 < 9)
cds <- order_cells(cds, root_cells = root.cell)
set.seed(42)
tmp.cell.id <- sample(rownames(pData(cds)))
plot_cells(cds[,tmp.cell.id], 
           color_cells_by = "CellType", 
           cell_size = 1.5,
           group_label_size = 9,
           label_principal_points = T,
           label_groups_by_cluster = FALSE, 
           label_leaves = FALSE, 
           label_roots = F,
           label_branch_points = T)+
  scale_color_manual(values = tmp.colors)+
  theme_cowplot(font_size = 28) +
  theme(legend.position = "none",
        legend.justification = "center",
        aspect.ratio = 1,
        strip.text = element_text(face = "bold"),
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5,
                                  face = "bold"),
        panel.border = element_rect(size = 1,
                                    color = "black"))


cds <- order_cells(cds, root_pr_nodes = c("Y_96"))


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


####save trajectory
saveRDS(object = cds,file = myFileName(prefix = "res/R/B_blastoid_monocle3_cds",suffix = ".rds"))





####---------2.3 construct TRN-------------

####---------2.3.1 load data-------
####load TF
gene.TF <- read.delim(file = "database/Animal_TF_DB/Mus_musculus_TF.txt",stringsAsFactors = F) %>%
  pull(Symbol) %>%
  unique()
####load sc object
seu <- readRDS(file = "res/R/B_blastoid_seurat_inter_2_data_2022042218.rds")
cds <- readRDS(file = "res/R/B_blastoid_monocle3_cds_2022042219.rds")
####check cell types
seu$CellType <- Idents(seu)
####-------2.3.2 run preprocess steps--------
####-------2.3.2.1 run imputation by ALRA method----------
seu <- RunALRA(seu)
saveRDS(seu,file = myFileName(prefix = "res/R/B_blastoid_seurat_RunALRA",suffix = ".rds"))
####Now the default assay is
####DefaultAssay(seu)
####[1] "alra"

# Check the impuation result
# DefaultAssay(seu) <- "RNA"
# p1 <- FeaturePlot(seu,features = c("Cnbp","Zscan4d"),order = T) &
#   scale_color_gradientn(colours = rdwhbu(256))
# DefaultAssay(seu) <- "alra"
# p2 <- FeaturePlot(seu,features = c("Cnbp","Zscan4d"),order = T) &
#   scale_color_gradientn(colours = rdwhbu(256))
# p1 / p2


####-------2.3.2.2 get all de gene in order to filter the master regulators-------
####FindAllMarkers
####This one is required to refiene the master regulator
####Default assay, RNA

tmp.all.de.RNA <- FindAllMarkers(object = seu,
                             assay = "RNA",
                             only.pos = T)
#head(tmp.all.de)
table(tmp.all.de$cluster)

tmp.all.de.ALRA <- FindAllMarkers(object = seu,
                             assay = "alra",
                             only.pos = T)



####


####show all the cell type composition
tmp.df <- seu[[]]
####table(tmp.df$CellType,tmp.df$orig.ident)
table(tmp.df$CellType,tmp.df$orig.ident)

tmp.df.plot <- tmp.df %>%
  group_by(orig.ident,CellType) %>%
  summarise(n = n()) %>%
  mutate(tmp.ratio = n/sum(n)) %>%
  ungroup()

### plot cell type composition
tmp.colors <- c("#F2BE58",divergentcolor(7)[1:4],
                grey.colors(3,start = 0.5,rev = T))
ggplot(tmp.df.plot,aes(orig.ident,tmp.ratio,
                       fill=CellType))+
  geom_bar(position = "stack",stat = "identity")+
  scale_fill_manual(name = "Cell Type",values = tmp.colors)+
  scale_y_continuous(expand = c(0,0))+
  ylab("ratio") +
  xlab(NULL) +
  theme_cowplot(font_size = 28)+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1))
ggsave(filename = myFileName(prefix = "res/fig/fig5_B_blastoid_CellType_Composition",suffix = ".png"),
       width = 8,height = 6,dpi = 350)

####-------2.3.2.3 build cell type specific TRN-------

####get the variable gene along the trajectory
pr_graph_test_res <- graph_test(cds = cds,
                                neighbor_graph="knn", 
                                cores=8)
gene.use <- pr_graph_test_res %>%
  dplyr::filter(q_value < 1e-8)
pr_deg_ids <- gene.use %>%
  pull(gene_short_name) %>%
  as.vector() %>%
  as.character()

####use these variable gene as phenotype
tmp.hits <- pr_graph_test_res %>%
  dplyr::filter(q_value < 1e-8) %>%
  pull(gene_short_name) %>%
  as.vector() %>%
  as.character()
tmp.phenotype <- pr_graph_test_res %>%
  dplyr::filter(q_value < 1e-8) %>%
  pull(morans_test_statistic)
names(tmp.phenotype) <- tmp.hits


####get the imputation data
tmp.data <- GetAssayData(object = seu,slot = "data",assay = "alra")
tmp.data.alra <- as.matrix(tmp.data)


### test build matrix 
tmp.levels <- levels(tmp.df$CellType)
tmp.ident <- "2C_like"

tmp.res.list <- list()
for(ii in levels(seu)){
  cat(ii,sep = "\n")
  tmp.ident <- ii
  cat(paste0(tmp.ident,":start"),sep = "\n")
  tmp.df <- seu[[]]
  tmp.cell.id <- tmp.df %>%
    filter(CellType == tmp.ident) %>%
    rownames()
  ### run tni pipeline
  rtni <- tni.constructor(expData = tmp.data.alra[,tmp.cell.id],
                          regulatoryElements = gene.TF)
  # Please set nPermutations >= 1000
  options(cluster=snow::makeCluster(spec=10, "SOCK"))
  rtni <- tni.permutation(rtni,nPermutations = 1000,pValueCutoff = 0.01)
  rtni <- tni.bootstrap(rtni)
  stopCluster(getOption("cluster"))
  # Compute the DPI-filtered regulatory network
  rtni <- tni.dpi.filter(rtni)
  saveRDS(rtni,file = myFileName(prefix = paste0("res/R/",tmp.ident,"RTN_rtni_object"),suffix = ".rds"))
  
  ### run tna pipeline
  rtna <- tni2tna.preprocess(rtni, phenotype = tmp.phenotype,hits = tmp.hits)
  # Run the MRA method
  rtna <- tna.mra(rtna)
  saveRDS(rtna,file = myFileName(prefix = paste0("res/R/",tmp.ident,"RTN_rtna_object"),suffix = ".rds"))
  
  # # Get MRA results;
  # #..setting 'ntop = -1' will return all results, regardless of a threshold
  # mra <- tna.get(rtna, what="mra", ntop = -1)
  # 
  # # use differential gene to filter
  # tmp.gene.de.use <- tmp.all.de %>%
  #   filter(cluster == tmp.ident) %>%
  #   filter(p_val_adj < 0.05) %>%
  #   rownames()
  # mra.res <- mra %>%
  #   filter(Pvalue < 0.05) %>%
  #   filter(Regulon %in% tmp.gene.de.use)
  # 
  # tmp.res.list <- c(tmp.res.list,list(mra.res))
}


####read result

tmp.res.list <- lapply(levels(seu),FUN = function(ii){
  cat(ii,sep = "\n")
  tmp.ident <- ii
  tmp.files <- list.files(path = "res/R",pattern = paste0(tmp.ident,"RTN_rtna_object"))
  rtna <- readRDS(file = paste0("res/R/",tmp.files))
  mra <- tna.get(rtna, what="mra", ntop = -1)
  tmp.gene.use <- tmp.all.de.RNA %>%
    filter(cluster == tmp.ident) %>%
    filter(p_val_adj < 0.05) %>%
    pull(gene)
  mra.res <- mra %>%
    filter(Pvalue < 0.05) %>%
    filter(Regulon %in% tmp.gene.use) %>% 
    mutate(group = tmp.ident)
  return(mra.res)
})


tmp.mra.res <- Reduce(rbind,tmp.res.list)
saveRDS(object = tmp.mra.res,file = "res/R/B_blastoid_RTN_mra_result.rds")
table(tmp.mra.res$group)

# rtna <- readRDS(file = "res/R/TERTN_rtna_object_2022042600.rds")
# mra <- tna.get(rtna, what="mra", ntop = -1)
# tf.net <- tna.get(rtna, what="tnet", ntop = -1)
# tf.regulon <- tna.get(rtna,what = "regulons")
# tf.use <- tna.get(rtna,what = "regulatoryElements")
# which(tf.use == "")
# x <- tf.net["Klf6",]
# x[x!=0]
# intersect(tf.regulon$Klf6,tf.use)

####---------2.3.2.4 infer temporal variable LR-----------------

####---------2.3.2.4.1 check by pseudotime-----------
####prepare the psedotime object
tmp.df <- data.frame(pseudotime = pseudotime(cds,reduction_method = "UMAP"),
                     stringsAsFactors = F)
###pseudotime() was the function to acesss the pseudotime 
####range(cds@principal_graph_aux$UMAP$pseudotime)
tmp.df <- cbind(seu[[]],tmp.df)

####Let's plot pseudotime pheatmap monocle3

genes <- pr_graph_test_res %>%
  filter(morans_I > 0.25 & q_value == 0) %>%
  pull(gene_short_name) %>%
  as.vector() %>%
  as.character()
  
# use normalized_cousnt matrix 
pt.matrix <- normalized_counts(cds, norm_method = "log")[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#or exprs(cds)
pt.matrix <- as.matrix(pt.matrix)
tmp.colnames <- colnames(pt.matrix)
# smooth the matrix
pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes
colnames(pt.matrix) <- tmp.colnames
head(pt.matrix[1:10,1:10])

range(pt.matrix)
tmp.colannotation <- tmp.df %>%
  select(orig.ident,CellType,pseudotime)


tmp.colors.ident <- c("#2369ad","#b9d9e9")
names(tmp.colors.ident) <- c("blastocyst","EPS_blastoid")
tmp.colors.celltype <- c("#F2BE58",divergentcolor(7)[1:4],
                         grey.colors(3,start = 0.5,rev = T))
names(tmp.colors.celltype) <- levels(seu)
ann_colors = list(
  pseudotime = warm(100),
  orig.ident = tmp.colors.ident,
  CellType =  tmp.colors.celltype)

pt.matrix <- pt.matrix[rownames(pt.matrix),]
tmp_labels <-  c("Fgf4","Fgfr1","Gata6")
tmp_idx <- which(rownames(pt.matrix) %in% tmp_labels)
rownames(pt.matrix)[tmp_idx]

tmp_anno <- anno_mark(at=tmp_idx,
                      labels=rownames(pt.matrix)[tmp_idx],
                      which = "row",
                      labels_gp = gpar(fontsize=18,
                                       fontface="italic"))


ph <- pheatmap_fixed(mat = pt.matrix,
               color = coldwarm(100),
               show_rownames = F,
               na_col = "black",
               annotation_colors = ann_colors,
               show_colnames = F,
               fontsize = 16,
               cluster_cols = F,
               name = "zscore",
               annotation_col = tmp.colannotation,
               clustering_method = "ward.D2")+
  rowAnnotation(mark = tmp_anno)



png(filename = myFileName(prefix = "res/fig/fig5_test_pseudotime",
                          suffix = ".png"),
    width = 6,height = 8,units = "in",res = 300)
ph
dev.off()

pdf(file = myFileName(prefix = "res/fig/fig5_test_pseudotime",
                          suffix = ".pdf"),
    width = 6,height = 8)
ph
dev.off()

#####----2.3.2.4.2 check by correlation--------------
##### This step was used to set up the correlation of eLR with the core TF
tmp.df <- data.frame(pseudotime = pseudotime(cds,reduction_method = "UMAP"),
                     stringsAsFactors = F)
tmp.df <- cbind(seu[[]],tmp.df)
tmp.data <-  GetAssayData(object = seu,slot = "data",assay = "RNA")

tmp.orig.ident <- unique(tmp.df$orig.ident)[1]
tmp.ident.1 <- "EPI"
tmp.ident.2 <- "PE"

tmp.cell.meta.1 <- tmp.df %>%
  rownames_to_column("cell_id") %>%
  filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.1) %>% 
  #sample_n(size = 100,replace = T) %>%
  arrange(pseudotime)

tmp.cell.meta.2 <- tmp.df %>%
  rownames_to_column("cell_id") %>%
  filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.2) %>%
  #sample_n(size = 100,replace = T) %>%
  arrange(pseudotime)


LRpairs.df <- read.delim(file = "database/Ligand-Receptor-Pairs/Mouse/Mouse-2020-Shao-LR-pairs.txt",
                         stringsAsFactors = F)
LRpairs <- LRpairs.df$lr_pair
Lgenelist <- LRpairs.df$ligand_gene_symbol
Rgenelist <- LRpairs.df$receptor_gene_symbol 

gene_symbols <- rownames(tmp.data)
l.remove <- setdiff(Lgenelist,gene_symbols)
r.remove <- setdiff(Rgenelist,gene_symbols)
index.remove <- c(which(Lgenelist %in% l.remove),which(Rgenelist %in% r.remove))
LRpairs <- LRpairs[-index.remove]
Lgenelist <- Lgenelist[-index.remove]
Rgenelist <- Rgenelist[-index.remove]



#####--------2.3.2.4.3 using cell-----------
tmp.scale.function <- function(x){
  (x-min(x))/(max(x)-min(x))
}

tmp.scale.function(1:3)



tmp.orig.ident <- unique(tmp.df$orig.ident)[1]
#### test principle
tmp.ident.1 <- "EPI"
tmp.ident.2 <- "PE"
tmp.cell.meta.1 <- tmp.df %>%
  rownames_to_column("cell_id") %>%
  dplyr::filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.1) %>% 
  arrange(pseudotime) 
x <- tmp.cell.meta.1$pseudotime 
tmp.cell.meta.1$pseudotime <- tmp.scale.function(x)
range(tmp.cell.meta.1$pseudotime)


tmp.cell.meta.2 <- tmp.df %>%
  rownames_to_column("cell_id") %>%
  filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.2) %>%
  arrange(pseudotime)
x <- tmp.cell.meta.2$pseudotime 
range(x)
tmp.cell.meta.2$pseudotime <- tmp.scale.function(x)
range(tmp.cell.meta.2$pseudotime)

tmp.mat.1 <- tmp.data[,tmp.cell.meta.1$cell_id]
tmp.mat.2 <- tmp.data[,tmp.cell.meta.2$cell_id]

tmp.mat.pseudotime.1 <- tmp.cell.meta.1$pseudotime
tmp.mat.pseudotime.2 <- tmp.cell.meta.2$pseudotime

names(tmp.mat.pseudotime.1) <- tmp.cell.meta.1$cell_id
names(tmp.mat.pseudotime.2) <- tmp.cell.meta.2$cell_id


Thresh=0.2
numPts = 200
inter.tmp.mat.1 = interWeights(expDataBatch = tmp.mat.1, 
                               trajCond = tmp.mat.pseudotime.1, 
                               winSz = 0.1, 
                               numPts = numPts)
inter.tmp.mat.2 = interWeights(expDataBatch = tmp.mat.2, 
                               trajCond = tmp.mat.pseudotime.2, 
                               winSz = 0.1, 
                               numPts = numPts)
inter.tmp.mat.1 = cellAlign::scaleInterpolate(inter.tmp.mat.1)
inter.tmp.mat.2 = cellAlign::scaleInterpolate(inter.tmp.mat.2)
time <- inter.tmp.mat.1$traj


x <- inter.tmp.mat.1$scaledData["Fgf4",]
y <- inter.tmp.mat.2$scaledData["Fgfr2",]

tmp.data.plot <- data.frame(exp = c(x,y,sqrt(x*y)),
                            group=c(rep("L",length(x)),
                                    rep("R",length(y)),
                                    rep("IS",length(y))),
                            time=c(time,time,time),
                            stringsAsFactors = F) %>%
  mutate(group = factor(group,levels = c("L","R","IS")))

ggplot(tmp.data.plot,aes(time,exp,color=group))+
  geom_point(size=1.5)+
  scale_color_manual(values = ggsci::pal_npg()(3),
                     guide=guide_legend(title = NULL))+
  theme_cowplot(font_size = 38)+
  NoAxes()+
  theme(legend.position = "top",legend.justification = "center")
  

ggsave(filename = myFileName(prefix = "res/fig/fig5_Talkmodel_cor_var_LR_case",suffix = ".png"),
       width = 6,height = 6,dpi = 350,bg = "white")

tmp.levels <- tmp.df %>%
  filter(orig.ident == tmp.orig.ident) %>%
  pull(CellType) %>%
  as.character() %>%
  as.vector() %>%
  unique()
tmp.cell.comb.df <- expand.grid(cell.L = tmp.levels, 
                                cell.R = tmp.levels,
                                stringsAsFactors = F)
tmp.res.TimeTalk.list <- lapply(1:nrow(tmp.cell.comb.df), function(ii){
  
  tmp.ident.1 <- tmp.cell.comb.df$cell.L[ii]
  tmp.ident.2 <- tmp.cell.comb.df$cell.R[ii]
  cat(paste0(tmp.ident.1,"-",tmp.ident.2),sep = "\n")
  tmp.cell.meta.1 <- tmp.df %>%
    rownames_to_column("cell_id") %>%
    dplyr::filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.1) %>% 
    arrange(pseudotime) 
  x <- tmp.cell.meta.1$pseudotime 
  tmp.cell.meta.1$pseudotime <- tmp.scale.function(x)
  range(tmp.cell.meta.1$pseudotime)
  
  
  tmp.cell.meta.2 <- tmp.df %>%
    rownames_to_column("cell_id") %>%
    filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.2) %>%
    arrange(pseudotime)
  x <- tmp.cell.meta.2$pseudotime 
  range(x)
  tmp.cell.meta.2$pseudotime <- tmp.scale.function(x)
  range(tmp.cell.meta.2$pseudotime)
  
  tmp.mat.1 <- tmp.data[,tmp.cell.meta.1$cell_id]
  tmp.mat.2 <- tmp.data[,tmp.cell.meta.2$cell_id]
  
  tmp.mat.pseudotime.1 <- tmp.cell.meta.1$pseudotime
  tmp.mat.pseudotime.2 <- tmp.cell.meta.2$pseudotime
  
  names(tmp.mat.pseudotime.1) <- tmp.cell.meta.1$cell_id
  names(tmp.mat.pseudotime.2) <- tmp.cell.meta.2$cell_id
  
  
  Thresh=0.2
  numPts = 200
  inter.tmp.mat.1 = interWeights(expDataBatch = tmp.mat.1, 
                                 trajCond = tmp.mat.pseudotime.1, 
                                 winSz = 0.1, 
                                 numPts = numPts)
  inter.tmp.mat.2 = interWeights(expDataBatch = tmp.mat.2, 
                                 trajCond = tmp.mat.pseudotime.2, 
                                 winSz = 0.1, 
                                 numPts = numPts)
  inter.tmp.mat.1 = cellAlign::scaleInterpolate(inter.tmp.mat.1)
  inter.tmp.mat.2 = cellAlign::scaleInterpolate(inter.tmp.mat.2)
  time <- inter.tmp.mat.1$traj
  
  
  plan("multisession",workers = 10)
  tmp.res.list <- future_lapply(seq_along(Lgenelist),FUN = function(ii){
    cat(ii,sep = "\n")
    x <- inter.tmp.mat.1$scaledData[Lgenelist[ii],]
    y <- inter.tmp.mat.2$scaledData[Rgenelist[ii],]
    tmp.res <- cor(x,y)
    return(tmp.res)
  })
  plan("sequential")
  names(tmp.res.list) <- paste0(Lgenelist,"-",Rgenelist)
  
  tmp.res.cor <- data.frame(PCC = unlist(tmp.res.list),stringsAsFactors = F) %>%
    rownames_to_column("LRpairs") %>%
    mutate(PCC = ifelse(is.na(PCC),0,PCC)) %>%
    arrange(-PCC) %>%
    mutate(Rank = row_number())
  
  # tmp.res.cor %>%
  #   filter(abs(PCC) > 0.5) %>%
  #   nrow()
  
  # ggplot(tmp.res.cor,aes(Rank,PCC))+
  #   geom_point()+
  #   theme_cowplot(font_size = 28)
  # ggsave(filename = myFileName(prefix =  "res/fig/test_TimeTalk_PCC_rank_plot",
  #                              suffix = ".png"),
  #        width = 8,height = 8,dpi = 350)
  
  
  #### read master regulator analysis
  tmp.mra.res <- readRDS(file = "res/R/B_blastoid_RTN_mra_result.rds")
  tmp.TF.gene <- tmp.mra.res %>% 
    filter(group == tmp.ident.2) %>%
    pull(Regulon)
  # #### Let's show the TF levels
  # ph <- pheatmap_fixed(inter.tmp.mat.2$scaledData[tmp.TF.gene,],
  #                      show_colnames = F,
  #                      show_rownames = T,
  #                      cluster_cols = F,
  #                      color = coldwarm(100),
  #                      name = "scaled values")
  # ph
  
  
  tmp.LR.list <- tmp.res.cor %>%
    filter(abs(PCC) > 0.5) %>%
    pull(LRpairs)
  
  plan("multisession",workers = 10)
  tmp.ttt.res <- future_lapply(tmp.LR.list, function(tmp.LR){
    cat(tmp.LR,sep = "\n")
    tmp.L.gene <- unlist(lapply(strsplit(tmp.LR,split = "-"),FUN = function(ii){ ii[1]}))
    tmp.R.gene <- unlist(lapply(strsplit(tmp.LR,split = "-"),FUN = function(ii){ ii[2]}))
    
    x <- inter.tmp.mat.1$scaledData[tmp.L.gene,]
    y <- inter.tmp.mat.2$scaledData[tmp.R.gene,]
    
    tmp.res.list <- lapply(tmp.TF.gene, function(ii){
      cat(paste0(tmp.ident.2,"_TF:",ii),sep = "\n")
      tmp.TF.level <- inter.tmp.mat.2$scaledData[ii,]
      tmp.IS <- sqrt(x*y)
      tmp.res <- grangertest(tmp.IS,tmp.TF.level)
      tmp.res.LRtoTF.pvalue <- tmp.res$`Pr(>F)`[2]
      tmp.res <- grangertest(tmp.TF.level,tmp.IS)
      tmp.res.TFtoLR.pvalue <- tmp.res$`Pr(>F)`[2]
      tmp.PCC <- cor(tmp.IS,tmp.TF.level,method = "pearson")
      tmp.SCC <- cor(tmp.IS,tmp.TF.level,method = "spearman")
      tmp.res.df <- data.frame(LR=tmp.LR,
                               L=tmp.L.gene,
                               R=tmp.R.gene,
                               TF=ii,
                               PCC=tmp.PCC,
                               SCC=tmp.SCC,
                               LRtoTF=tmp.res.LRtoTF.pvalue,
                               TFtoLR=tmp.res.TFtoLR.pvalue,
                               cell.L=tmp.ident.1,
                               cell.R=tmp.ident.2,
                               stringsAsFactors = F)
      return(tmp.res.df)
    })
    tmp.res.df <- Reduce(rbind,tmp.res.list)
    tmp.res.df <- tmp.res.df %>%
      mutate(category=ifelse(LRtoTF < 1e-6 | TFtoLR < 1e-6,"PASS","SKIP"))
    return(tmp.res.df)
  })
  plan("sequential")
  
  tmp.ttt.res.df <- Reduce(rbind,tmp.ttt.res)
  
  tmp.res.df <- tmp.ttt.res.df %>%
    filter(category != "SKIP") %>%
    filter(abs(PCC) > 0.5)
  tmp.res.df$orig.ident <- tmp.orig.ident
  
  return(tmp.res.df)
})
TimeTalk.res <- Reduce(rbind,tmp.res.TimeTalk.list)
TimeTalk.res.EPS_blastoid <- TimeTalk.res



tmp.orig.ident <- "blastocyst"
# tmp.ident.1 <- "EPI"
# tmp.ident.2 <- "PE"
tmp.levels <- tmp.df %>%
  filter(orig.ident == tmp.orig.ident) %>%
  pull(CellType) %>%
  as.character() %>%
  as.vector() %>%
  unique()

tmp.cell.comb.df <- expand.grid(cell.L = tmp.levels, 
                                cell.R = tmp.levels,
                                stringsAsFactors = F)

tmp.cell.comb.df <- tmp.cell.comb.df %>%
  filter(cell.L == "EPI" & cell.R == "PE")


tmp.res.TimeTalk.list <- lapply(1:nrow(tmp.cell.comb.df), function(ii){
  tmp.ident.1 <- tmp.cell.comb.df$cell.L[ii]
  tmp.ident.2 <- tmp.cell.comb.df$cell.R[ii]
  cat(paste0(tmp.ident.1,"-",tmp.ident.2),sep = "\n")
  tmp.cell.meta.1 <- tmp.df %>%
    rownames_to_column("cell_id") %>%
    dplyr::filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.1) %>% 
    arrange(pseudotime) 
  x <- tmp.cell.meta.1$pseudotime 
  tmp.cell.meta.1$pseudotime <- tmp.scale.function(x)
  range(tmp.cell.meta.1$pseudotime)
  
  tmp.cell.meta.2 <- tmp.df %>%
    rownames_to_column("cell_id") %>%
    dplyr::filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.2) %>%
    arrange(pseudotime)
  x <- tmp.cell.meta.2$pseudotime 
  range(x)
  tmp.cell.meta.2$pseudotime <- tmp.scale.function(x)
  range(tmp.cell.meta.2$pseudotime)
  
  tmp.mat.1 <- tmp.data[,tmp.cell.meta.1$cell_id]
  tmp.mat.2 <- tmp.data[,tmp.cell.meta.2$cell_id]
  
  tmp.mat.pseudotime.1 <- tmp.cell.meta.1$pseudotime
  tmp.mat.pseudotime.2 <- tmp.cell.meta.2$pseudotime
  
  names(tmp.mat.pseudotime.1) <- tmp.cell.meta.1$cell_id
  names(tmp.mat.pseudotime.2) <- tmp.cell.meta.2$cell_id
  
  
  Thresh=0.2
  numPts = 200
  inter.tmp.mat.1 = interWeights(expDataBatch = tmp.mat.1, 
                                 trajCond = tmp.mat.pseudotime.1, 
                                 winSz = 0.1 , 
                                 numPts = numPts)
  inter.tmp.mat.2 = interWeights(expDataBatch = tmp.mat.2, 
                                 trajCond = tmp.mat.pseudotime.2, 
                                 winSz = 0.1 ,
                                 numPts = numPts)
  inter.tmp.mat.1 = cellAlign::scaleInterpolate(inter.tmp.mat.1)
  inter.tmp.mat.2 = cellAlign::scaleInterpolate(inter.tmp.mat.2)
  time <- inter.tmp.mat.1$traj
  
  inter.tmp.mat.1 <- myRemoveNA(inter.tmp.mat.1$scaledData)
  inter.tmp.mat.2 <- myRemoveNA(inter.tmp.mat.2$scaledData)
  

  plan("multisession",workers = 10)
  tmp.res.list <- future_lapply(seq_along(Lgenelist),FUN = function(ii){
    cat(ii,sep = "\n")
    x <- inter.tmp.mat.1[Lgenelist[ii],]
    y <- inter.tmp.mat.2[Rgenelist[ii],]
    tmp.res <- cor(x,y)
    return(tmp.res)
  })
  plan("sequential")
  names(tmp.res.list) <- paste0(Lgenelist,"-",Rgenelist)
  
  tmp.res.cor <- data.frame(PCC = unlist(tmp.res.list),stringsAsFactors = F) %>%
    rownames_to_column("LRpairs") %>%
    mutate(PCC = ifelse(is.na(PCC),0,PCC)) %>%
    arrange(-PCC) %>%
    mutate(Rank = row_number())
  
  # tmp.res.cor %>%
  #   filter(abs(PCC) > 0.5) %>%
  #   nrow()
  
  # ggplot(tmp.res.cor,aes(Rank,PCC))+
  #   geom_point()+
  #   theme_cowplot(font_size = 28)
  # ggsave(filename = myFileName(prefix =  "res/fig/test_TimeTalk_PCC_rank_plot",
  #                              suffix = ".png"),
  #        width = 8,height = 8,dpi = 350)
  
  
  #### read master regulator analysis
  tmp.mra.res <- readRDS(file = "res/R/B_blastoid_RTN_mra_result.rds")
  tmp.TF.gene <- tmp.mra.res %>% 
    filter(group == tmp.ident.2) %>%
    pull(Regulon)
  # #### Let's show the TF levels
  # ph <- pheatmap_fixed(inter.tmp.mat.2$scaledData[tmp.TF.gene,],
  #                      show_colnames = F,
  #                      show_rownames = T,
  #                      cluster_cols = F,
  #                      color = coldwarm(100),
  #                      name = "scaled values")
  # ph
  
  
  tmp.LR.list <- tmp.res.cor %>%
    filter(abs(PCC) > 0.5) %>%
    pull(LRpairs)
  
  plan("multisession",workers = 10)
  tmp.ttt.res <- future_lapply(tmp.LR.list, function(tmp.LR){
    cat(tmp.LR,sep = "\n")
    tmp.L.gene <- unlist(lapply(strsplit(tmp.LR,split = "-"),FUN = function(ii){ ii[1]}))
    tmp.R.gene <- unlist(lapply(strsplit(tmp.LR,split = "-"),FUN = function(ii){ ii[2]}))
    
    x <- inter.tmp.mat.1[tmp.L.gene,]
    y <- inter.tmp.mat.2[tmp.R.gene,]
    
    tmp.res.list <- lapply(tmp.TF.gene, function(tmp.TF.gene.use){
      cat(paste0(tmp.ident.2,"_TF:",tmp.TF.gene.use),sep = "\n")
      tmp.TF.level <- inter.tmp.mat.2[tmp.TF.gene.use,]
      tmp.IS <- sqrt(x*y)
      tmp.res.LRtoTF.pvalue <- tryCatch(
            expr = {
              tmp.res <- grangertest(tmp.IS,tmp.TF.level)
              tmp.res$`Pr(>F)`[2]},
            error = function(e) {
              1
            })
      
      tmp.res.TFtoLR.pvalue <- tryCatch(
            expr = {
              tmp.res <- grangertest(tmp.TF.level,tmp.IS)
              tmp.res$`Pr(>F)`[2]},
            error = function(e) {
              1
            })
      tmp.PCC <- cor(tmp.IS,tmp.TF.level,method = "pearson")
      tmp.PCC <- ifelse(is.na(tmp.PCC),0,tmp.PCC)
      tmp.SCC <- cor(tmp.IS,tmp.TF.level,method = "spearman")
      tmp.SCC <- ifelse(is.na(tmp.SCC),0,tmp.SCC)
      tmp.res.df <- data.frame(LR=tmp.LR,
                               L=tmp.L.gene,
                               R=tmp.R.gene,
                               TF=tmp.TF.gene.use,
                               PCC=tmp.PCC,
                               SCC=tmp.SCC,
                               LRtoTF=tmp.res.LRtoTF.pvalue,
                               TFtoLR=tmp.res.TFtoLR.pvalue,
                               cell.L=tmp.ident.1,
                               cell.R=tmp.ident.2,
                               stringsAsFactors = F)
      return(tmp.res.df)
    })
    tmp.res.df <- Reduce(rbind,tmp.res.list)
    tmp.res.df <- tmp.res.df %>%
      mutate(category=ifelse(LRtoTF < 1e-6 | TFtoLR < 1e-6,"PASS","SKIP"))
    return(tmp.res.df)
  })
  plan("sequential")
  tmp.ttt.res.df <- Reduce(rbind,tmp.ttt.res)
  tmp.res.df <- tmp.ttt.res.df %>%
    filter(category != "SKIP") %>%
    filter(abs(PCC) > 0.5)
  return(tmp.res.df)
})

TimeTalk.res <- Reduce(rbind,tmp.res.TimeTalk.list)
TimeTalk.res.blastocyst <- TimeTalk.res
TimeTalk.res.blastocyst$orig.ident <- tmp.orig.ident
TimeTalk.res.all <- rbind(TimeTalk.res.blastocyst,TimeTalk.res.EPS_blastoid)
saveRDS(TimeTalk.res.all,
        file = myFileName(prefix = "res/R/TimeTalk_res_blastocyst",suffix = ".rds"))


tmp.df <- data.frame(pseudotime = pseudotime(cds,reduction_method = "UMAP"),
                     stringsAsFactors = F)
tmp.df <- cbind(seu[[]],tmp.df)

#####---------2.3.2.4.4 show EPI_PE blastocyst-------------------------

tmp.orig.ident <- "blastocyst"
tmp.data <-  GetAssayData(object = seu,slot = "data",assay = "RNA")
tmp.cell.comb.df <- tmp.cell.comb.df %>%
  filter(cell.L == "EPI" & cell.R == "PE")

tmp.res.TimeTalk.list <- lapply(1:nrow(tmp.cell.comb.df), function(ii){
  tmp.ident.1 <- tmp.cell.comb.df$cell.L[ii]
  tmp.ident.2 <- tmp.cell.comb.df$cell.R[ii]
  cat(paste0(tmp.ident.1,"-",tmp.ident.2),sep = "\n")
  tmp.cell.meta.1 <- tmp.df %>%
    rownames_to_column("cell_id") %>%
    dplyr::filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.1) %>% 
    arrange(pseudotime) 
  x <- tmp.cell.meta.1$pseudotime 
  tmp.cell.meta.1$pseudotime <- tmp.scale.function(x)
  range(tmp.cell.meta.1$pseudotime)
  
  tmp.cell.meta.2 <- tmp.df %>%
    rownames_to_column("cell_id") %>%
    dplyr::filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.2) %>%
    arrange(pseudotime)
  x <- tmp.cell.meta.2$pseudotime 
  range(x)
  tmp.cell.meta.2$pseudotime <- tmp.scale.function(x)
  range(tmp.cell.meta.2$pseudotime)
  
  tmp.mat.1 <- tmp.data[,tmp.cell.meta.1$cell_id]
  tmp.mat.2 <- tmp.data[,tmp.cell.meta.2$cell_id]
  
  tmp.mat.pseudotime.1 <- tmp.cell.meta.1$pseudotime
  tmp.mat.pseudotime.2 <- tmp.cell.meta.2$pseudotime
  
  names(tmp.mat.pseudotime.1) <- tmp.cell.meta.1$cell_id
  names(tmp.mat.pseudotime.2) <- tmp.cell.meta.2$cell_id
  
  
  Thresh=0.2
  numPts = 200
  inter.tmp.mat.1 = interWeights(expDataBatch = tmp.mat.1, 
                                 trajCond = tmp.mat.pseudotime.1, 
                                 winSz = 0.1 , 
                                 numPts = numPts)
  inter.tmp.mat.2 = interWeights(expDataBatch = tmp.mat.2, 
                                 trajCond = tmp.mat.pseudotime.2, 
                                 winSz = 0.1 ,
                                 numPts = numPts)
  inter.tmp.mat.1 = cellAlign::scaleInterpolate(inter.tmp.mat.1)
  inter.tmp.mat.2 = cellAlign::scaleInterpolate(inter.tmp.mat.2)
  time <- inter.tmp.mat.1$traj
  
  inter.tmp.mat.1 <- myRemoveNA(inter.tmp.mat.1$scaledData)
  inter.tmp.mat.2 <- myRemoveNA(inter.tmp.mat.2$scaledData)
  
  
  plan("multisession",workers = 10)
  tmp.res.list <- future_lapply(seq_along(Lgenelist),FUN = function(ii){
    cat(ii,sep = "\n")
    x <- inter.tmp.mat.1[Lgenelist[ii],]
    y <- inter.tmp.mat.2[Rgenelist[ii],]
    tmp.res <- cor(x,y,method = "pearson")
    return(tmp.res)
  })
  names(tmp.res.list) <- paste0(Lgenelist,"-",Rgenelist)
  tmp.res.list.PCC <- tmp.res.list
  
  tmp.res.list <- future_lapply(seq_along(Lgenelist),FUN = function(ii){
    cat(ii,sep = "\n")
    x <- inter.tmp.mat.1[Lgenelist[ii],]
    y <- inter.tmp.mat.2[Rgenelist[ii],]
    tmp.res <- cor(x,y,method = "spearman")
    return(tmp.res)
  })
  names(tmp.res.list) <- paste0(Lgenelist,"-",Rgenelist)
  tmp.res.list.SCC <- tmp.res.list
  plan("sequential")
  
  
  tmp.res.cor <- data.frame(PCC = unlist(tmp.res.list),
                            SCC = unlist(tmp.res.list),
                            stringsAsFactors = F) %>%
    rownames_to_column("LRpairs") %>%
    mutate(PCC = ifelse(is.na(PCC),0,PCC)) %>%
    mutate(SCC = ifelse(is.na(SCC),0,SCC)) %>%
    arrange(-SCC) %>%
    mutate(Rank = row_number())
  
  # tmp.res.cor %>%
  #   filter(abs(PCC) > 0.5) %>%
  #   nrow()
  
  # ggplot(tmp.res.cor,aes(Rank,PCC))+
  #   geom_point()+
  #   theme_cowplot(font_size = 28)
  # ggsave(filename = myFileName(prefix =  "res/fig/test_TimeTalk_PCC_rank_plot",
  #                              suffix = ".png"),
  #        width = 8,height = 8,dpi = 350)
  
  
  #### read master regulator analysis
  tmp.mra.res <- readRDS(file = "res/R/B_blastoid_RTN_mra_result.rds")
  tmp.TF.gene <- tmp.mra.res %>% 
    filter(group == tmp.ident.2) %>%
    pull(Regulon)
  # #### Let's show the TF levels
  # ph <- pheatmap_fixed(inter.tmp.mat.2$scaledData[tmp.TF.gene,],
  #                      show_colnames = F,
  #                      show_rownames = T,
  #                      cluster_cols = F,
  #                      color = coldwarm(100),
  #                      name = "scaled values")
  # ph
  
  
  tmp.LR.list <- tmp.res.cor %>%
    filter(abs(SCC) > 0.2) %>%
    pull(LRpairs)
  
  plan("multisession",workers = 10)
  tmp.ttt.res <- future_lapply(tmp.LR.list, function(tmp.LR){
    cat(tmp.LR,sep = "\n")
    tmp.L.gene <- unlist(lapply(strsplit(tmp.LR,split = "-"),FUN = function(ii){ ii[1]}))
    tmp.R.gene <- unlist(lapply(strsplit(tmp.LR,split = "-"),FUN = function(ii){ ii[2]}))
    
    x <- inter.tmp.mat.1[tmp.L.gene,]
    y <- inter.tmp.mat.2[tmp.R.gene,]
    
    tmp.res.list <- lapply(tmp.TF.gene, function(tmp.TF.gene.use){
      cat(paste0(tmp.ident.2,"_TF:",tmp.TF.gene.use),sep = "\n")
      tmp.TF.level <- inter.tmp.mat.2[tmp.TF.gene.use,]
      tmp.IS <- sqrt(x*y)
      tmp.res.LRtoTF.pvalue <- tryCatch(
        expr = {
          tmp.res <- grangertest(tmp.IS,tmp.TF.level)
          tmp.res$`Pr(>F)`[2]},
        error = function(e) {
          1
        })
      
      tmp.res.TFtoLR.pvalue <- tryCatch(
        expr = {
          tmp.res <- grangertest(tmp.TF.level,tmp.IS)
          tmp.res$`Pr(>F)`[2]},
        error = function(e) {
          1
        })
      tmp.PCC <- cor(tmp.IS,tmp.TF.level,method = "pearson")
      tmp.PCC <- ifelse(is.na(tmp.PCC),0,tmp.PCC)
      tmp.SCC <- cor(tmp.IS,tmp.TF.level,method = "spearman")
      tmp.SCC <- ifelse(is.na(tmp.SCC),0,tmp.SCC)
      tmp.res.df <- data.frame(LR=tmp.LR,
                               L=tmp.L.gene,
                               R=tmp.R.gene,
                               TF=tmp.TF.gene.use,
                               PCC=tmp.PCC,
                               SCC=tmp.SCC,
                               LRtoTF=tmp.res.LRtoTF.pvalue,
                               TFtoLR=tmp.res.TFtoLR.pvalue,
                               cell.L=tmp.ident.1,
                               cell.R=tmp.ident.2,
                               stringsAsFactors = F)
      return(tmp.res.df)
    })
    tmp.res.df <- Reduce(rbind,tmp.res.list)
    tmp.res.df <- tmp.res.df %>%
      mutate(category=ifelse(LRtoTF < 1e-2 | TFtoLR < 1e-2,"PASS","SKIP"))
    return(tmp.res.df)
  })
  plan("sequential")
  tmp.ttt.res.df <- Reduce(rbind,tmp.ttt.res)
  tmp.res.df <- tmp.ttt.res.df 
  return(tmp.res.df)
})


TimeTalk.res <- Reduce(rbind,tmp.res.TimeTalk.list)
TimeTalk.res.blastocyst.EPI_PE <- TimeTalk.res %>%
  filter(category == "PASS")
TimeTalk.res.blastocyst.EPI_PE$orig.ident <- "blastocyst"




####-------2.3.2.4.5 EPS_blastoid-------------
table(tmp.df$orig.ident)
tmp.orig.ident <- "EPS_blastoid"
tmp.data <-  GetAssayData(object = seu,slot = "data",assay = "RNA")
tmp.cell.comb.df <- tmp.cell.comb.df %>%
  filter(cell.L == "EPI" & cell.R == "PE")
tmp.res.TimeTalk.list <- lapply(1:nrow(tmp.cell.comb.df), function(ii){
  tmp.ident.1 <- tmp.cell.comb.df$cell.L[ii]
  tmp.ident.2 <- tmp.cell.comb.df$cell.R[ii]
  cat(paste0(tmp.ident.1,"-",tmp.ident.2),sep = "\n")
  tmp.cell.meta.1 <- tmp.df %>%
    rownames_to_column("cell_id") %>%
    dplyr::filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.1) %>% 
    arrange(pseudotime) 
  x <- tmp.cell.meta.1$pseudotime 
  tmp.cell.meta.1$pseudotime <- tmp.scale.function(x)
  range(tmp.cell.meta.1$pseudotime)
  
  tmp.cell.meta.2 <- tmp.df %>%
    rownames_to_column("cell_id") %>%
    dplyr::filter(orig.ident == tmp.orig.ident & CellType == tmp.ident.2) %>%
    arrange(pseudotime)
  x <- tmp.cell.meta.2$pseudotime 
  range(x)
  tmp.cell.meta.2$pseudotime <- tmp.scale.function(x)
  range(tmp.cell.meta.2$pseudotime)
  
  tmp.mat.1 <- tmp.data[,tmp.cell.meta.1$cell_id]
  tmp.mat.2 <- tmp.data[,tmp.cell.meta.2$cell_id]
  
  tmp.mat.pseudotime.1 <- tmp.cell.meta.1$pseudotime
  tmp.mat.pseudotime.2 <- tmp.cell.meta.2$pseudotime
  
  names(tmp.mat.pseudotime.1) <- tmp.cell.meta.1$cell_id
  names(tmp.mat.pseudotime.2) <- tmp.cell.meta.2$cell_id
  
  
  Thresh=0.2
  numPts = 200
  inter.tmp.mat.1 = interWeights(expDataBatch = tmp.mat.1, 
                                 trajCond = tmp.mat.pseudotime.1, 
                                 winSz = 0.1 , 
                                 numPts = numPts)
  inter.tmp.mat.2 = interWeights(expDataBatch = tmp.mat.2, 
                                 trajCond = tmp.mat.pseudotime.2, 
                                 winSz = 0.1 ,
                                 numPts = numPts)
  inter.tmp.mat.1 = cellAlign::scaleInterpolate(inter.tmp.mat.1)
  inter.tmp.mat.2 = cellAlign::scaleInterpolate(inter.tmp.mat.2)
  time <- inter.tmp.mat.1$traj
  
  inter.tmp.mat.1 <- myRemoveNA(inter.tmp.mat.1$scaledData)
  inter.tmp.mat.2 <- myRemoveNA(inter.tmp.mat.2$scaledData)
  
  
  plan("multisession",workers = 10)
  tmp.res.list <- future_lapply(seq_along(Lgenelist),FUN = function(ii){
    cat(ii,sep = "\n")
    x <- inter.tmp.mat.1[Lgenelist[ii],]
    y <- inter.tmp.mat.2[Rgenelist[ii],]
    tmp.res <- cor(x,y,method = "pearson")
    return(tmp.res)
  })
  names(tmp.res.list) <- paste0(Lgenelist,"-",Rgenelist)
  tmp.res.list.PCC <- tmp.res.list
  
  tmp.res.list <- future_lapply(seq_along(Lgenelist),FUN = function(ii){
    cat(ii,sep = "\n")
    x <- inter.tmp.mat.1[Lgenelist[ii],]
    y <- inter.tmp.mat.2[Rgenelist[ii],]
    tmp.res <- cor(x,y,method = "spearman")
    return(tmp.res)
  })
  names(tmp.res.list) <- paste0(Lgenelist,"-",Rgenelist)
  tmp.res.list.SCC <- tmp.res.list
  plan("sequential")
  
  
  tmp.res.cor <- data.frame(PCC = unlist(tmp.res.list),
                            SCC = unlist(tmp.res.list),
                            stringsAsFactors = F) %>%
    rownames_to_column("LRpairs") %>%
    mutate(PCC = ifelse(is.na(PCC),0,PCC)) %>%
    mutate(SCC = ifelse(is.na(SCC),0,SCC)) %>%
    arrange(-SCC) %>%
    mutate(Rank = row_number())
  
  # tmp.res.cor %>%
  #   filter(abs(PCC) > 0.5) %>%
  #   nrow()
  
  # ggplot(tmp.res.cor,aes(Rank,PCC))+
  #   geom_point()+
  #   theme_cowplot(font_size = 28)
  # ggsave(filename = myFileName(prefix =  "res/fig/test_TimeTalk_PCC_rank_plot",
  #                              suffix = ".png"),
  #        width = 8,height = 8,dpi = 350)
  
  
  #### read master regulator analysis
  tmp.mra.res <- readRDS(file = "res/R/B_blastoid_RTN_mra_result.rds")
  tmp.TF.gene <- tmp.mra.res %>% 
    filter(group == tmp.ident.2) %>%
    pull(Regulon)
  # #### Let's show the TF levels
  # ph <- pheatmap_fixed(inter.tmp.mat.2$scaledData[tmp.TF.gene,],
  #                      show_colnames = F,
  #                      show_rownames = T,
  #                      cluster_cols = F,
  #                      color = coldwarm(100),
  #                      name = "scaled values")
  # ph
  
  
  tmp.LR.list <- tmp.res.cor %>%
    filter(abs(SCC) > 0.2) %>%
    pull(LRpairs)
  
  plan("multisession",workers = 10)
  tmp.ttt.res <- future_lapply(tmp.LR.list, function(tmp.LR){
    cat(tmp.LR,sep = "\n")
    tmp.L.gene <- unlist(lapply(strsplit(tmp.LR,split = "-"),FUN = function(ii){ ii[1]}))
    tmp.R.gene <- unlist(lapply(strsplit(tmp.LR,split = "-"),FUN = function(ii){ ii[2]}))
    
    x <- inter.tmp.mat.1[tmp.L.gene,]
    y <- inter.tmp.mat.2[tmp.R.gene,]
    
    tmp.res.list <- lapply(tmp.TF.gene, function(tmp.TF.gene.use){
      cat(paste0(tmp.ident.2,"_TF:",tmp.TF.gene.use),sep = "\n")
      tmp.TF.level <- inter.tmp.mat.2[tmp.TF.gene.use,]
      tmp.IS <- sqrt(x*y)
      tmp.res.LRtoTF.pvalue <- tryCatch(
        expr = {
          tmp.res <- grangertest(tmp.IS,tmp.TF.level)
          tmp.res$`Pr(>F)`[2]},
        error = function(e) {
          1
        })
      
      tmp.res.TFtoLR.pvalue <- tryCatch(
        expr = {
          tmp.res <- grangertest(tmp.TF.level,tmp.IS)
          tmp.res$`Pr(>F)`[2]},
        error = function(e) {
          1
        })
      tmp.PCC <- cor(tmp.IS,tmp.TF.level,method = "pearson")
      tmp.PCC <- ifelse(is.na(tmp.PCC),0,tmp.PCC)
      tmp.SCC <- cor(tmp.IS,tmp.TF.level,method = "spearman")
      tmp.SCC <- ifelse(is.na(tmp.SCC),0,tmp.SCC)
      tmp.res.df <- data.frame(LR=tmp.LR,
                               L=tmp.L.gene,
                               R=tmp.R.gene,
                               TF=tmp.TF.gene.use,
                               PCC=tmp.PCC,
                               SCC=tmp.SCC,
                               LRtoTF=tmp.res.LRtoTF.pvalue,
                               TFtoLR=tmp.res.TFtoLR.pvalue,
                               cell.L=tmp.ident.1,
                               cell.R=tmp.ident.2,
                               stringsAsFactors = F)
      return(tmp.res.df)
    })
    tmp.res.df <- Reduce(rbind,tmp.res.list)
    tmp.res.df <- tmp.res.df %>%
      mutate(category=ifelse(LRtoTF < 1e-2 | TFtoLR < 1e-2,"PASS","SKIP"))
    return(tmp.res.df)
  })
  plan("sequential")
  tmp.ttt.res.df <- Reduce(rbind,tmp.ttt.res)
  tmp.res.df <- tmp.ttt.res.df 
  return(tmp.res.df)
})

TimeTalk.res <- Reduce(rbind,tmp.res.TimeTalk.list)
TimeTalk.res.EPS_blastoid.EPI_PE <- TimeTalk.res %>%
  filter(category == "PASS")
TimeTalk.res.EPS_blastoid.EPI_PE$orig.ident <- "EPS_blastoid"

length(unique(TimeTalk.res.EPS_blastoid.EPI_PE$LR))


TimeTalk.res_EPI_PE <- rbind(TimeTalk.res.blastocyst.EPI_PE,
                             TimeTalk.res.EPS_blastoid.EPI_PE)

saveRDS(TimeTalk.res_EPI_PE,
        file = myFileName(prefix = "res/R/TimeTalk_res_B_blastoid_EPI_PE",suffix = ".rds"))



####--------2.4 Brief analysis of TimeTalk result---------

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

tmp.plot.2 <- tmp.TimeTalk.res %>%
  filter(orig.ident == "EPS_blastoid") %>%
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

write.table(tmp.TimeTalk.res,
            file = myFileName(prefix = "res/txt/TimeTalk_res_blastoids_EPI_PE_raw",suffix = ".txt"),
            quote = F,row.names = F,
            sep = "\t")

write.table(tmp.plot.1,
            file = myFileName(prefix = "res/txt/TimeTalk_res_blastoids_EPI_PE_filter_blastocyst",suffix = ".txt"),
            quote = F,row.names = F,
            sep = "\t")

write.table(tmp.plot.2,
            file = myFileName(prefix = "res/txt/TimeTalk_res_blastoids_EPI_PE_filter_blastoid",suffix = ".txt"),
            quote = F,row.names = F,
            sep = "\t")


boxplot(log10(tmp.plot.1$LRtoTF.ens))

tmp.list.plot <- list(blastocyst = tmp.plot.1$LR,
                      blastoid = tmp.plot.2$LR)

set.seed(666)

intersect(tmp.plot.1$LR,tmp.plot.2$LR)

png(filename = myFileName(prefix = "res/fig/fig5_B_blastoid_EPI-PE_LR",suffix = ".png"),
    width = 8,height = 8,res = 350,units = "in")
plot(eulerr::euler(tmp.list.plot,shape = "circle"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("#2369ad","#b9d9e9"),
     col="black",
     main=list(label="EPI-PE LR pairs overlap",fontsize = 24),
     bg="grey")
dev.off()

?eulerr::euler

length(tmp.list.plot$blastocyst)

#####---------2.5 show KEGG---------------
tmp.gene.list <- lapply(tmp.list.plot,
                        FUN = function(LRpairs){
  tmp.Lgene <- unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[1] }))
  tmp.Rgene <- unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[2] }))
  tmp.gene <- unique(c(tmp.Lgene,tmp.Rgene))
  eg = bitr(tmp.gene, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)
  tmp.gene <- eg$ENTREZID
  return(tmp.gene)
})


Compare <- compareCluster(geneCluster=tmp.gene.list,
                          fun="enrichKEGG",
                          organism = "mmu",
                          pvalueCutoff=0.05)

p <- dotplot(Compare,showCategory=20,font.size = 28)+
  scale_x_discrete(labels=c("blastocyst","EPS_blastoid"))+
  scale_color_gradientn(colors  = rev(coolwarm(100)))+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))

p
myggsave(p = p,prefix = "res/fig/fig5_TimeTalk_res",suffix = ".png",width = 16,height = 12,dpi=350)

###--------3.Rewrite the code by TimeTalk package-------------
####-------3.1 load data---------------------
tmp.mra.res <- readRDS(file = "res/R/B_blastoid_RTN_mra_result.rds")
LRpairs.df <- read.delim(file = "database/Ligand-Receptor-Pairs/Mouse/Mouse-2020-Shao-LR-pairs.txt",stringsAsFactors = F)

####load sc object
seu <- readRDS(file = "res/R/B_blastoid_seurat_inter_2_data_2022042218.rds")
cds <- readRDS(file = "res/R/B_blastoid_monocle3_cds_2022042219.rds")
####check cell types
seu$CellType <- Idents(seu)


DimPlot(seu)

TimeTalk.result <- RunTimeTalk(tmp.cds=cds,
                               tmp.seu=seu,
                               tmp.orig.ident = "blastocyst",
                               tmp.ident.1="EPI",
                               tmp.ident.2 = "PE",
                               LRpairs.df = LRpairs.df,
                               tmp.mra.res = tmp.mra.res,
                               tmp.winsz = 0.1,
                               numPts = 200,
                               tmp.SCC.cutoff = 0.2,
                               tmp.granger.cutoff = 1e-2)






