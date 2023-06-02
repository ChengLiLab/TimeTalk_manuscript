####TimeTalk manuscript
####Author:wlt
####Fig3
####20211003
####eLR validation
####edit this in 20211021
####20220318
####20220411: will modify this code to stable version

rm(list = ls())
####----0.1 load package--------
library(EnrichedHeatmap)
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
library(ggrepel)
library(VennDiagram)
#library(nichenetr)
library(igraph)
library(DiagrammeR)
library(DiagrammeRsvg)
library(patchwork)
library(GenomicFeatures)
library(ggbeeswarm)
source(file = "code/myUtils.R")
library(readxl)
#library(cellcall)
library(ggprism)
library(ggrepel)
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

####-----------1. visualization of eLR screen---------------------

####prepare data
tmp.color <- c("brown3","navy","#f2be58","grey")
names(tmp.color) <-  c("forward","backward","feedback","background")
data.plot <- readRDS(file = "res/R/putative_eLR_pairs_1094_2022032109.rds")
data.plot <- data.plot %>%
  filter(abs(PCC) > 0.1) %>%
  #filter(LRtoTF < 0.01 | TFtoLR < 0.01) %>%
  mutate(LRpairs = paste0(Lgene,"-",Rgene)) %>%
  mutate(group = "background") %>%
  mutate(group = ifelse(LRtoTF < 0.01 & TFtoLR > 0.01,"forward",group)) %>%
  mutate(group = ifelse(LRtoTF > 0.01 & TFtoLR < 0.01,"backward",group)) %>%
  mutate(group = ifelse(LRtoTF < 0.01 & TFtoLR < 0.01,"feedback",group)) %>%
  mutate(group = factor(group,levels = c("forward","backward","feedback","background"))) %>%
  mutate(LRtoTF=-log10(LRtoTF)) %>%
  mutate(TFtoLR=-log10(TFtoLR))


tmp.label <- c("Fgf4-Fgfr1","Fgf4-Fgfr2","Bmp6-Bmpr2",
               "Gdf9-Bmpr1b","Ntn1-Dcc",
               "Rspo2-Lgr5","Rspo2-Lgr6",'Hbegf-Erbb4',
               "Ctsd-Lrp1")


p <- ggplot(data.plot,aes(LRtoTF,TFtoLR,color=group,label=ifelse(LRpairs %in% tmp.label,LRpairs,"")))+
  geom_point(size=3,alpha=0.8)+
  geom_text_repel(box.padding = 0.5,
                  size=6, 
                  segment.curvature = 0.3,
                  segment.size  = 1,
                  force = 400,
                  fontface = "italic",
                  seed = 46,
                  arrow = arrow(length = unit(0.03, "npc")),
                  max.overlaps = Inf)+
  geom_hline(yintercept = 2,linetype="dotdash")+
  geom_vline(xintercept = 2,linetype="dotdash")+
  scale_color_manual(values = tmp.color,name=NULL,
                     guide=guide_legend(nrow = 2))+
  guides(color = guide_legend(override.aes = aes(label = ""),nrow = 2))+
  scale_x_continuous(breaks = seq(0,16,by=2))+
  scale_y_continuous(breaks = seq(0,16,by=2))+
  xlab(expression(-log[10](P[LRtoTF])))+
  ylab(expression(-log[10](P[TFtoLR])))+
  theme_cowplot(font_size = 28)+
  theme(axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        legend.position = "top",legend.justification = "center")
p  
myggsave(p = p,prefix = "res/fig/Fig3a_show_eLR_result",suffix = ".png",width = 8,height = 6,dpi = 350)
myggsave(p = p,prefix = "res/fig/Fig3a_show_eLR_result",suffix = ".pdf",width = 8,height = 6,dpi = 350)

colnames(data.plot)[8:9] <- c("-log10(LRtoTF_ens)","-log10(TFtoLR_ens)")

write.table(data.plot,
            file = myFileName(prefix = "res/txt/candidate_eLR",suffix = ".txt"),
            quote = F,
            row.names = F,
            sep = "\t")


#####--------2. nichenet----------
#####there are some problem here, so this section is ommit here



######---------2.1 load data------------
early.scRNAseq.exp <- readRDS(file = "res/R/early.scRNAseq.exp_20210114.rds")
early.scRNAseq.meta <- readRDS(file = "res/R/early.scRNAseq.meta_20210114.rds")
ligand_target_matrix = readRDS("res/R/nichenet_ligand_target_matrix.rds")
lr_network <- readRDS(file = "res/R/nichenet_lr_network_mouse.rds")
mouse_sig_network <- readRDS(file = "res/R/nichenet_mouse_sig_network.rds")
weighted_networks <- readRDS(file = "res/R/nichenet_mouse_weighted_network_20210131.rds")
ligand_tf_matrix = readRDS(file = "res/R/nichenet_ligand_tf_matrix.rds")

#####---------3.2.2 define expressed genes---------
tmp.id <- early.scRNAseq.meta %>%
  dplyr::filter(Stage == "early2cell") %>%
  pull(cell_id)
expressed_genes_sender <- early.scRNAseq.exp[,tmp.id] %>%
  apply(2,function(x){2^x-1}) %>% 
  apply(1,function(x){log2( mean(x) +1)}) %>% 
  .[. >= 4] %>% 
  names()
expressed_genes_receiver <- expressed_genes_sender



######-------3.2.3 define gene sets, ligands------
zga_gene_zhangyi <- read_delim(file = "database/other_ZGA/Nfya_KD_zhangyi.txt",delim = "\t") %>%
  pull(gene) 
geneset_oi <- zga_gene_zhangyi
background_expressed_genes = expressed_genes_receiver  %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
expressed_ligands = intersect(ligands,expressed_genes_sender)

receptors = lr_network %>% pull(to) %>% unique()
expressed_receptors = intersect(receptors,expressed_genes_receiver)

lr_network_expressed = lr_network %>% 
  dplyr::filter(from %in% expressed_ligands & to %in% expressed_receptors) 
head(lr_network_expressed)

potential_ligands = lr_network_expressed %>% pull(from) %>% unique()

########--------3.2.4 Perform NicheNet’s ligand activity analysis on the gene set of interest--------
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


######------------3.2.5 Infer target genes of top-ranked ligands and visualize in a heatmap ----------

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
                      "ZGA gene", 
                      color = "purple",
                      x_axis = T,
                      legend_position = "top", 
                      x_axis_position = "top",
                      legend_title = "Regulatory potential") + 
  scale_fill_gradientn(colours = hic.orage(256),
                        breaks = c(0,0.009),
                        labels = c(0,0.009))+ 
  theme(text = element_text(size = 24),
        axis.text.x.top= element_blank(),
        axis.ticks.x = element_blank())

p_ligand_target_network
round(range(vis_ligand_target),2)
ggsave(filename = myFileName(prefix = "res/fig/Fig3_ZGA_ligand_nichenet",suffix = ".png"),
       width = 24,height = 8,dpi = 350)



#####-----------3.2.6 : Ligand-receptor network inference for top-ranked ligands------------------------


# get the ligand-receptor network of the top-ranked ligands
lr_network_top = lr_network %>% 
  dplyr::filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

# get the weights of the ligand-receptor interactions as used in the NicheNet model
lr_network_top_df = weighted_networks$lr_sig %>% 
  dplyr::filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

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
p_ligand_pearson = vis_ligand_pearson %>% make_heatmap_ggplot("Prioritized CAF-ligands","Ligand activity", color = "darkorange",legend_position = "top", 
                                                              x_axis_position = "top", legend_title = "activity")
p_ligand_pearson = p_ligand_pearson+scale_fill_gradientn(colors  = blues(256),
                                      breaks = round(range(vis_ligand_pearson),2),
                                      labels = round(range(vis_ligand_pearson),2))

#?make_heatmap_ggplot

# (p_ligand_pearson | p_ligand_target_network) / p_ligand_receptor_network
# ggsave(filename = "res/fig/nichenet_res_test_early2cell_20210130.jpg",width = 8,height = 8,dpi = 350)

plot_grid(p_ligand_pearson+ylab("prioritized Ligands")+theme(text = element_text(size = 30),
                                                             legend.text = element_text(size = 24),
                                                             legend.text.align = 0.5,
                                                             axis.text.x.top = element_blank(),
                                                             axis.ticks.x = element_blank()),
          p_ligand_target_network+ylab(NULL)+
            theme(axis.text.y = element_blank(),
                  legend.text = element_text(size = 24),
                  axis.ticks.y = element_blank()),rel_widths = c(2,8))
ggsave(filename = myFileName(prefix = "res/fig/fig3_nichenet_res_ZGA",suffix = ".png"),
       width = 16,height = 9,dpi = 350,units = "in")

tmp.mat <- vis_ligand_pearson
ComplexHeatmap::pheatmap(mat = tmp.mat[rev(1:nrow(tmp.mat)),],
                         color = rdbu(100),
                         cluster_rows = F,cluster_cols = F)


######----------3.2.7 get graph------
ligands_all = c("Tgfb2") # this can be a list of multiple ligands if required
targets_all = c("Cdkn1a")

ligand_target_matrix = readRDS("res/R/nichenet_ligand_target_matrix.rds")
lr_network <- readRDS(file = "res/R/nichenet_lr_network_mouse.rds")
mouse_sig_network <- readRDS(file = "res/R/nichenet_mouse_sig_network.rds")
weighted_networks <- readRDS(file = "res/R/nichenet_mouse_weighted_network_20210131.rds")
ligand_tf_matrix = readRDS(file = "res/R/nichenet_ligand_tf_matrix.rds")
head(lr_network)
# weighted_networks = construct_weighted_networks(lr_network, mouse_sig_network, mouse_gr_network,
#                                                 source_weights_df)

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
jpeg(filename = myFileName(prefix = "res/fig/fig3_test_nichenet_path_width",suffix = ".jpg"),
     width = 8,height = 8,res = 350,units = "in")
plot.igraph(test,
            layout=tmp.layout,
            vertex.size = 30,
            vertex.label.font=18,
            vertex.label.cex=1,
            vertex.frame.color=tmp.vertext$fillcolor,
            vertex.color=tmp.vertext$fillcolor,
            vertex.label.color=tmp.vertext$fontcolor,
            edge.width = 10*tmp.edge$penwidth)
dev.off()


pdf(file  = myFileName(prefix = "res/fig/fig3_test_nichenet_path_width",suffix = ".pdf"),
    width = 8,height = 8)
plot.igraph(test,
            layout=tmp.layout,
            vertex.size = 30,
            vertex.label.font=18,
            vertex.label.cex=1,
            vertex.frame.color=tmp.vertext$fillcolor,
            vertex.color=tmp.vertext$fillcolor,
            vertex.label.color=tmp.vertext$fontcolor,
            edge.width = 10*tmp.edge$penwidth)
dev.off()


###-------4.gene sets overlap---------
mm9KG_txdb <- makeTxDbFromGFF(file = "D://Ubuntu/wlt/igenomes/Mus_musculus/UCSC/mm9/Annotation/genes.gtf")
mm9.gene <- genes(mm9KG_txdb)
mm9.gene.list <- mm9.gene$gene_id
test <- transcripts(mm9KG_txdb)
test <- exons(mm9KG_txdb)

columns(mm9KG_txdb)
library(TxDb.Mmusculus.UCSC.mm9.knownGene)
columns(TxDb.Mmusculus.UCSC.mm9.knownGene)



seqinfo(mm9KG_txdb)
seqlengths(mm9KG_txdb)



###-------4.1 essential gene---------

###-------4.1.1 LR--------

### load essential gene list
tmp.dir <- "database/deg_annotation_e.csv/deg_annotation_e.csv"
tmp.df <- read.delim(file = tmp.dir,sep = ";",stringsAsFactors = F,header = F)
essential.gene <- tmp.df %>%
  dplyr::filter(V8 == "Mus musculus") %>%
  pull(V3)
head(essential.gene)
essential.gene <- intersect(essential.gene,mm9.gene.list)

### LR gene
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_2022031117.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)
LRgene.all <- LRgene

#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(essential.gene)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,essential.gene))
n <- length(LRgene)
x <- M/N
y <- k/n
p <- phyper(k-1,M,N-M,n,lower.tail = F)
tmp.df.1 <- data.frame(value = c(x,y),group=c("all","LR"),stringsAsFactors = F)
pvalue.list <- p

####-------4.1.2 eLR------------------
eLR.df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds") %>%
  rownames_to_column("LRpairs") %>%
  mutate(Lgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[2]})))

Lgenelist <- eLR.df$Lgene 
Rgenelist <- eLR.df$Rgene
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)
eLRgene <- LRgene
#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(essential.gene)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,essential.gene))
n <- length(LRgene)
x <- M/N
y <- k/n


pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(LRgene.all,essential.gene)),
                                                  N = length(LRgene.all),
                                                  k = length(intersect(eLRgene,essential.gene)),
                                                  n = length(eLRgene)))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(mm9.gene.list,essential.gene)),
                                                  N = length(mm9.gene.list),
                                                  k = length(intersect(eLRgene,essential.gene)),
                                                  n = length(eLRgene)))


tmp.df.2 <- data.frame(value = y,
                       group=c("eLR"),
                       stringsAsFactors = F)
tmp.df <- rbind(tmp.df.1,tmp.df.2)
tmp.df$group <- factor(tmp.df$group,levels = c("all","LR","eLR"))

ggplot(tmp.df,aes(group,value,fill=group))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(value,2),vjust=-0.5),size=12)+
  geom_signif(comparisons = list(c("all","LR"),
                                 c("LR","eLR"),
                                 c("all","eLR")),
              annotations = paste0("p=",signif(pvalue.list,4)),
              y_position = c(0.4,0.55,0.65),size = 2,textsize = 10,vjust = - 0.5)+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.75),breaks = seq(0,0.7,by=0.1))+
  labs(x = NULL,y="ratio")+
  ggtitle("Essential gene enrichment")+
  scale_fill_manual(values = c("darkgrey",blues(3)[2:3]))+
  theme_cowplot(font_size = 40)+
  theme(legend.position = "none",
        plot.title = element_text(size = 36,hjust = 0.5),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2))
ggsave(filename = myFileName(prefix = "res/fig/Fig3_essential_gene_enrichment",
                             suffix = ".png"),
       width = 8,height = 8,dpi = 350)

ggsave(filename = myFileName(prefix = "res/fig/Fig3_essential_gene_enrichment",
                             suffix = ".pdf"),
       width = 8,height = 8,dpi = 350)

#####---------4.1.3 euler plot---------
tmp.list <- list(all = essential.gene,
                 LR = LRgene.all,
                 eLR = eLRgene)
names(tmp.list)[1] <- "Essential gene"
png(filename = myFileName(prefix = "res/fig/Fig3_essential_gene_euler",suffix = ".png"),
    width = 8,height = 8,units = "in",bg = "grey",res = 350)
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="Essential gene overlap",fontsize = 24),
     bg="grey")
dev.off()

pdf(file = myFileName(prefix = "res/fig/Fig3_essential_gene_euler",suffix = ".pdf"),
    width = 8,height = 8,bg = "grey")
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="Essential gene overlap",fontsize = 24),
     bg="grey")
dev.off()

###-------4.2 housekeeping gene---------

###-------4.2.1 LR--------

### load HK gene
tmp.dir <- "database/HRT_housekeeping_gene_database/Housekeeping_GenesMouse.csv"
tmp.df <- read.delim(file = tmp.dir,sep = ";",stringsAsFactors = F,header = T)
HK.gene <- tmp.df %>%
  pull(Gene) %>%
  unique()
HK.gene <- intersect(HK.gene,mm9.gene.list)

### LR gene
### LR gene
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_2022031117.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)
LRgene.all <- LRgene

#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(HK.gene)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,HK.gene))
n <- length(LRgene)
x <- M/N
y <- k/n
p <- phyper(k-1,M,N-M,n,lower.tail = T)
tmp.df.1 <- data.frame(value = c(x,y),group=c("all","LR"),stringsAsFactors = F)
pvalue.list <- p

####-------4.2.2 eLR------------------
eLR.df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds") %>%
  rownames_to_column("LRpairs") %>%
  mutate(Lgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[2]})))
Lgenelist <- eLR.df$Lgene 
Rgenelist <- eLR.df$Rgene

Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)
eLRgene <- LRgene
#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(HK.gene)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,HK.gene))
n <- length(LRgene)
x <- M/N
y <- k/n

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(LRgene.all,HK.gene)),
                                                  N = length(LRgene.all),
                                                  k = length(intersect(eLRgene,HK.gene)),
                                                  n = length(eLRgene)))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(mm9.gene.list,HK.gene)),
                                                  N = length(mm9.gene.list),
                                                  k = length(intersect(eLRgene,HK.gene)),
                                                  n = length(eLRgene),lower.tail = T))


tmp.df.2 <- data.frame(value = y,
                       group=c("eLR"),
                       stringsAsFactors = F)
tmp.df <- rbind(tmp.df.1,tmp.df.2)
tmp.df$group <- factor(tmp.df$group,levels = c("all","LR","eLR"))

ggplot(tmp.df,aes(group,value,fill=group))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(value,2),vjust=-0.5),size=12)+
  geom_signif(comparisons = list(c("all","LR"),
                                 c("LR","eLR"),
                                 c("all","eLR")),
              annotations = paste0("p=",signif(pvalue.list,4)),
              y_position = c(0.16,0.05,0.2),
              size = 2,textsize = 10,vjust = -0.5)+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.25),breaks = seq(0,0.2,by=.1))+
  labs(x = NULL,y="ratio")+
  ggtitle("Housekeeping gene \n enrichment")+
  scale_fill_manual(values = c("darkgrey",blues(3)[2:3]))+
  theme_cowplot(font_size = 40)+
  theme(legend.position = "none",
        plot.title = element_text(size = 36,hjust = 0.5),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2))
ggsave(filename = myFileName(prefix = "res/fig/Fig3_HK_gene_enrichment",
                             suffix = ".png"),
       width = 8,height = 8,dpi = 350)

ggsave(filename = myFileName(prefix = "res/fig/Fig3_HK_gene_enrichment",
                             suffix = ".pdf"),
       width = 8,height = 8,dpi = 350)
####--------4.2.3 euler plot---------------
tmp.list <- list(all = HK.gene,
                 LR = LRgene.all,
                 eLR = eLRgene)

names(tmp.list)[1] <- "Housekeeping gene"

png(filename = myFileName(prefix = "res/fig/Fig3_HK_gene_euler",suffix = ".png"),
    width = 8,height = 8,units = "in",bg = "grey",res = 350)
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="Housekeeping gene overlap",fontsize = 24),
     bg="grey")
dev.off()

pdf(file = myFileName(prefix = "res/fig/Fig3_HK_gene_euler",suffix = ".pdf"),
    width = 8,height = 8,bg = "grey")
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="Housekeeping gene overlap",fontsize = 24),
     bg="grey")
dev.off()




####--------4.3 ZGA gene overlap--------------

###-------4.3.1 LR--------

### load ZGA gene
zga_gene_zhangyi <- read_delim(file = "database/other_ZGA/Nfya_KD_zhangyi.txt",delim = "\t") %>%
  pull(gene)
zga_gene_zhangyi <- intersect(zga_gene_zhangyi,mm9.gene.list)

### LR gene
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_2022031117.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)
LRgene.all <- LRgene

#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(zga_gene_zhangyi)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,zga_gene_zhangyi))
n <- length(LRgene)
x <- M/N
y <- k/n
p <- phyper(k-1,M,N-M,n,lower.tail = F)
tmp.df.1 <- data.frame(value = c(x,y),group=c("all","LR"),stringsAsFactors = F)
pvalue.list <- p

####-------4.3.2 eLR------------------
eLR.df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds") %>%
  rownames_to_column("LRpairs") %>%
  mutate(Lgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[2]})))
Lgenelist <- eLR.df$Lgene 
Rgenelist <- eLR.df$Rgene
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)
eLRgene <- LRgene
#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(zga_gene_zhangyi)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,zga_gene_zhangyi))
n <- length(LRgene)
x <- M/N
y <- k/n

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(LRgene.all,zga_gene_zhangyi)),
                                                  N = length(LRgene.all),
                                                  k = length(intersect(eLRgene,zga_gene_zhangyi)),
                                                  n = length(eLRgene)))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(mm9.gene.list,zga_gene_zhangyi)),
                                                  N = length(mm9.gene.list),
                                                  k = length(intersect(eLRgene,zga_gene_zhangyi)),
                                                  n = length(eLRgene)))


tmp.df.2 <- data.frame(value = y,
                       group=c("eLR"),
                       stringsAsFactors = F)
tmp.df <- rbind(tmp.df.1,tmp.df.2)
tmp.df$group <- factor(tmp.df$group,levels = c("all","LR","eLR"))

ggplot(tmp.df,aes(group,value,fill=group))+
  geom_bar(stat="identity")+
  geom_text(aes(label=round(value,2),vjust=-0.5),size=12)+
  geom_signif(comparisons = list(c("all","LR"),
                                 c("LR","eLR"),
                                 c("all","eLR")),
              annotations = paste0("p=",signif(pvalue.list,4)),
              y_position = c(0.11,0.15,0.18),
              size = 2,textsize = 10,vjust = -0.5)+
  scale_y_continuous(expand = c(0,0),limits = c(0,0.25),breaks = seq(0,0.25,by=0.1))+
  labs(x = NULL,y="ratio")+
  ggtitle("ZGA gene  enrichment")+
  scale_fill_manual(values = c("darkgrey",blues(3)[2:3]))+
  theme_cowplot(font_size = 40)+
  theme(legend.position = "none",
        plot.title = element_text(size = 36,hjust = 0.5),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2))
ggsave(filename = myFileName(prefix = "res/fig/Fig3_ZGA_gene_enrichment",
                             suffix = ".png"),
       width = 8,height = 8,dpi = 350)
ggsave(filename = myFileName(prefix = "res/fig/Fig3_ZGA_gene_enrichment",
                             suffix = ".pdf"),
       width = 8,height = 8,dpi = 350)

####--------4.3.3 euler plot---------------
tmp.list <- list(all = zga_gene_zhangyi,
                 LR = LRgene.all,
                 eLR = eLRgene)
names(tmp.list)[1] <- "ZGA gene"
png(filename = myFileName(prefix = "res/fig/Fig3_ZGA_gene_euler",suffix = ".png"),
    width = 8,height = 8,units = "in",bg = "grey",res = 350)
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="ZGA gene overlap",fontsize = 24),
     bg="grey")
dev.off()


pdf(file = myFileName(prefix = "res/fig/Fig3_ZGA_gene_euler",suffix = ".pdf"),
    width = 8,height = 8,bg = "grey")
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="ZGA gene overlap",fontsize = 24),
     bg="grey")
dev.off()


# ####-----(drpreacted)4.3 gene age---------
# ####20211021
# ####The following code will remove in the future version
# 
# tmp.path <- "database/essential_gene_chart.xlsx"
# readxl::excel_sheets(tmp.path)
# tmp.df <- readxl::read_xlsx(tmp.path,sheet = "Sheet1",skip = 1)
# 
# colnames(tmp.hom)
# tmp.hom <- read.delim(file = "database/HOM_MouseHumanSequence.rpt",stringsAsFactors = F)
# tmp.hg <- tmp.hom %>%
#   dplyr::filter( Common.Organism.Name == "human" ) %>%
#   dplyr::filter( Symbol %in% tmp.df$gene) %>%
#   dplyr::select( HomoloGene.ID,Symbol)
# tmp.mm <- tmp.hom %>%
#   dplyr::filter( Common.Organism.Name == "mouse, laboratory" ) %>%
#   dplyr::filter( HomoloGene.ID %in% tmp.hg$HomoloGene.ID) %>%
#   dplyr::select( HomoloGene.ID,Symbol)
# tmp.hom <- merge(tmp.hg,tmp.mm,by="HomoloGene.ID") 
# colnames(tmp.hom)[2:3] <- c("human","mouse")
# 
# colnames(tmp.df)
# 
# not_eLR_gene <- setdiff(LRgene.all,eLRgene)
# #### obtain gene_age table
# gene_age.df <- tmp.df %>%
#   dplyr::select(gene,`gene_age(26)`) %>%
#   dplyr::filter(gene %in% tmp.hom$human) %>%
#   mutate(gene=plyr::mapvalues(gene,from = tmp.hom$human,to = tmp.hom$mouse)) %>%
#   mutate(group="not_LR") %>%
#   mutate(group=ifelse(gene %in% eLRgene,"eLR",group)) %>%
#   mutate(group=ifelse(gene %in% not_eLR_gene,"other_LR",group))
# table(gene_age.df$group)
# gene_age.df$`gene_age(26)` <- as.numeric(gene_age.df$`gene_age(26)`)
# ggplot(gene_age.df,aes(group,`gene_age(26)`,fill=group))+
#   geom_boxplot(size=2)+
#   theme_cowplot(font_size = 30)+
#   xlab(NULL)+
#   theme(axis.line = element_line(size = 2),
#         legend.position = "none",
#         axis.ticks = element_line(size = 2))
# ggsave(filename = myFileName(prefix = "res/fig/Fig3_gene_age_LRgene",
#                              suffix = ".png"),width = 8,height = 8,dpi = 350)





####-------5.blastoid analysis------------

####-------5.1 load data and correct bath-------

####-------5.1.1 load data------------
####EPS blastoid
tmp.dir <- "data/blastoid/GSE135701/EPS_blastoid/"
tmp.mat <- Read10X(data.dir = tmp.dir)
tttt <- DropletUtils::barcodeRanks(tmp.mat)
br.out <- tttt
tmp.plot <- br.out[!is.na(br.out$fitted),]
tmp.cell.id <- rownames(tmp.plot)
tmp.mat.1 <- tmp.mat[,tmp.cell.id]
seu <- CreateSeuratObject(tmp.mat.1,min.features = 2000,min.cells = 0,project = "EPS_blastoid")
colnames(seu@meta.data)
VlnPlot(seu,features = c("nCount_RNA","nFeature_RNA"))
seu.1 <- seu

####blastocyst
tmp.dir <- "data/blastoid/GSE135701/Blastocyst/"
tmp.mat <- Read10X(data.dir = tmp.dir)
tttt <- DropletUtils::barcodeRanks(tmp.mat)
br.out <- tttt
tmp.plot <- br.out[!is.na(br.out$fitted),]
tmp.cell.id <- rownames(tmp.plot)
tmp.mat.1 <- tmp.mat[,tmp.cell.id]
seu <- CreateSeuratObject(tmp.mat.1,min.features = 2000,min.cells = 0,project = "blastocyst")
colnames(seu@meta.data)
seu.2 <- seu

###-----5.1.2 direct merge----------
seu <- merge(seu.1,seu.2)
seu <- NormalizeData(seu)
seu <- ScaleData(seu)
seu <- FindVariableFeatures(seu,selection.method = "vst",nfeatures = 3000)
seu <- RunPCA(seu,features = VariableFeatures(seu))
#### The following code doesn't work small samples
# tmp.mat <- as.matrix(seu@assays$RNA@data)
# tmp.genes <- VariableFeatures(seu)
# seu.dist <- 1- cor(tmp.mat[tmp.genes,],method = "pearson")
# seu <- RunTSNE(seu,
#                distance.matrix = seu.dist,
#                dims = 1:20,
#                seed.use = 42,
#                perplixity=100)
seu <- RunTSNE(seu,dims = 1:20,seed.use = 42)
seu <- RunUMAP(seu,dims = 1:20,umap.method = "umap-learn")
seu <- FindNeighbors(seu,dims = c(1:20))
seu <- FindClusters(seu,resolution = 0.2)
UMAPPlot(seu,group.by="orig.ident")
colnames(seu@meta.data)
FeaturePlot(seu,features = c("Ascl2"))


####------5.1.3 batch effect correct-------------
tmp.list <- list(seu.1,seu.2)
names(tmp.list) <- c("EPS_blastoid","blastocyst")
tmp.list <- lapply(names(tmp.list),FUN = function(x){
  tmp.list[[x]] <- NormalizeData(tmp.list[[x]], verbose = FALSE)
  tmp.list[[x]] <- FindVariableFeatures(tmp.list[[x]], selection.method = "vst", 
                                             nfeatures = 2000, verbose = FALSE)
  return(tmp.list[[x]])
})

#reference.list <- pancreas.list[c("PT081", "PT089")]
#?FindIntegrationAnchors
tmp.list <- FindIntegrationAnchors(object.list = tmp.list, 
                                   dims = 1:20,
                                   k.anchor = 5,
                                   k.filter = 30)
#?IntegrateData
seu.integrated <- IntegrateData(anchorset = tmp.list, dims = 1:20)


###RunSeurat pipeline
DefaultAssay(seu.integrated) <- "integrated"
seu.integrated <- ScaleData(seu.integrated, verbose = FALSE)
seu.integrated <- RunPCA(seu.integrated, npcs = 30, verbose = FALSE)
seu.integrated <- RunUMAP(seu.integrated, reduction = "pca", dims = 1:30,umap.method = "umap-learn")
seu.integrated <- FindNeighbors(seu.integrated,dims = c(1:20))
seu.integrated <- FindClusters(seu.integrated,resolution = 0.2)

####-------5.2 assign clusters-------------

####-------5.2.1 assign markers----------
TE.markers <- c("Cdx2","Krt8","Krt18","Ascl2","Tacstd2")
ICM.markers <- c("Pou5f1","Nanog","Sox2","Esrrb","Sox15")
PE.markers <- c("Gata4","Gata6","Sox17","Pdgfra","Col4a1")
p1 <- UMAPPlot(seu.integrated) 
p1
tmp.markers <- c(TE.markers,
                 "Gapdh",ICM.markers,
                 "Ppia",PE.markers)

# p2 <- FeaturePlot(seu.integrated,features = c(TE.markers,
#                                    "Gapdh",ICM.markers,
#                                    "Ppia",PE.markers),combine = F) 

p2 <- lapply(tmp.markers, function(x){
  res.p <- FeaturePlot(seu.integrated,features = x)+
    scale_color_gradientn(colours = rdwhbu(256))
  return(res.p)
})

p <- wrap_plots(c(list(p1),p2),nrow = 3) & 
  theme_cowplot(font_size = 24) &
  NoAxes()

ggsave(plot = p,filename = myFileName(prefix = "res/fig/Fig3_EPS_blastoid_markers_raw",
                             suffix = ".png"),
       dpi = 350,width = 18,height = 8)

seu.integrated <- RenameIdents(seu.integrated,
                               "0"="intermediate",
                               "1"="PE",
                               "2"="intermediate",
                               "3"="PE",
                               "4"="intermediate",
                               "5"="EPI",
                               "6"="TE",
                               "7"="EPI")
levels(seu.integrated) <- c("EPI","TE","PE","intermediate")
p1 <- UMAPPlot(seu.integrated) + 
  scale_color_manual(values = c("#e90b8e","#3b53a5","#6ebe45","grey"))
p1
p <- wrap_plots(c(list(p1),p2),nrow = 3) & 
  theme_cowplot(font_size = 24) &
  NoAxes()
p
ggsave(plot = p,filename = myFileName(prefix = "res/fig/Fig3_EPS_blastoid_markers",
                                      suffix = ".png"),
       dpi = 350,width = 18,height = 8)

p3 <- UMAPPlot(seu.integrated,pt.size=1.5,
               split.by="orig.ident") +
  scale_color_manual(values = c("#e90b8e","#3b53a5","#6ebe45","grey"))+
  theme(axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2))
p3
ggsave(filename = myFileName(prefix = "res/fig/Fig3_EPS_split_plot",
                             suffix = ".png"),
       width = 8,height = 4,dpi = 350)
saveRDS(seu.integrated,file = myFileName(prefix = "res/R/EPS_blastoid",suffix = ".rds"))


####-------------5.2.2 show markers-------------------
seu.integrated <- readRDS(file = "res/R/EPS_blastoid_20211017.rds")
TE.markers <- c("Cdx2","Krt8","Krt18","Ascl2","Tacstd2")
ICM.markers <- c("Pou5f1","Nanog","Sox2","Esrrb","Sox15")
PE.markers <- c("Gata4","Gata6","Sox17","Pdgfra","Col4a1")
tmp.markers <- c("Gapdh","Ppia",
                 TE.markers,
                 ICM.markers,
                 PE.markers)
tmp.seu <- seu.integrated
levels(tmp.seu) <- c("TE","EPI","PE","intermediate")
tmp.seu$cell.type <- Idents(tmp.seu)
DefaultAssay(tmp.seu) <- "RNA"
head(tmp.seu[[]])
p1 <- DimPlot(tmp.seu,group.by = "orig.ident")
p2 <- DimPlot(tmp.seu,group.by = "cell.type")
p1 | p2
UMAPPlot(tmp.seu)
myscAreaPlot(tmp.seu = tmp.seu,
             features = tmp.markers,
             tmp.color = divergentcolor(length(levels(tmp.seu))))





###------------5.3 IS -----------

###---------5.3.1 define function-------
MyAverageSeurat <- function(seu.obj){
  cat("Calcualte Average Expression",sep = "\n")
  tmp.exp.data  <-  GetAssayData(seu.obj,slot = "data",assay = "RNA")
  tmp.exp.meta  <- seu.obj[[]]
  tmp.cell.type <- levels(seu.obj)
  tmp.res <- lapply(tmp.cell.type,FUN = function(x){
    cat(x,sep = "\n")
    tmp.id <- tmp.exp.meta %>% 
      dplyr::filter(cell_type == x) %>%
      pull(cell_id)
    tmp.exp.data <- rowMeans(expm1(tmp.exp.data[,tmp.id]))
    return(tmp.exp.data)
  })
  tmp.res <- Reduce("cbind",tmp.res)
  colnames(tmp.res) <- tmp.cell.type
  return(tmp.res)
}


myInteraction_Score <- function(gene.exp.mat=tmp.mat,
                                cell.type=tmp.cell.type,
                                Lgenelist=Lgenelist,
                                Rgenelist=Rgenelist){
  cat("Calcualte IS",sep = "\n")
  tmp.res.IS <- lapply(cell.type,FUN = function(x){
    cat(x,sep = "\n")
    res <- gene.exp.mat[Lgenelist,x] * gene.exp.mat[Rgenelist,]
    colnames(res) <- paste0(x,"_",cell.type)
    return(res)
  })
  tmp.res.IS <- Reduce(cbind,tmp.res.IS)
  rownames(tmp.res.IS) <- paste0(Lgenelist,"_",Rgenelist)
  return(tmp.res.IS)
}

myRandomIS <- function(seu.obj=seu,Lgenelist=Lgenelist,Rgenelist=Rgenelist){
  cat("random_sample",sep = "\n")
  tmp.exp.data  <-  GetAssayData(seu.obj,slot = "data",assay = "RNA")
  tmp.cell.type <- levels(seu.obj)
  tmp.label <- unname(seu.obj$cell_type)
  tmp.label <- sample(x = as.character(tmp.label),size = length(tmp.label),replace = F)
  seu.obj$random_label <- tmp.label
  tmp.exp.meta  <- seu.obj[[]]
  
  tmp.mat <- lapply(tmp.cell.type,FUN = function(x){
    cat(x,sep = "\n")
    tmp.id <- tmp.exp.meta %>% 
      dplyr::filter(random_label == x) %>%
      pull(cell_id)
    tmp.exp.data <- rowMeans(expm1(tmp.exp.data[,tmp.id]))
    return(tmp.exp.data)
  })
  tmp.mat <- Reduce("cbind",tmp.mat)
  colnames(tmp.mat) <- tmp.cell.type
  
  tmp.res <- myInteraction_Score(gene.exp.mat = tmp.mat,
                                 cell.type = colnames(tmp.mat),
                                 Lgenelist = Lgenelist,
                                 Rgenelist = Rgenelist)
  return(tmp.res)
}

myISwithPvalue <- function(seu=seu,Lgenelist = Lgenelist,Rgenelist = Rgenelist,calc.p=FALSE){
  
  tmp.seu <- seu
  tmp.res <- MyAverageSeurat(seu.obj = tmp.seu)
  tmp.res.IS <- myInteraction_Score(gene.exp.mat = tmp.res,
                                    cell.type = colnames(tmp.res),
                                    Lgenelist = Lgenelist,
                                    Rgenelist = Rgenelist)
  tmp.res.list <- list(IS=tmp.res.IS)
  
  if(calc.p){
    #### use future to do parallel
    options(future.globals.maxSize= 1024^3)
    plan(multisession, workers = 10)
    nboot = 1000
    tmp.res.IS.boot <- future_replicate(nboot, myRandomIS(seu.obj = seu,Lgenelist = Lgenelist,Rgenelist = Rgenelist))
    .tempGetPvalueFromOnePairs <- function(xx,tmp.res.IS,tmp.res.IS.boot,nboot=1000){
      test.p <- lapply(1:nrow(tmp.res.IS), 
                       FUN = function(x) sum(tmp.res.IS[x,xx] >= tmp.res.IS.boot[x,xx,])/nboot)
      test.p <- unlist(test.p)
      return(test.p)
    }
    GetPvalueFromISBoot <- function(tmp.res.IS,tmp.res.IS.boot){
      tmp.res <- lapply(1:ncol(tmp.res.IS),
                        FUN = function(x) .tempGetPvalueFromOnePairs(xx = x,tmp.res.IS = tmp.res.IS,tmp.res.IS.boot = tmp.res.IS.boot))
      tmp.res.df <- Reduce(cbind,tmp.res)
      return(tmp.res.df)
    }
    tmp.res.IS.Pvalue <- GetPvalueFromISBoot(tmp.res.IS = tmp.res.IS,tmp.res.IS.boot = tmp.res.IS.boot)
    rownames(tmp.res.IS.Pvalue) <- rownames(tmp.res.IS)
    colnames(tmp.res.IS.Pvalue) <- colnames(tmp.res.IS)
  }
  tmp.res.list <- c(tmp.res.list,list(IS.p=tmp.res.IS.Pvalue))
  return(tmp.res.list)
}


####--------5.3.2 load data----------

####load seurat
seu.integrated <- readRDS(file = "res/R/EPS_blastoid_20211017.rds")
seu.integrated$cell_type = Idents(seu.integrated)
seu.integrated$cell_id = colnames(seu.integrated)
seu.integrated <- subset(seu.integrated,idents = "intermediate", invert = TRUE)

seu.list <- SplitObject(seu.integrated,split.by = 'orig.ident')

###load pairs
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)



####filter genes
DefaultAssay(seu.integrated) <- "RNA"
gene_symbols <- rownames(seu.integrated)
l.remove <- setdiff(Lgenelist,gene_symbols)
r.remove <- setdiff(Rgenelist,gene_symbols)
index.remove <- c(which(Lgenelist %in% l.remove),which(Rgenelist %in% r.remove))
LRpairs <- LRpairs[-index.remove]
Lgenelist <- Lgenelist[-index.remove]
Rgenelist <- Rgenelist[-index.remove]
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)

####load cell.type
# tmp.cell.type <- c("EPI","TE","PE")
# 
# 
# head(colnames(seu.integrated))
# tmp.seu <- seu.list$EPS_blastoid
# tmp.res <- MyAverageSeurat(seu.obj = tmp.seu)
# tmp.res.IS <- myInteraction_Score(gene.exp.mat = tmp.res[,tmp.cell.type],
#                                   cell.type = tmp.cell.type,
#                                   Lgenelist = Lgenelist,
#                                   Rgenelist = Rgenelist)
# rownames(tmp.res.IS) <- paste0(Lgenelist,"_",Rgenelist)


tmp.seu <- seu.list$EPS_blastoid
cell_id <- sample(x = colnames(tmp.seu),size = ncol(seu.list$blastocyst))
tmp.seu <- subset(tmp.seu,cells=cell_id)
UMAPPlot(tmp.seu)
UMAPPlot(seu.list$blastocyst)

tmp.res.list <- myISwithPvalue(seu = tmp.seu,
                               Lgenelist = Lgenelist,
                               Rgenelist = Rgenelist,
                               calc.p=T)
tmp.res.list.EPS_blastoid <- tmp.res.list

tmp.res.list <- myISwithPvalue(seu = seu.list$blastocyst,
                               Lgenelist = Lgenelist,
                               Rgenelist = Rgenelist,
                               calc.p=T)
tmp.res.list.blastocyst <- tmp.res.list


saveRDS(object = tmp.res.list.blastocyst, file = "res/R/GSE135701_blastocyst.rds")
saveRDS(object = tmp.res.list.EPS_blastoid, file = "res/R/GSE135701_EPS_blastoid.rds")


####------5.3.3 buble plot--------
head(tmp.res.list$IS)
head(tmp.res.list$IS.p)

tmp.IS.p.buble.plot <- function(tmp.res.list,title=NULL){
  
  data.plot <- as.data.frame(tmp.res.list$IS)
  data.plot <- data.plot %>%
    rownames_to_column("LRpairs") %>%
    gather(key = "cell_pairs",value="IS",-LRpairs)
  
  data.plot.p <- as.data.frame(tmp.res.list$IS.p)
  data.plot.p <- data.plot.p %>%
    rownames_to_column("LRpairs") %>%
    gather(key = "cell_pairs",value="IS.p",-LRpairs)
  
  
  data.plot$IS.p <- data.plot.p$IS.p
  
  #### load eLR
  
  tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
  tmp.cutoff <- 0.1
  tmp.LR.pairs.1 <- tmp.res.cor.1 %>%
    dplyr::filter(PCC > tmp.cutoff) %>%
    head(20) %>%
    pull(LRpairs)
  
  tmp.LR.pairs.2 <- tmp.res.cor.1 %>%
    dplyr::filter(PCC < tmp.cutoff) %>%
    tail(20) %>%
    pull(LRpairs)
  
  tmp.LR.pairs <- c(tmp.LR.pairs.1,tmp.LR.pairs.2)
  
  tmp.LR.pairs <- gsub(pattern = "_",replacement = "-",tmp.LR.pairs)
  data.plot$LRpairs <- gsub(pattern = "_",replacement = "-",data.plot$LRpairs)
  data.plot$cell_pairs <- gsub(pattern = "_",replacement = "-",data.plot$cell_pairs)
  
  #tmp.LR.pairs <- c(tmp.LR.pairs,"Bmp4_Bmpr2")
  
  # tmp.LR.pairs.1
  # tmp.LR.pairs.2
  
  plot.data <- data.plot %>%
    filter(LRpairs %in% tmp.LR.pairs)
  plot.data$LRpairs <- factor(plot.data$LRpairs,levels = rev(tmp.LR.pairs))
  
  ###my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  tmp.min <- 0.001
  #log10(tmp.min)
  plot.data$IS.p[plot.data$IS.p <= tmp.min] = tmp.min
  
  tmp.min <- 0.000001
  plot.data$IS[plot.data$IS <= tmp.min] = tmp.min
  
  
  p <- ggplot(plot.data,aes(x=cell_pairs,y=LRpairs)) +
    geom_point(aes(size=-log10(IS.p),color=log10(IS))) +
    scale_color_gradientn('Log10(IS)', colors=rdbu(100)) +
    theme_bw(base_size = 24) +
    ggtitle(title)+
    guides(size = guide_legend(order = 1),
           color = guide_colorbar(order = 2))+
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text = element_text(size = 24, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1),
          axis.text.y = element_text(face = "italic"),
          axis.title=element_blank(),
          plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(size = 0.7, 
                                      linetype = "solid", 
                                      colour = "black"))
  
  return(p)
  
}

###plot dot
tmp.res.list.blastocyst <- readRDS(file = "res/R/GSE135701_blastocyst.rds")
p1 <- tmp.IS.p.buble.plot(tmp.res.list = tmp.res.list.blastocyst,
                          title = "blastocyst")
p1
ggsave(filename = myFileName(prefix = "res/fig/fig3_scInt_IS_blastocyst_buble",
                             suffix = ".png"),
       width = 8,height = 16,dpi = 350)

tmp.res.list.EPS_blastoid <- readRDS(file = "res/R/GSE135701_EPS_blastoid.rds")
p2 <- tmp.IS.p.buble.plot(tmp.res.list = tmp.res.list.EPS_blastoid,
                          title = "EPS blastoid")
p2
ggsave(filename = myFileName(prefix = "res/fig/fig3_scInt_IS_EP_buble",suffix = ".png"),
       width = 8,height = 16,dpi = 350)

p1+theme(legend.position = "none",
         axis.text.x = element_text(size = 18),
         axis.text.y = element_text(size = 18,face = "italic")) | 
  p2+theme(axis.text.x = element_text(size = 18),
           axis.text.y = element_blank(),
           axis.ticks.y = element_blank())

ggsave(filename = myFileName(prefix = "res/fig/fig3_scInt_IS",suffix = ".png"),
       width = 9,height = 12,dpi = 350,units = "in")


####------5.3.4 diff---------
tmp.res.list <- tmp.res.list.blastocyst
data.plot <- as.data.frame(tmp.res.list$IS)
data.plot <- data.plot %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "cell_pairs",value="IS",-LRpairs)

data.plot.p <- as.data.frame(tmp.res.list$IS.p)
data.plot.p <- data.plot.p %>%
  rownames_to_column("LRpairs") %>%
  gather(key = "cell_pairs",value="IS.p",-LRpairs)


data.plot$IS.p <- data.plot.p$IS.p

#### load eLR

tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
tmp.cutoff <- 0.1
tmp.LR.pairs.1 <- tmp.res.cor.1 %>%
  dplyr::filter(PCC > tmp.cutoff) %>%
  pull(LRpairs)

tmp.LR.pairs.2 <- tmp.res.cor.1 %>%
  dplyr::filter(PCC < tmp.cutoff) %>%
  pull(LRpairs)

tmp.LR.pairs <- c(tmp.LR.pairs.1,tmp.LR.pairs.2)

#tmp.LR.pairs <- c(tmp.LR.pairs,"Bmp4_Bmpr2")

# tmp.LR.pairs.1
# tmp.LR.pairs.2

plot.data <- data.plot %>%
  filter(LRpairs %in% tmp.LR.pairs)



tmp.plot <- log10(tmp.res.list.EPS_blastoid$IS+1) - log10(tmp.res.list.blastocyst$IS+1)
tmp.LR.pairs <- intersect(rownames(tmp.plot),tmp.LR.pairs)
tmp.plot <- tmp.plot[tmp.LR.pairs,]
tmp.sd <- -sort(-apply(tmp.plot, 1, sd))
tmp.sd <- tmp.sd[tmp.sd>0.1]

tmp.plot <- tmp.plot[names(tmp.sd),]
tmp.plot[tmp.plot>1] <- 1
tmp.plot[tmp.plot < -1] <- -1

png(filename = myFileName(prefix = "res/fig/fig3_scInt_diff",suffix = ".png"),
    width = 8,height = 16,units = "in",res = 350)
ComplexHeatmap::pheatmap(tmp.plot,
         color = rdbu(100),
         cluster_cols = F,
         show_rownames = T,fontsize = 20,
         heatmap_legend_param = list(title = expression(Delta(IS)),
                                     title_gp = gpar(col = "black", fontsize = 20),
                                     labels_gp = gpar(col = "black", fontsize = 20)),
         show_colnames = T)
dev.off()




#####----------5.3.5 diff pheatmap---------------
tmp.IS.bubble.data.get <- function(tmp.res.list,title=NULL){
  
  data.plot <- as.data.frame(tmp.res.list$IS)
  data.plot <- data.plot %>%
    rownames_to_column("LRpairs") %>%
    gather(key = "cell_pairs",value="IS",-LRpairs)
  
  data.plot.p <- as.data.frame(tmp.res.list$IS.p)
  data.plot.p <- data.plot.p %>%
    rownames_to_column("LRpairs") %>%
    gather(key = "cell_pairs",value="IS.p",-LRpairs)
  
  
  data.plot$IS.p <- data.plot.p$IS.p
  
  #### load eLR
  
  tmp.res.cor.1 <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
  tmp.cutoff <- 0.1
  tmp.LR.pairs.1 <- tmp.res.cor.1 %>%
    dplyr::filter(PCC > tmp.cutoff) %>%
    head(20) %>%
    pull(LRpairs)
  
  tmp.LR.pairs.2 <- tmp.res.cor.1 %>%
    dplyr::filter(PCC < tmp.cutoff) %>%
    tail(20) %>%
    pull(LRpairs)
  
  tmp.LR.pairs <- c(tmp.LR.pairs.1,tmp.LR.pairs.2)
  
  tmp.LR.pairs <- gsub(pattern = "_",replacement = "-",tmp.LR.pairs)
  data.plot$LRpairs <- gsub(pattern = "_",replacement = "-",data.plot$LRpairs)
  data.plot$cell_pairs <- gsub(pattern = "_",replacement = "-",data.plot$cell_pairs)
  
  #tmp.LR.pairs <- c(tmp.LR.pairs,"Bmp4_Bmpr2")
  
  # tmp.LR.pairs.1
  # tmp.LR.pairs.2
  
  plot.data <- data.plot %>%
    filter(LRpairs %in% tmp.LR.pairs)
  plot.data$LRpairs <- factor(plot.data$LRpairs,levels = rev(tmp.LR.pairs))
  
  ###my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  tmp.min <- 0.001
  #log10(tmp.min)
  plot.data$IS.p[plot.data$IS.p <= tmp.min] = tmp.min
  
  tmp.min <- 0.000001
  plot.data$IS[plot.data$IS <= tmp.min] = tmp.min
  
  return(plot.data)
  
}

p1 <- tmp.IS.p.buble.plot(tmp.res.list = tmp.res.list.blastocyst,title = "blastocyst")
ggsave(filename = myFileName(prefix = "res/fig/fig3_scInt_IS_blastocyst_buble",suffix = ".png"),
       width = 8,height = 16,dpi = 350)

p2 <- tmp.IS.p.buble.plot(tmp.res.list = tmp.res.list.EPS_blastoid,title = "EPS blastoid")
ggsave(filename = myFileName(prefix = "res/fig/fig3_scInt_IS_EP_buble",suffix = ".png"),
       width = 8,height = 16,dpi = 350)

tmp.data.get.1 <- tmp.IS.bubble.data.get(tmp.res.list = tmp.res.list.blastocyst)
tmp.data.get.2 <- tmp.IS.bubble.data.get(tmp.res.list = tmp.res.list.EPS_blastoid)
all(tmp.data.get.1$LRpairs==tmp.data.get.2$LRpairs)
tmp.plot <- tmp.data.get.1 %>%
  mutate(IS.diff=tmp.data.get.1$IS-tmp.data.get.2$IS) %>%
  dplyr::select(LRpairs,cell_pairs,IS.diff) %>%
  spread(key = cell_pairs,value = IS.diff) %>%
  column_to_rownames("LRpairs") 

tmp.plot <- tmp.plot[rev(1:nrow(tmp.plot)),]
tmp.plot[tmp.plot>1] <- 1
tmp.plot[tmp.plot < -1] <- -1

png(filename = myFileName(prefix = "res/fig/fig3_scInt_diff",suffix = ".png"),
    width = 8,height = 16,units = "in",res = 350)
ComplexHeatmap::pheatmap(tmp.plot,
                         color = rdbu(100),
                         cluster_cols = F,
                         cluster_rows = F,
                         show_rownames = T,
                         fontsize = 20,
                         heatmap_legend_param = list(title = expression(Delta(IS)),
                                                     title_gp = gpar(col = "black", fontsize = 20),
                                                     labels_gp = gpar(col = "black",fontface="italic",fontsize = 20)),
                         show_colnames = T)
dev.off()


#### The folllowing code wiil be removed in the near future
# ####----5.3.6 circos plot--------
# 
# 
# tmp.res.list <- tmp.res.list.blastocyst
# 
# tmp.circos.plot <- function(tmp.res.list){
#   data.plot <- as.data.frame(tmp.res.list$IS)
#   data.plot <- data.plot %>%
#     rownames_to_column("LRpairs") %>%
#     gather(key = "cell_pairs",value="IS",-LRpairs)
#   
#   data.plot.p <- as.data.frame(tmp.res.list$IS.p)
#   data.plot.p <- data.plot.p %>%
#     rownames_to_column("LRpairs") %>%
#     gather(key = "cell_pairs",value="IS.p",-LRpairs)
#   data.plot$IS.p <- data.plot.p$IS.p
#   
#   data.plot$Ligand <- unlist(lapply(strsplit(data.plot$LRpairs,split = "_"),FUN = function(x) x[1])) 
#   data.plot$Receptor <- unlist(lapply(strsplit(data.plot$LRpairs,split = "_"),FUN = function(x) x[2])) 
#   data.plot$cell_L <- unlist(lapply(strsplit(data.plot$cell_pairs,split = "_"),FUN = function(x) x[1])) 
#   data.plot$cell_R <- unlist(lapply(strsplit(data.plot$cell_pairs,split = "_"),FUN = function(x) x[2]))
#   
#   
#   data.plot.gene <- data.frame(cell=c(data.plot$cell_L,data.plot.circos$cell_R),
#                                gene=c(data.plot$Ligand,data.plot.circos$Receptor),
#                                lr=c(rep("ligand",length(data.plot$Ligand)),
#                                     rep("receptor",length(data.plot$Receptor))),
#                                stringsAsFactors = F) %>%
#     unique()
#   data.plot.gene$gene_id <- paste0("gene",1:nrow(data.plot.gene))
#   data.plot.gene$map_id <- paste0(data.plot.gene$cell,"_",data.plot.gene$gene,"_",data.plot.gene$lr)
#   
#   
#   
#   library(circlize)
#   df <- data.plot.gene %>%
#     arrange(cell)
#   
#   tmp.cell.type <- unique(df$cell)
#   tmp.cell.colors <- c("#e90b8e","#3b53a5","#6ebe45","grey")
#   names(tmp.cell.colors) <- c("EPI","TE","PE","intermediate")
#   
#   circos.clear()
#   circos.par(canvas.xlim =c(-1.1,1.1),
#              canvas.ylim = c(-1.1,1.1),
#              gap.degree=0.001,
#              cell.padding = c(0.02, 0, 0.02, 0))
#   fa = df$gene_id
#   fa = factor(fa,levels = fa)
#   circos.initialize(factors = fa, xlim = c(0,1))
#   
#   circos.trackPlotRegion(
#     
#     ylim = c(0, 1), track.height = 0.15, bg.border = NA,
#     
#     panel.fun = function(x, y) {
#       
#       sector.index = get.cell.meta.data("sector.index")
#       
#       xlim = get.cell.meta.data("xlim")
#       
#       ylim = get.cell.meta.data("ylim")
#       
#     } )
#   
#   
#   
#   lapply(1:length(tmp.cell.type),FUN = function(ii){
#     x <- tmp.cell.type[ii]
#     tmp.gene.id <- df %>%
#       dplyr::filter(cell==x) %>%
#       pull(gene_id)
#     highlight.sector(tmp.gene.id, 
#                      track.index = 1,
#                      text = x, 
#                      niceFacing = F, 
#                      font = 2, 
#                      col = tmp.cell.colors[x])
#   })
#   
#   circos.trackPlotRegion(
#     
#     ylim = c(0, 1), track.height = 0.08, bg.border = NA,
#     
#     panel.fun = function(x, y) {
#       
#       sector.index = get.cell.meta.data("sector.index")
#       
#       xlim = get.cell.meta.data("xlim")
#       
#       ylim = get.cell.meta.data("ylim")
#       
#     } )
#   
#   
#   
#   tmp.cat <- unique(df$lr)
#   tmp.cat.color <- c('#ffc000',"#00af50")
#   names(tmp.cat.color) <- c("ligand","receptor")
#   lapply(1:length(tmp.cat),FUN = function(ii){
#     x <- tmp.cat[ii]
#     tmp.gene.id <- df %>%
#       dplyr::filter(lr==x) %>%
#       pull(gene_id)
#     highlight.sector(tmp.gene.id, 
#                      track.index = 2,
#                      text = x, niceFacing = F,
#                      col = tmp.cat.color[x],
#                      text.col = 'black')
#   })
#   
#   
#   lrid <- data.plot %>%
#     dplyr::filter(IS.p < 0.01) %>%
#     dplyr::mutate(xx = paste0(cell_L,"_",Ligand,"_","ligand"),
#                   yy = paste0(cell_R,"_",Receptor,"_","receptor")) %>%
#     dplyr::mutate(xx=plyr::mapvalues(xx,
#                                      from = data.plot.gene$map_id,
#                                      to = data.plot.gene$gene_id),
#                   yy=plyr::mapvalues(yy,
#                                      from = data.plot.gene$map_id,
#                                      to = data.plot.gene$gene_id))
#   
#   
#   
#   for(i in 1:nrow(lrid)){
#     circos.link(sector.index1 = lrid$xx[i], 
#                 point1 = 0, 
#                 sector.index2 = lrid$yy[i], 
#                 point2 = 0 ,
#                 directional = 0,
#                 h=0.5, 
#                 lwd=0.1, 
#                 col="grey",
#                 lty=1)
#   }
#   
# }
# 
# tmp.circos.plot(tmp.res.list = tmp.res.list.blastocyst)
# 
# p <- recordPlot()
# 
# png(filename = myFileName(prefix = "res/fig/fig3_circos_blastocyst",suffix = ".png"),width = 8,height = 8,units = "in",res = 350)
# p
# title("blastocyst")
# dev.off()
# 
# 
# tmp.circos.plot(tmp.res.list = tmp.res.list.EPS_blastoid)
# p <- recordPlot()
# png(filename = myFileName(prefix = "res/fig/fig3_circos_EPS_blastoid",suffix = ".png"),width = 8,height = 8,units = "in",res = 350)
# p
# title("EPS_blastoid")
# dev.off()


####The following code is important
####------5.3.7 eLR and regulon-------------

####------5.3.7.1 load data--------

####load seurat
seu.integrated <- readRDS(file = "res/R/EPS_blastoid_20211017.rds")


####-----6.Calculate cellcall--------

####-----6.1 read rds-----------

####load seurat
seu.integrated <- readRDS(file = "res/R/EPS_blastoid_20211017.rds")
tmp.seu <- seu.integrated
levels(tmp.seu) <- c("TE","EPI","PE","intermediate")
tmp.seu$cell.type <- Idents(tmp.seu)
DefaultAssay(tmp.seu) <- "RNA"

####get count data
data.mat <- GetAssayData(tmp.seu,
                         slot = "counts",
                         assay = "RNA")
data.meta <- tmp.seu[[]] %>%
  rownames_to_column("cell.id") %>%
  group_by(cell.type) %>%
  mutate(tmp.id = row_number()) %>%
  arrange(cell.type) %>%
  ungroup() %>%
  mutate(tmp.id = paste0(tmp.id,"_",cell.type))

colnames(data.mat) <- plyr::mapvalues(x = colnames(data.mat),
                                      from = data.meta$cell.id,
                                      to = data.meta$tmp.id)

####---------6.2 perform cellcall analysis-----------
####Create CellCall object
mt <- CreateNichConObject(data = as.data.frame(data.mat), 
                          min.feature = 0,
                          names.field = 2,
                          names.delim = "_",
                          source = "UMI",
                          scale.factor = 10^6,
                          Org = "Mus musculus",
                          project = "Microenvironment")
#### Infer the cell-cell communication score
mt <- TransCommuProfile(object = mt,
                        pValueCor = 0.05,
                        CorValue = 0.1,
                        topTargetCor=1,
                        p.adjust = 0.05,
                        use.type="median",
                        probs = 0.9,
                        method="weighted",
                        IS_core = TRUE,
                        Org = 'Mus musculus')
#### Pathway activity analysis
n <- mt@data$expr_l_r_log2_scale
pathway.hyper.list <- lapply(colnames(n), function(i){
  print(i)
  tmp <- getHyperPathway(data = n, 
                         object = mt, 
                         cella_cellb = i, 
                         Org="Mus musculus")
  return(tmp)
})

myPub.df <- getForBubble(pathway.hyper.list, cella_cellb=colnames(n))
p <- plotBubble(myPub.df)+
  theme(axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18))
p
ggsave(plot = p,
       filename = myFileName(prefix = "res/fig/cellcall_Pathway_Bublleplot",
                             suffix = ".png"),
       width = 12,height = 8,dpi = 350)

####------6.3 visualization-------


####------6.3.1 Circos----------
cell_color <- data.frame(color = divergentcolor(length(levels(tmp.seu))), 
                         stringsAsFactors = FALSE)
rownames(cell_color) <- levels(tmp.seu)



png(filename = myFileName(prefix = "res/fig/cellcal_LR_circos_plot",suffix = ".png"),
    width = 8,height = 6,units = "in",res = 350)
ViewInterCircos(object = mt, 
                font = 2, 
                cellColor = cell_color, 
                lrColor = c("#F16B6F", "#84B1ED"),
                arr.type = "big.arrow",
                arr.length = 0.04,
                trackhight1 = 0.05, 
                slot="expr_l_r_log2_scale",
                linkcolor.from.sender = TRUE,
                linkcolor = NULL, 
                gap.degree = 2,
                order.vector = levels(tmp.seu),
                trackhight2 = 0.032, 
                track.margin2 = c(0.01,0.12), 
                DIY = FALSE)
dev.off()

#####--------6.3.2 pheatmap------------
### present detailed communication scores by pheatmap
dat <- mt@data[["expr_l_r_log2_scale"]]
ph <- pheatmap::pheatmap(dat, 
                   color = rdbu(100),
                   show_rownames = T, 
                   show_colnames = T, 
                   cluster_rows = T, 
                   cluster_cols = F, 
                   fontsize = 20,
                   cellwidth = 20,
                   cellheight = 15,
                   angle_col = "45", 
                   main = "score")
tmp.size <- get_pheatmap_dims(ph)
tmp.size

png(filename = myFileName(prefix = "res/fig/cellcall_pheatmap",suffix = ".png"),
    width = tmp.size$width+2.5,
    height = tmp.size$height+0.5,
    units = "in",
    res = 350)
ph
dev.off()

####-------6.3.3 sankey plot ------------
levels(tmp.seu)
mt <- LR2TF(object = mt, 
            sender_cell="intermediate", 
            recevier_cell="TE",
            slot="expr_l_r_log2_scale", 
            org="Mus musculus")
head(mt@reductions$sankey)

library(magrittr)
library(dplyr)
tmp <- mt@reductions$sankey
tmp1 <- dplyr::filter(tmp, weight1>0) 
## filter triple relation with weight1 (LR score)
tmp.df <- trans2tripleScore(tmp1)  
## transform weight1 and weight2 to one value (weight)
head(tmp.df)
## set the color of node in sankey graph
mycol.vector = c('#5d62b5','#29c3be','#f2726f','#62b58f','#bc95df', '#67cdf2', '#ffc533', '#5d62b5', '#29c3be')  
elments.num <-  tmp.df %>% unlist %>% unique %>% length()
mycol.vector.list <- rep(mycol.vector, times=ceiling(elments.num/length(mycol.vector)))
sankey_graph(df = tmp.df, axes=1:3, mycol = mycol.vector.list[1:elments.num], nudge_x = NULL, font.size = 4, boder.col="white", isGrandSon = F)
ggsave(filename = myFileName(prefix = "res/fig/cellcall_sankey_plot",
                             suffix = ".png"),
       width = 8,height = 8,dpi = 350)

####-------6.3.4 TF enrichment plot--------------

mt@data$gsea.list$intermediate@geneSets
names(mt@data$gsea.list$intermediate@geneSets)

getGSEAplot(gsea.list=mt@data$gsea.list, 
            geneSetID=c("Smad1", "Sox9", "Lef1"), 
            myCelltype="intermediate",
            fc.list=mt@data$fc.list,  
            selectedGeneID = mt@data$gsea.list$intermediate@geneSets$Smad1[1:10],
            mycol = NULL)

## gsea object
egmt <- mt@data$gsea.list$intermediate

## filter TF
egmt.df <- data.frame(egmt)
head(egmt.df[,1:6])
flag.index <- which(egmt.df$p.adjust < 0.05)
ridgeplot.DIY(x=egmt, 
              fill="p.adjust", 
              showCategory=flag.index, 
              core_enrichment = T,
              orderBy = "NES", 
              decreasing = FALSE)
ggsave(filename = myFileName(prefix = "res/fig/cellcall_show_Rigeplot",
                             suffix = ".png"),
       width = 8,height = 8,dpi = 350)

####--------7. CellChat ----------

#### using CellChat compare the two dataset

####--------7.1 load data----------
seu.integrated <- readRDS(file = "res/R/EPS_blastoid_20211017.rds")
tmp.seu <- seu.integrated
levels(tmp.seu) <- c("TE","EPI","PE","intermediate")
tmp.seu$cell.type <- Idents(tmp.seu)
DefaultAssay(tmp.seu) <- "RNA"

####-------7.2 compared the two data set----------
table(tmp.seu$orig.ident)
tmp.seu.list <- SplitObject(tmp.seu, split.by = "orig.ident")




####-------8. validate on ZGA gene---------------
####20230409

####------8.1 regenerate the 
#### load expression data
tmp.df.1 <- read.delim(file = "data/GSE71434_FPKM_inhibition.txt",sep = "\t",stringsAsFactors = F)
tmp.df.2 <- read.delim(file = "data/GSE71434_FPKM_stage.txt",sep = "\t",stringsAsFactors = F)
tmp.df <- merge(tmp.df.1,tmp.df.2,by="Gene")
zga_gene_zhangyi <- read_delim(file = "database/other_ZGA/Nfya_KD_zhangyi.txt",delim = "\t") %>%
  pull(gene)
#zga_gene_zhangyi <- intersect(zga_gene_zhangyi,mm9.gene.list)

#### load eLR data
LRpairs.df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds") %>%
  rownames_to_column("LRpairs") %>%
  mutate(Lgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ii[2]})))

tmp.id.1 <- which(LRpairs.df$Lgene %in% zga_gene_zhangyi)
LRpairs.df$Lgene[tmp.id.1]
tmp.id.2 <- which(LRpairs.df$Rgene %in% zga_gene_zhangyi)
tmp.id <- unique(c(tmp.id.1,tmp.id.2))
tmp.LR.df <- LRpairs.df[tmp.id,]

head(tmp.LR.df)
### all LR pairs ### these are expressed pairs


LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_2022031117.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 


#### Calculate ZGA
tmp.df <- tmp.df[!duplicated(tmp.df$Gene),]
data.plot <-  tmp.df %>% 
  dplyr::select(Gene,
                X2cell_early,
                X2cell_late,
                X2cell_alpha_amanitin)
rownames(data.plot) <- data.plot$Gene
data.plot <- data.plot[,-1]
colnames(data.plot) <- c("early_2cell","late_2cell","2cell_alpha_amanitin")
#data.plot <- log10(data.plot+1)

data.plot.L <- data.plot[tmp.LR.df$Lgene,]
data.plot.R <- data.plot[tmp.LR.df$Rgene,]
data.plot.IS <- sqrt(data.plot.L * data.plot.R)
rownames(data.plot.IS) <- paste0(tmp.LR.df$Lgene,"-",tmp.LR.df$Rgene)
data.plot.IS <- pheatmap:::scale_rows(data.plot.IS)
data.plot.IS[is.na(data.plot.IS)] <- 0

tmp.data.plot <- data.plot.IS


ph <- pheatmap::pheatmap(tmp.data.plot,
                         scale = "none",
                         color = rdbu(100),
                         show_rownames = T,
                         show_colnames = T,
                         cellwidth = 16,
                         cellheight = 16,
                         fontsize = 16,
                         angle_col = "45")

get_plot_dims <- function(heat_map){
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}
tmp.size <- get_plot_dims(ph)



pdf(file = myFileName(prefix = "res/fig/fig3F_ZGA_inhibit_zoomout",suffix = ".pdf"),
    width = tmp.size$width+2.5,
    height = tmp.size$height+0.5)
ComplexHeatmap::Heatmap(matrix = tmp.data.plot,
                        col = rdbu(100),
                        row_names_gp = gpar(fontsize=20,col=ifelse(rownames(tmp.data.plot) %in% c("Fgf4-Fgfr1","Fgf4-Fgfr2"),"red","black"),
                                            fontface="italic"),
                        rect_gp = gpar(col = "black"),
                        column_names_gp = gpar(fontsize=20),
                        column_names_rot = 45,
                        column_names_centered = F,
                        heatmap_legend_param = list(title="zscore",
                                                    title_gp=gpar(fontsize=20,fontface="bold"),
                                                    labels_gp=gpar(fontsize=16,fontface="bold")))
dev.off()


png(file = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_zoomout",suffix = ".png"),
    width = tmp.size$width+2,
    height = tmp.size$height+2,
    units = "in",res = 350)
ComplexHeatmap::Heatmap(matrix = tmp.data.plot,
                        col = rdbu(100),
                        row_names_gp = gpar(fontsize=20,
                                            col=ifelse(rownames(tmp.data.plot) %in% c("Fgf4-Fgfr1","Fgf4-Fgfr2"),"red","black"),
                                            fontface="italic"),
                        rect_gp = gpar(col = "black"),
                        column_names_gp = gpar(fontsize=20),
                        column_names_rot = 45,
                        column_names_centered = F,
                        heatmap_legend_param = list(title="zscore",
                                                    title_gp=gpar(fontsize=20,fontface="bold"),
                                                    labels_gp=gpar(fontsize=20)))
dev.off()


png(filename = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_pheatmap_filter",suffix = ".png"),
    width = 6,height = 8,units = "in",res = 350)
tmp.idx <- annotation_row %>%
  dplyr::filter(Cluster %in% c(4,6)) %>%
  rownames()
num_clusters <- 6
ph <- pheatmap::pheatmap(data.plot.IS[tmp.idx,],
               scale = "none",
               color = rdbu(100),
               show_rownames = T,
               show_colnames = T,
               cutree_rows = num_clusters,
               annotation_row = annotation_row,
               silent = T,
               angle_col = 45)
annotation_row.1 <- data.frame(Cluster = factor(cutree(ph$tree_row,num_clusters)))
ph <- pheatmap::pheatmap(data.plot.IS[tmp.idx,],
                         scale = "none",
                         color = rdbu(100),
                         show_rownames = F,
                         show_colnames = T,
                         cutree_rows = num_clusters,
                         annotation_row = annotation_row.1,
                         silent = T,
                         fontsize = 28,
                         angle_col = 45)
ph
dev.off()


tmp.idx <- annotation_row.1 %>%
  dplyr::filter(Cluster %in% c(2,3,6)) %>%
  rownames()

png(filename = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_zoomout",suffix = ".png"),
    width = 3,height = 9,units = "in",res = 350)
pheatmap::pheatmap(data.plot.IS[tmp.idx,],
         scale = "none",
         color = rdbu(100),
         show_rownames = T,
         show_colnames = T,
         silent = F,
         angle_col = 45)
dev.off()

pdf(file = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_zoomout",suffix = ".pdf"),
    width = 6,height = 12)
ComplexHeatmap::pheatmap(data.plot.IS[tmp.idx,],
                   scale = "none",
                   color = rdbu(100),
                   show_rownames = T,
                   show_colnames = T,
                   fontsize = 16,
                   name = "IS",
                   angle_col = "45")
dev.off()



png(file = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_zoomout",suffix = ".png"),
    width = tmp.size$width+0.5,
    height = tmp.size$height+0.5,
    units = "in",res = 350)
ComplexHeatmap::pheatmap(data.plot.IS[tmp.idx,],
                         scale = "none",
                         color = rdbu(100),
                         show_rownames = T,
                         show_colnames = T,
                         fontsize = 16,
                         heatmap_legend_param = list(title_gp=gpar(fontsize=16,fontface="bold"),
                                                     labels_gp=gpar(fontsize=16,fontface="bold")),
                         name = "IS",
                         angle_col = "45")
dev.off()

pdf(file = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_zoomout",suffix = ".pdf"),
    width = tmp.size$width+0.5,
    height = tmp.size$height+0.5)
ComplexHeatmap::pheatmap(data.plot.IS[tmp.idx,],
                         scale = "none",
                         color = rdbu(100),
                         show_rownames = T,
                         show_colnames = T,
                         fontsize = 16,
                         heatmap_legend_param = list(title_gp=gpar(fontsize=16,fontface="bold"),
                                                     labels_gp=gpar(fontsize=16,fontface="bold")),
                         name = "IS",
                         angle_col = "45")
dev.off()



tmp.data.plot <- data.plot.IS[tmp.idx,]
rownames(tmp.data.plot) <- gsub(rownames(tmp.data.plot),pattern = "_",replacement = "-")
colnames(tmp.data.plot) <- c("early 2-cell","late 2-cell","2-cell alpha amanitin")

tmp.data.plot[tmp.data.plot > 1] <- 1
pdf(file = myFileName(prefix = "res/fig/fig3F_ZGA_inhibit_zoomout",suffix = ".pdf"),
    width = tmp.size$width+2.5,
    height = tmp.size$height+0.5)
ComplexHeatmap::Heatmap(matrix = tmp.data.plot,
                        col = rdbu(100),
                        row_names_gp = gpar(fontsize=20,col=ifelse(rownames(tmp.data.plot) %in% c("Fgf4-Fgfr1","Fgf4-Fgfr2"),"red","black"),
                                            fontface="italic"),
                        rect_gp = gpar(col = "black"),
                        column_names_gp = gpar(fontsize=20),
                        column_names_rot = 45,
                        column_names_centered = F,
                        heatmap_legend_param = list(title="IS",
                                                    title_gp=gpar(fontsize=20,fontface="bold"),
                                                    labels_gp=gpar(fontsize=16,fontface="bold")))
dev.off()


png(file = myFileName(prefix = "res/fig/fig3_ZGA_inhibit_zoomout",suffix = ".png"),
    width = tmp.size$width+2,
    height = tmp.size$height+2,
    units = "in",res = 350)
ComplexHeatmap::Heatmap(matrix = tmp.data.plot,
                        col = rdbu(100),
                        row_names_gp = gpar(fontsize=20,
                                            col=ifelse(rownames(tmp.data.plot) %in% c("Fgf4-Fgfr1","Fgf4-Fgfr2"),"red","black"),
                                            fontface="italic"),
                        rect_gp = gpar(col = "black"),
                        column_names_gp = gpar(fontsize=20),
                        column_names_rot = 45,
                        column_names_centered = F,
                        heatmap_legend_param = list(title="IS",
                                                    title_gp=gpar(fontsize=20,fontface="bold"),
                                                    labels_gp=gpar(fontsize=20)))
dev.off()

# ####---------export data--------
# tmp.df <- readRDS(file = "res/R/putative_eLR_correlation_20210129.Rds")
# LRpairs <- tmp.df %>%
#   dplyr::filter(abs(PCC) > 0.1)
# write.table(LRpairs,file = "res/txt/Additional_file_2_eLR_PCC_value.txt",quote = F,sep = "\t",row.names = F)



###-------9.Clearly checked the result of ZGA by boxplot--------

####20230409

#### load expression data
tmp.df.1 <- read.delim(file = "data/GSE71434_FPKM_inhibition.txt",sep = "\t",stringsAsFactors = F)
tmp.df.2 <- read.delim(file = "data/GSE71434_FPKM_stage.txt",sep = "\t",stringsAsFactors = F)
tmp.df <- merge(tmp.df.1,tmp.df.2,by="Gene")
zga_gene_zhangyi <- read_delim(file = "database/other_ZGA/Nfya_KD_zhangyi.txt",delim = "\t") %>%
  pull(gene)
#zga_gene_zhangyi <- intersect(zga_gene_zhangyi,mm9.gene.list)

#### load eLR data
eLRpairs.df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds") %>%
  rownames_to_column("LRpairs") %>%
  mutate(Lgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ii[2]})))

# tmp.id.1 <- which(LRpairs.df$Lgene %in% zga_gene_zhangyi)
# LRpairs.df$Lgene[tmp.id.1]
# tmp.id.2 <- which(LRpairs.df$Rgene %in% zga_gene_zhangyi)
# tmp.id <- unique(c(tmp.id.1,tmp.id.2))
# tmp.LR.df <- LRpairs.df[tmp.id,]
tmp.LR.df <- eLRpairs.df

head(tmp.LR.df)

### all LR pairs ### these are expressed pairs
 

#### Calculate ZGA
tmp.df <- tmp.df[!duplicated(tmp.df$Gene),]
data.plot <-  tmp.df %>% 
  dplyr::select(Gene,
                X2cell_early,
                X2cell_late,
                X2cell_alpha_amanitin)
rownames(data.plot) <- data.plot$Gene
data.plot <- data.plot[,-1]
colnames(data.plot) <- c("early_2cell","late_2cell","2cell_alpha_amanitin")
#data.plot <- log10(data.plot+1)

tmp.LR.df <- eLRpairs.df
data.plot.L <- data.plot[tmp.LR.df$Lgene,]
data.plot.R <- data.plot[tmp.LR.df$Rgene,]
data.plot.IS <- sqrt(data.plot.L * data.plot.R)
rownames(data.plot.IS) <- paste0(tmp.LR.df$Lgene,"-",tmp.LR.df$Rgene)
# data.plot.IS <- pheatmap:::scale_rows(data.plot.IS)
# data.plot.IS[is.na(data.plot.IS)] <- 0

data.plot.IS <- log10(data.plot.IS+1)
head(data.plot.IS)
boxplot(data.plot.IS,main = "eLR (430 LR pairs)",ylab = "log10(IS+1)")

data.plot.IS <- data.plot.IS %>%
  dplyr::mutate(delta_res = `2cell_alpha_amanitin` - late_2cell)
plot(density(data.plot.IS$delta_res))

boxplot(data.plot.IS$delta_res,
        main = "eLR (430 LR pairs)",
        ylab = "log10(IS+1)")

#min(data.plot.IS$cell_alpha_amanitin)
### not eLR all 

data.plot.IS.1 <- data.plot.IS



LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_2022031117.rds")
LRpairs <- gsub(pattern = "_",replacement = "-",LRpairs)

not_eLR <- setdiff(LRpairs,eLRpairs.df$LRpairs)

Lgenelist <- unlist(lapply(strsplit(not_eLR,split = "-"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(not_eLR,split = "-"),FUN = function(x) x[2])) 

data.plot.L <- data.plot[Lgenelist,]
data.plot.R <- data.plot[Rgenelist,]
data.plot.IS <- sqrt(data.plot.L * data.plot.R)
rownames(data.plot.IS) <- paste0(Lgenelist,"-",Rgenelist)
# data.plot.IS <- pheatmap:::scale_rows(data.plot.IS)
# data.plot.IS[is.na(data.plot.IS)] <- 0

#boxplot(data.plot.IS,main="1409 not eLR",ylab = "IS")
boxplot(log10(data.plot.IS+1),
        main = "not eLR (1409 LR pairs)",
        ylab = "log10(IS+1)")

data.plot.IS <- data.plot.IS %>%
  dplyr::mutate(delta_res = `2cell_alpha_amanitin` - late_2cell)
plot(density(data.plot.IS$delta_res))

boxplot(data.plot.IS$delta_res)
###

data.plot.IS.2 <- data.plot.IS



boxplot(data.plot.IS.1$delta_res,
        data.plot.IS.2$delta_res,ylim = c(-1,1))


boxplot(data.plot.IS.1$delta_res)
boxplot(data.plot.IS.2$delta_res)




### test  which group of delta value is higher,
### This is really important, it can check which one is better;
### The left is really important for the right ones
sd(data.plot.IS.1$delta_res)
mean(data.plot.IS.1$delta_res)
#### -0.1786005
sd(data.plot.IS.2$delta_res)
mean(data.plot.IS.2$delta_res)
#### -0.3211708



tmp.p.res <- t.test(abs(data.plot.IS.1$delta_res),
                    abs(data.plot.IS.2$delta_res),
                    alternative = "g")
tmp.p.res <- wilcox.test(abs(data.plot.IS.1$delta_res),
                    abs(data.plot.IS.2$delta_res),
                    alternative = "g")

tmp.p.res

tmp.p.res <- signif(tmp.p.res$p.value,digits = 5)

df_p_val <- data.frame(group1 = "eLR",
                       group2 = "not_eLR",
                       label = tmp.p.res,
                       y.position = "y.position")

tmp.data.plot <- data.frame(delta_value = c(abs(data.plot.IS.1$delta_res),
                                            abs(data.plot.IS.2$delta_res)),
                            group = c(rep("eLR",length(data.plot.IS.1$delta_res)),
                                      rep("not_eLR",length(data.plot.IS.2$delta_res))),
                            stringsAsFactors = F)


range(tmp.data.plot$delta_value)


tmp.bar.df <- data.frame(delta_value = c())


# This version of plot will be perturb by range y, so, it not robust
# p <- ggplot(tmp.data.plot,aes(group,delta_value))+
#   geom_violin(aes(group,delta_value,fill = group))+
#   # geom_quasirandom(alpha = 0.1)+
#   # stat_summary(fun.y=mean,
#   #              geom = "point",
#   #              shape = 20,
#   #              size = 12,
#   #              color = "black",
#   #              fill = "black") +
#   stat_summary(tmp.data.plot,
#                fun.y = mean,
#                geom = "crossbar",
#                color = "red",
#                width = 0.2) +
#   geom_signif(comparisons = list(c("eLR","not_eLR")),
#               annotations = paste0("p=",tmp.p.res),
#               y_position = 3,
#               size = 1,
#               tip_length = 0.001,
#               textsize = 10)+
#   xlab(NULL)+
#   ylab("delta value")+
#   scale_fill_manual(values = c("#8FD3F7","#E8E0D1"))+
#   scale_y_continuous(limits = c(-1,1))+
#   theme_cowplot(font_size = 22)+
#   theme(legend.position = "none")
# p

# This version of function is mannual calculate the linerange
head(tmp.data.plot)

# A very important idea, calculate the summary of characteristic of value
df.sum <- tmp.data.plot %>%
  group_by(group) %>%
  summarise(mean = mean(delta_value),
            sd = sd(delta_value)) %>%
  ungroup()

df.sum
df.sum %>%
  group_by(group) %>%
  mutate(ymin = (mean-sd),
         ymax = (mean+sd),
         range = ymax - ymin) %>%
  ungroup() %>%
  dplyr::filter(range == max(range))


p <- ggplot()+
  geom_violin(data = tmp.data.plot,
              mapping = aes(group,delta_value,
                            fill = group))+
  geom_pointrange(data = df.sum,
                  color = "black",
                  mapping = aes(x = group,
                                y = mean,
                               ymin = (mean -sd),
                               ymax = (mean + sd)))+
  geom_text_repel(data = df.sum,
                  point.padding = 0.2,
                  nudge_x = .5,
                  nudge_y = 0,
                  size = 6,
                  segment.curvature = 0,
                  arrow = arrow(length = unit(0.015, "npc")),
                  mapping = aes(x = group,
                                y = mean,
                                label = signif(mean,3)))+
  add_pvalue(data.frame(group1 = "eLR",
                        group2 =  "not_eLR",
                        label = tmp.p.res,
                        y.position = 2.3),
             label = "p={label}",
             tip.length = 0,
             label.size = 6 )+
  xlab(NULL)+
  ylab("delta value")+
  scale_fill_manual(values = c("#8FD3F7","#E8E0D1"))+
  scale_y_continuous(limits = c(-2.1,2.3))+
  theme_cowplot(font_size = 22)+
  theme(legend.position = "none",
        axis.ticks = element_line(size = 1),
        axis.line = element_line(size = 1))

p

ggsave(filename = myFileName(prefix = "res/fig/ZGA_inhibit_boxplot",
                             suffix = ".jpg"),
       width = 6,height = 6,dpi = 350)
ggsave(filename = myFileName(prefix = "res/fig/ZGA_inhibit_boxplot",
                             suffix = ".pdf"),
       width = 6,height = 6)


