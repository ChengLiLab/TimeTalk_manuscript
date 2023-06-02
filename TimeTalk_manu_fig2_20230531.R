####Fig2
####20220311
####author:wlt
####TimeTalk method

for(ii in 1:100000){
  cat(ii,sep = "\n")
  Sys.sleep(time = ii/100)
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
#library(monocle)
library(netSmooth)
library(clusterProfiler)
library(org.Mm.eg.db)
library(entropy)
library(RTN)
library(RTNduals)
library(snow)
library(ComplexHeatmap)
library(msigdbr)
library(patchwork)
library(TimeTalk)
source(file = "code/myUtils.R")


####--------1. Ligand, Receptor expression-----------
#####-------1.1 show Ligand, Receptor gene expression landscape--------
#####--------1.1.1 prepare for Ligand, Receptor database-----------

LRpairs.df <- read.delim(file = "database/Ligand-Receptor-Pairs/Mouse/Mouse-2020-Shao-LR-pairs.txt",
                         stringsAsFactors = F)
LRpairs <- LRpairs.df$lr_pair
Lgenelist <- LRpairs.df$ligand_gene_symbol
Rgenelist <- LRpairs.df$receptor_gene_symbol 

data.merge.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
data.merge.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")

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

#####------1.1.2 prepare data for plot-------
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_2022031117.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)

data.merge.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
data.merge.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
colnames(data.merge.meta)[1] <- "Cell"

#### don't forget normalized data
data.merge.exp  <- myRemoveNA(data.merge.exp)
data.merge.exp  <- myNormalize(data.merge.exp)

tmp.row.names <- rownames(data.merge.exp)
tmp.col.names <- colnames(data.merge.exp)

data.merge.exp <- preprocessCore::normalize.quantiles(as.matrix(data.merge.exp))
rownames(data.merge.exp) <- tmp.row.names
colnames(data.merge.exp) <- tmp.col.names

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
  dplyr::filter(is.na(gene.exp)) %>%
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
                                     "fibroblast"))

saveRDS(data.plot,file = myFileName(prefix = "res/R/data_merge_plot",suffix = ".rds"))


####--------------1.2 gene mean expression changes------------------

####-------------1.2.1 calc mean expression-------------
####mean expression
####load data
data.plot <- readRDS(file = "res/R/data_merge_plot_2022031117.rds")

#####Background
data.plot.mean.Background <- data.plot %>%
  dplyr::filter(group=="Background") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp),
            gene.exp.sd=sd(gene.exp),
            gene.exp.cv=ifelse(mean(gene.exp)!=0,sd(gene.exp)/mean(gene.exp),0)) %>%
  ungroup()
#####Ligand
data.plot.mean.Ligand <- data.plot %>%
  dplyr::filter(group=="Lgene") %>%
  group_by(Gene,Stage) %>%
  summarise(gene.exp.mean=mean(gene.exp),
            gene.exp.sd=sd(gene.exp),
            gene.exp.cv=ifelse(mean(gene.exp)!=0,sd(gene.exp)/mean(gene.exp),0)) %>%
  ungroup()

#####Receptor
data.plot.mean.Receptor <- data.plot %>%
  dplyr::filter(group=="Rgene") %>%
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

data.plot.mean <- data.plot.mean %>%
  dplyr::filter(Stage != "fibroblast")

#####------1.2.2 plot-------


#####-----1.2.2.1 mean--------
tmp.data.plot <- data.plot.mean
levels(tmp.data.plot$group)
mypalette <- c("#ffc000","#00b050","lightgrey")



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
               fill = mypalette[1], 
               color=mypalette[1],
               outlier.shape = NA,
               width = 0.6)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.mean=gene.exp.mean), 
               geom = "point", 
               fun = mean)+
  stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
                 mutate(gene.exp.mean=gene.exp.mean), 
               geom = "line", 
               color="black",
               aes(group=""),
               fun = mean)+
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
  geom_hline(yintercept = 0,linetype="solid",size=1)+
  xlab(NULL)+
  ylab("gene.exp.mean")+
  scale_y_continuous(limits = c(-6,6),
                     breaks = seq(-6,6,2),
                     labels = as.character(c(-seq(-6,-2,by = 2),0,seq(2,6,by = 2))))+
  scale_x_discrete(labels=tmp.stage <- c("MII oocyte","zygote",
                                         "early 2-cell","mid 2-cell","late 2-cell",
                                         "4-cell","8-cell","16-cell",
                                         "early blastocyst","mid blastocyst","late blastocyst"))+
  theme_cowplot(font_size  = 28)+
  theme(legend.position = "right",
        axis.ticks.y = element_line(size = 2),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line( size = 2 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))


p2
myggsave(p = p2,prefix = "res/fig/Fig2_landscape_gene_mean_exp",
         suffix = ".png",
         width = 8,height = 6,dpi=350)
myggsave(p = p2,prefix = "res/fig/Fig2_landscape_gene_mean_exp",
         suffix = ".pdf",
         width = 8,height = 6)



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
  geom_hline(yintercept = 0,linetype="solid",size=1)+
  scale_x_discrete(labels=tmp.stage <- c("MII oocyte","zygote",
                                         "early 2-cell","mid 2-cell","late 2-cell",
                                         "4-cell","8-cell","16-cell",
                                         "early blastocyst","mid blastocyst","late blastocyst"))+
  scale_y_continuous(limits = c(-5,5),
                     breaks = seq(-5,5,1),
                     labels = as.character(c(-seq(-5,-1,by = 1),0,seq(1,5,by = 1))))+
  ylab("gene.exp.sd")+
  xlab(NULL)+
  theme_cowplot(font_size  = 28)+
  theme(legend.position = "right",
        axis.ticks.y = element_line(size = 2),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line( size = 2 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))

myggsave(p = p3,prefix = "res/fig/FigS3_landscape_gene_sd",suffix = ".png",width = 8,height = 6,dpi=350)
myggsave(p = p3,prefix = "res/fig/FigS3_landscape_gene_sd",suffix = ".pdf",width = 8,height = 6,dpi=350)



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
  geom_hline(yintercept = 0,linetype="solid",size=1)+
  scale_x_discrete(labels=tmp.stage <- c("MII oocyte","zygote",
                                         "early 2-cell","mid 2-cell","late 2-cell",
                                         "4-cell","8-cell","16-cell",
                                         "early blastocyst","mid blastocyst","late blastocyst"))+
  scale_y_continuous(limits = c(-30,30),
                     breaks = seq(-30,30,10),
                     labels = as.character(c(-seq(-30,-10,by = 10),0,seq(10,30,by = 10))))+
  ylab("gene.exp.cv")+
  xlab(NULL)+
  theme_cowplot(font_size  = 28)+
  theme(legend.position = "right",
        axis.ticks.y = element_line(size = 2),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_line( size = 2 ),
        axis.text.x = element_text(hjust = 0.8,
                                   angle = 30))
#p4 
myggsave(p = p4,prefix = "res/fig/FigS3_gene_cv",suffix = ".png",width = 8,height = 6,dpi=300)
myggsave(p = p4,prefix = "res/fig/FigS3_gene_cv",suffix = ".pdf",width = 8,height = 6,dpi=300)


# #####---------1.2.2.2 rank---------------
# 
# LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_20210115.rds")
# Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
# Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
# 
# Lgene <- unique(Lgenelist)
# Rgene <- unique(Rgenelist)
# 
# data.plot.rank.Background <- data.plot.mean.Background %>%
#   group_by(Stage) %>%
#   mutate(rank=row_number(gene.exp.mean)) %>%
#   ungroup() %>%
#   dplyr::select(-group,-gene.exp.mean,-gene.exp.sd,-gene.exp.cv) %>%
#   spread(key = "Stage",value = "rank") %>%
#   column_to_rownames("Gene")
# 
# data.plot.rank.L <- data.plot.rank.Background[Lgene,] %>%
#   rownames_to_column("Gene") %>%
#   gather(key = "Stage",value = "rank",-Gene) %>%
#   mutate(group="Lgene")
# 
# data.plot.rank.R <- data.plot.rank.Background[Rgene,] %>%
#   rownames_to_column("Gene") %>%
#   gather(key = "Stage",value = "rank",-Gene) %>%
#   mutate(group="Rgene")
# 
# data.plot.rank <- data.plot.rank.Background %>%
#   rownames_to_column("Gene") %>%
#   gather(key = "Stage",value = "rank",-Gene) %>%
#   mutate(group="Background")
# 
# data.plot.rank <- rbind(data.plot.rank,data.plot.rank.L,data.plot.rank.R)
# 
# data.plot.rank$Stage <- factor(data.plot.rank$Stage,
#                                levels = c("Oocyte","zygote","early2cell","mid2cell","late2cell",
#                                           "4cell","8cell","16cell","earlyblast","midblast","lateblast",
#                                           "E5.25","E5.5","E6.25","E6.5","fibroblast"))
# 
# data.plot.rank$group <- factor(data.plot.rank$group,
#                                levels = c("Lgene","Rgene","Background"))
# 
# 
# 
# tmp.data.plot <- data.plot.rank
# head(tmp.data.plot)
# mypalette <- c("#ffc000","#00b050","lightgrey")
# names(mypalette) <- c("Ligand","Receptor","Background")
# 
# 
# tmp.data.plot.bar <- subset(tmp.data.plot,group == "Background") %>%
#   mutate(rank=rank) %>%
#   group_by(Stage) %>%
#   summarise(rank=n()) %>%
#   ungroup()
# 
# levels(tmp.data.plot$Stage)[1:11]
# 
# tmp.stage <- c("MII oocyte","zygote",
#   "early 2-cell","mid 2-cell","late 2-cell",
#   "4-cell","8-cell","16-cell",
#   "early blastocyst","mid blastocyst","late blastocyst",
#   "E5.25","E5.5","E6.25","E6.5","fibroblast")
# p5 <- ggplot(data = tmp.data.plot,aes(Stage,rank))+
#   # geom_bar(data = tmp.data.plot.bar,
#   #          mapping = aes(x = Stage,y = rank),fill = mypalette[3], color=mypalette[3],
#   #          stat = "identity",width = 0.8)+
#   stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
#                  mutate(rank=rank), geom = "point",shape = 15, fun = mean)+
#   stat_summary(data = subset(tmp.data.plot,group == "Lgene") %>%
#                  mutate(rank=rank), geom = "line", color=mypalette[1],aes(group=""),fun = mean)+
#   stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
#                  mutate(rank=rank), geom = "point",shape = 15, fun = mean)+
#   stat_summary(data = subset(tmp.data.plot,group == "Rgene") %>%
#                  mutate(rank=rank), geom = "line", color=mypalette[2],aes(group=""),fun = mean)+
#   geom_vline(xintercept = "lateblast",linetype="dashed",size=1)+
#   scale_y_continuous(breaks = seq(6500,8500,500))+
#   scale_x_discrete(labels=tmp.stage)+
#   ylab("rank")+
#   xlab(NULL)+
#   theme_cowplot(font_size  = 28)+
#   theme(legend.position = "right",
#         axis.ticks = element_line(size = 2),
#         axis.line.x = element_line( size = 2 ),
#         axis.line.y = element_line( size = 2 ),
#         axis.text.x = element_text(hjust = 0.8,
#                                    angle = 35))
# p5
# myggsave(p = p5,prefix = "res/fig/Fig2_gene_exp.mean_rank",suffix = ".jpg",width = 8,height = 6,dpi=350)
# 
# 
# 
# x <- tmp.data.plot %>%
#   dplyr::filter(group=="Lgene") %>%
#   group_by(Stage) %>%
#   summarise(ttt=mean(rank)) %>%
#   pull(ttt)
# y <- tmp.data.plot %>%
#   dplyr::filter(group=="Rgene") %>%
#   group_by(Stage) %>%
#   summarise(ttt=mean(rank)) %>%
#   pull(ttt)
# tmp.res <- cor.test(x[1:11],y[1:11])
# tmp.res$p.value
# tmp.stage
# 
# 
# p <- p5+annotate("text",
#                  x = "4cell",
#                  y = 8500,
#                  label=paste0("PCC=",round(cor(x[1:11],y[1:11]),4),"\n",
#                               "p=",round(cor.test(x[1:11],y[1:11])$p.value,4)),
#                  size=10)
# p
# myggsave(p = p,prefix = "res/fig/Fig2_gene_exp.mean_rank",suffix = ".jpg",
#          width = 8,height = 6,dpi=350)
# myggsave(p = p,prefix = "res/fig/Fig2B_gene_exp.mean_rank",suffix = ".pdf",
#          width = 8,height = 6,dpi=350)
# 





# beforeEPI.exp.mean <- beforeEPI.exp %>%
#   rownames_to_column("Gene") %>%
#   gather(key = "cell_id",value = "gene.exp",-Gene) %>%
#   left_join(beforeEPI.meta,"cell_id") %>%
#   group_by(Gene,Stage) %>%
#   summarise(gene.exp.mean=mean(gene.exp)) %>%
#   ungroup() %>%
#   mutate(Stage=factor(Stage,levels = Stage.level)) %>%
#   spread(key="Stage",value = "gene.exp.mean") %>%
#   column_to_rownames("Gene")
# 
# boxplot(beforeEPI.exp.mean)


####-------2. visualize gene expression---------

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
mat.df <- beforeEPI.exp[,tmp.cell.id] %>%
  rownames_to_column("gene_symbol") %>%
  gather(key = "cell_id",value = "gene.exp",-gene_symbol)

mypalette <- c("#ffc000","#00b050","lightgrey")



tmp.LR.plot.function <- function(tmp.Lgene,tmp.Rgene,mat.df,tmp.meta){
  
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
  
  tmp.data.plot <- rbind(tmp.data.plot.1,tmp.data.plot.2)
  
  ###plot(tmp.data.plot.1$gene.exp,tmp.data.plot.2$gene.exp,col = warm(262))
  
  
  
  tmp.pcc <- cor(tmp.data.plot.1$gene.exp,
                 tmp.data.plot.2$gene.exp,method = "pearson")
  tmp.scc <- cor(tmp.data.plot.1$gene.exp,
                 tmp.data.plot.2$gene.exp,method = "spearman")
  
  p <- ggplot(tmp.data.plot,aes(tmp.rank,gene.exp,group=group,col=group))+
    geom_point(alpha=0.3)+
    geom_smooth(se = F,method = "gam", formula = y ~ s(x, bs = "cs"),alpha=0.3)+
    xlab("rank by pseudotime")+
    ylab("log2(RPKM+1)")+
    ggtitle(label = paste0(tmp.Lgene,"-",tmp.Rgene,"\n",
                           "PCC=",signif(tmp.pcc,digits = 4),",",
                           "SCC=",signif(tmp.scc,digits = 4)))+
    scale_color_manual(values = mypalette[1:2])+
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

# tmp.LR.plot.function <- function(tmp.Lgene,tmp.Rgene,mat.df,tmp.meta){
#   
#   tmp.data.plot.1 <- mat.df %>%
#     dplyr::filter(gene_symbol == tmp.Lgene) %>%
#     left_join(beforeEPI.meta,by = "cell_id") %>%
#     dplyr::select(cell_id,gene_symbol,gene.exp,Pseudotime,cell_type,Stage) %>%
#     mutate(tmp.rank = row_number(Pseudotime)) %>%
#     arrange(Stage,tmp.rank) %>%
#     mutate(tmp.rank = row_number()) %>%
#     mutate(group = "Ligand")
#   
#   tmp.data.plot.2 <- mat.df %>%
#     dplyr::filter(gene_symbol == tmp.Rgene) %>%
#     left_join(beforeEPI.meta,by = "cell_id") %>%
#     dplyr::select(cell_id,gene_symbol,gene.exp,Pseudotime,cell_type,Stage) %>%
#     mutate(tmp.rank = row_number(Pseudotime)) %>%
#     arrange(Stage,tmp.rank) %>%
#     mutate(tmp.rank = row_number()) %>%
#     mutate(group = "Receptor")
#   
#   tmp.data.plot <- rbind(tmp.data.plot.1,tmp.data.plot.2)
#   
#   tmp.pcc <- cor(tmp.data.plot.1$gene.exp,
#                  tmp.data.plot.2$gene.exp,method = "pearson")
#   tmp.scc <- cor(tmp.data.plot.1$gene.exp,
#                  tmp.data.plot.2$gene.exp,method = "spearman")
#   
#   p <- ggplot(tmp.data.plot,aes(tmp.rank,gene.exp,group=group,col=group))+
#     geom_point(alpha=0.3)+
#     geom_smooth(se = F,method = "gam", formula = y ~ s(x, bs = "cs"),alpha=0.3)+
#     xlab("rank by pseduotime")+
#     ylab("log2(RPKM+1)")+
#     ggtitle(label = paste0(tmp.Lgene,"-",tmp.Rgene,"\n",
#                            "PCC=",signif(tmp.pcc,digits = 4),",",
#                            "SCC=",signif(tmp.scc,digits = 4)))+
#     scale_color_manual(values = mypalette[1:2])+
#     scale_x_continuous(expand = c(0,0))+
#     scale_y_continuous(expand = c(0,0))+
#     theme_cowplot(font_size = 28)+
#     theme(legend.position = "top",legend.justification = "center",
#           axis.line = element_line(size = 2),
#           axis.ticks = element_line(size = 2),
#           plot.title = element_text(hjust = 0.5))
#   
#   return(p)
# }




p1 <- tmp.LR.plot.function(tmp.Lgene = "Bmp6",
                           tmp.Rgene = "Bmpr1b",
                           mat.df = mat.df,
                           tmp.meta = beforeEPI.meta)

p1
myggsave(p = p1,prefix = "res/fig/fig2_Bmp6_Bmpr1b_show_eLR",suffix = ".png",width = 8,height = 8,dpi=350)
myggsave(p = p1,prefix = "res/fig/fig2_Bmp6_Bmpr1b_show_eLR",suffix = ".pdf",width = 8,height = 8,dpi=350)

p1 <- tmp.LR.plot.function(tmp.Lgene = "Adam17",
                           tmp.Rgene = "Cd9",
                           mat.df = mat.df,
                           tmp.meta = beforeEPI.meta)
myggsave(p = p1,prefix = "res/fig/fig2_Adam17_Cd9_show_eLR",suffix = ".png",width = 8,height = 8,dpi=350)
myggsave(p = p1,prefix = "res/fig/fig2_Adam17_Cd9_show_eLR",suffix = ".pdf",width = 8,height = 8,dpi=350)

p1 <- tmp.LR.plot.function(tmp.Lgene = "Fgf4",
                           tmp.Rgene = "Fgfr1",
                           mat.df = mat.df,
                           tmp.meta = beforeEPI.meta)
myggsave(p = p1,prefix = "res/fig/fig2_Fgf4_Fgfr1_show_eLR",suffix = ".png",width = 8,height = 8,dpi=350)
myggsave(p = p1,prefix = "res/fig/fig2_Fgf4_Fgfr1_show_eLR",suffix = ".pdf",width = 8,height = 8,dpi=350)


p1 <- tmp.LR.plot.function(tmp.Lgene = "Fgf4",
                           tmp.Rgene = "Fgfr2",
                           mat.df = mat.df,
                           tmp.meta = beforeEPI.meta)
myggsave(p = p1,prefix = "res/fig/fig2_Fgf4_Fgfr2_show_eLR",suffix = ".png",width = 8,height = 8,dpi=350)
myggsave(p = p1,prefix = "res/fig/fig2_Fgf4_Fgfr2_show_eLR",suffix = ".pdf",width = 8,height = 8,dpi=350)

p1 <- tmp.LR.plot.function(tmp.Lgene = "Igf2",
                           tmp.Rgene = "Igf2r",
                           mat.df = mat.df,
                           tmp.meta = beforeEPI.meta)
myggsave(p = p1,prefix = "res/fig/fig2_Igf2_Igf2r_show_eLR",suffix = ".png",width = 8,height = 8,dpi=350)
myggsave(p = p1,prefix = "res/fig/fig2_Igf2_Igf2r_show_eLR",suffix = ".pdf",width = 8,height = 8,dpi=350)


p1 <- tmp.LR.plot.function(tmp.Lgene = "Pdgfc",
                           tmp.Rgene = "Pdgfra",
                           mat.df = mat.df,
                           tmp.meta = beforeEPI.meta)
myggsave(p = p1,prefix = "res/fig/fig2_Pdgfc_Pdgfra_show_eLR",suffix = ".png",width = 8,height = 8,dpi=350)
myggsave(p = p1,prefix = "res/fig/fig2_Pdgfc_Pdgfra_show_eLR",suffix = ".pdf",width = 8,height = 8,dpi=350)






#####---------3. identify eLR directly using TimeTalk-------------

##### insert in 20220620-20220629

#####---------3.1 load data ---------------
tmp.seu <- readRDS(file = "res/R/early.scRNAseq.seurat_2022031015.rds")
tmp.seu <- RunALRA(tmp.seu)

tmp.cds <- readRDS(file = "res/R/early_embryo_monocle3_object_2022062221.rds")
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_2022031117.rds")
Idents(tmp.seu) <- "Stage"
levels(tmp.seu)

#####--------3.2 get variable gene-------------

tmp.GetVariableGeneInTrajectory <- function(tmp.cds, tmp.cores = 10, tmp.neighbor_graph = "knn", 
          tmp.q_value_cutoff = 0.05, tmp_moran_I_cutoff = 0.2) {
  
  pr_graph_test_res <- monocle3::graph_test(cds = tmp.cds, neighbor_graph = tmp.neighbor_graph, cores = tmp.cores)
  
  tmp.hits <- pr_graph_test_res %>% 
    dplyr::filter(q_value < tmp.q_value_cutoff & morans_I > tmp_moran_I_cutoff & status == "OK") %>% 
    pull(gene_short_name) %>% 
    as.vector() %>% 
    as.character()
  
  tmp.phenotype <- pr_graph_test_res %>% 
    dplyr::filter(q_value < tmp.q_value_cutoff & morans_I > tmp_moran_I_cutoff & status == "OK") %>% 
    pull(morans_test_statistic)
  
  names(tmp.phenotype) <- tmp.hits
  return(tmp.phenotype)
}

tmp.GetCellTypeSpecificTRN <- function (tmp.seu, tmp.ident, gene.TF, tmp.phenotype,use.impute = F) {
  tmp.cell.type <- tmp.ident
  cat(paste0(tmp.cell.type, " TRN reconsturct start:"), sep = "\n")
  if( tmp.ident == "all" ){
    tmp.ident <- levels(tmp.seu)
  }
  
  if (!"CellType" %in% colnames(tmp.seu[[]])) {
    stop("Please add CellType annotation in seurat object!")
  }
  object.assays <- Seurat:::FilterObjects(object = tmp.seu, 
                                          classes.keep = "Assay")
  if (use.impute == F){
    tmp.mat <- as.matrix(GetAssayData(object = tmp.seu,slot = "data", assay = "RNA"))
  }
  else{
    if (!"alra" %in% object.assays) {
      tmp.seu <- RunALRA(tmp.seu)
    }
    tmp.mat <- as.matrix(GetAssayData(object = tmp.seu,slot = "data", assay = "alra"))
  }
  
  tmp.df <- tmp.seu[[]]
  tmp.cell.id <- tmp.df %>% filter(CellType %in% tmp.ident) %>% 
    rownames()
  rtni <- tni.constructor(expData = tmp.mat[, tmp.cell.id], 
                          regulatoryElements = gene.TF)
  options(cluster = snow::makeCluster(spec = 10, "SOCK"))
  rtni <- tni.permutation(rtni, nPermutations = 1000, pValueCutoff = 0.01)
  rtni <- tni.bootstrap(rtni)
  stopCluster(getOption("cluster"))
  rtni <- tni.dpi.filter(rtni)
  rtni@summary
  rtna <- tni2tna.preprocess(rtni, phenotype = tmp.phenotype, 
                             hits = names(tmp.phenotype))
  rtna <- tna.mra(rtna)
  cat(paste0(tmp.cell.type, " TRN reconsturct end."), sep = "\n")
  return(rtna)
}




#####---------3.7.1 load data---------

#### gene expression data
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)

beforeEPI.cds <- readRDS(file = "res/R/early.scRNAseq.moncole_2022031116.rds")
beforeEPI.meta <- pData(beforeEPI.cds)

tmp.cell.id <- rownames(beforeEPI.meta)
beforeEPI.exp <- beforeEPI.exp[,tmp.cell.id]

#### TF
gene.TF <- read.delim(file = "database/Animal_TF_DB/Mus_musculus_TF.txt",stringsAsFactors = F) %>%
  pull(Symbol)
gene.TF.cofactor <- read.delim(file = "database/Animal_TF_DB/Mus_musculus_TF_cofactors.txt",stringsAsFactors = F) %>%
  pull(Symbol)
gene.TF <- intersect(gene.TF,rownames(beforeEPI.exp))

####------------3.7.2 Construct TNI object---------------

rtni <- tni.constructor(expData = as.matrix(beforeEPI.exp),
                        regulatoryElements = gene.TF)


# Please set nPermutations >= 1000
options(cluster=snow::makeCluster(spec=10, "SOCK"))
rtni <- tni.permutation(rtni,nPermutations = 1000,pValueCutoff = 1e-7)
rtni <- tni.bootstrap(rtni)
stopCluster(getOption("cluster"))

# Compute the DPI-filtered regulatory network
rtni <- tni.dpi.filter(rtni, eps = NA)
tni.regulon.summary(rtni)



saveRDS(rtni,file = "res/R/early.rtni_20220316.rds")


###------3.7.3 TRA reference----------


###load data
rtni <- readRDS(file = "res/R/early.rtni_20220316.rds")
tmp.cds <- readRDS(file = "res/R/early.scRNAseq.moncole_2022031116.rds")
head(fData(tmp.cds))
tmp.fData <- fData(tmp.cds) 


tmp.vst.variance <- tmp.fData$vst.variance
names(tmp.vst.variance) <- rownames(tmp.fData)
head(tmp.vst.variance)
tmp.vst.gene <- tmp.fData %>%
  filter(vst.variable==T) %>%
  rownames()

###

rtna <- tni2tna.preprocess(rtni, 
                           phenotype=tmp.vst.variance ,
                           hits=tmp.vst.gene)

?tna.mra
#run MRA analysis pipeline
rtna <- tna.mra(rtna)
?tna.mra
# options(cluster=snow::makeCluster(10, "SOCK"))
# rtna <- tna.gsea1(rtna, nPermutations=1000)
# stopCluster(getOption("cluster"))

# options(cluster=snow::makeCluster(10, "SOCK"))
# rtna <- tna.gsea2(rtna, nPermutations=1000)
# stopCluster(getOption("cluster"))

tmp.master.regulon <- tna.get(rtna,what = "mra")

saveRDS(tmp.master.regulon,
        file = myFileName(prefix = "res/R/putative_master_regulon",suffix = ".rds"))

saveRDS(rtna,
        file = "res/R/early.rtna_20220316.rds")

#### Illustrate cases

which(rownames(tmp.master.regulon) %in% 
        c("Pou5f1","Sox2","Klf4","Cdx2","Tead4"))

tmp.master.regulon[ c("Pou5f1","Zscan4f","Zscan4d","Sox2","Klf4","Cdx2","Tead4"),]

#####---------3.3------------


cat(paste0("prepare ligand and receptor genelist"),sep = "\n")
LRpairs <- LRpairs.df$lr_pair
Lgenelist <- LRpairs.df$ligand_gene_symbol
Rgenelist <- LRpairs.df$receptor_gene_symbol

### using data in RNA assay
tmp.data <-  GetAssayData(object = seu,slot = "data",assay = "RNA")
gene_symbols <- rownames(tmp.data)
l.remove <- setdiff(Lgenelist,gene_symbols)
r.remove <- setdiff(Rgenelist,gene_symbols)
index.remove <- c(which(Lgenelist %in% l.remove),which(Rgenelist %in% r.remove))
LRpairs <- LRpairs[-index.remove]
Lgenelist <- Lgenelist[-index.remove]
Rgenelist <- Rgenelist[-index.remove]


cat(paste0(tmp.ident.1,"-",tmp.ident.2," start:"),sep = "\n")
tmp.df <- data.frame(pseudotime = pseudotime(cds,reduction_method = "UMAP"),
                     stringsAsFactors = F)
tmp.df <- cbind(seu[[]],tmp.df)
### Check CellTypes
if (!c("CellType") %in% colnames(tmp.df)) {
  stop("Please add CellType annotation in seurat object!")
}

if (!c("orig.ident") %in% colnames(tmp.df)) {
  stop("Please add CellType annotation in seurat object!")
}

cat(paste0("scale pseudotime"),sep = "\n")
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

tmp.mat.1 <- tmp.data[,tmp.cell.meta.1$cell_id]
tmp.mat.2 <- tmp.data[,tmp.cell.meta.2$cell_id]

tmp.mat.pseudotime.1 <- tmp.cell.meta.1$pseudotime
tmp.mat.pseudotime.2 <- tmp.cell.meta.2$pseudotime

names(tmp.mat.pseudotime.1) <- tmp.cell.meta.1$cell_id
names(tmp.mat.pseudotime.2) <- tmp.cell.meta.2$cell_id

inter.tmp.mat.1 = cellAlign::interWeights(expDataBatch = tmp.mat.1,
                                          trajCond = tmp.mat.pseudotime.1,
                                          winSz = tmp.winsz,
                                          numPts = numPts)
inter.tmp.mat.2 = cellAlign::interWeights(expDataBatch = tmp.mat.2,
                                          trajCond = tmp.mat.pseudotime.2,
                                          winSz = tmp.winsz,
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

tmp.LR.list <- tmp.res.cor %>%
  filter(abs(SCC) > tmp.SCC.cutoff) %>%
  pull(LRpairs)

#### read master regulator analysis
tmp.TF.gene <- tmp.mra.res %>%
  filter(group == tmp.ident.2) %>%
  pull(Regulon)


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
    mutate(category=ifelse(LRtoTF < tmp.granger.cutoff | TFtoLR < tmp.granger.cutoff,"PASS","SKIP"))
  return(tmp.res.df)
})
plan("sequential")
tmp.ttt.res.df <- Reduce(rbind,tmp.ttt.res)
tmp.res.df <- tmp.ttt.res.df
return(tmp.res.df)








#####load LR pair
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_2022031117.rds")
Lgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgene.list <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
Lgene <- unique(Lgene.list)
Rgene <- unique(Rgene.list)
####test weather there exist Lgene.list == Rgene.list
sum(Lgene.list==Rgene.list)
####none, that's good!

####tmp.
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


#####--------3.2 calc cell stage correlation----------------------------------


tmp.L.exp <- beforeEPI.exp[Lgene.list,tmp.cell.id] 
tmp.R.exp <- beforeEPI.exp[Rgene.list,tmp.cell.id]
tmp.res.PCC <- sapply(1:nrow(tmp.L.exp), function(i) cor(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),
                                                         method = "pearson"))
tmp.res.SCC <- sapply(1:nrow(tmp.L.exp), function(i) cor(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),
                                                         method = "spearman"))


# tmp.res.granger.LtoR <- sapply(1:nrow(tmp.L.exp), function(i) {
#   tmp.res.res <-  tryCatch(
#     expr = {
#       tmp.res <- grangertest(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]))
#       return(tmp.res$`Pr(>F)`[2])},
#     error = function(e) {
#       return(1)
#     })
#   return(tmp.res.res)
# })
# 
# tmp.res.granger.RtoL <- sapply(1:nrow(tmp.L.exp), function(i) {
#   tmp.res.res <-  tryCatch(
#     expr = {
#       tmp.res <- grangertest(as.numeric(tmp.R.exp[i,]), as.numeric(tmp.L.exp[i,]))
#       return(tmp.res$`Pr(>F)`[2])},
#     error = function(e) {
#       return(1)
#     })
#   return(tmp.res.res)
# })

tmp.res.rate.L <- sapply(1:nrow(tmp.L.exp), function(i) sum(as.numeric(tmp.L.exp[i,])>0)/length(as.numeric(tmp.L.exp[i,])>0))
tmp.res.rate.R <- sapply(1:nrow(tmp.R.exp), function(i) sum(as.numeric(tmp.R.exp[i,])>0)/length(as.numeric(tmp.R.exp[i,])>0))

tmp.res.cor.1 <- data.frame(Lgene=Lgene.list,
                            Rgene=Rgene.list,
                            LRpairs=paste0(Lgene.list,"_",Rgene.list),
                            PCC=tmp.res.PCC,
                            SCC=tmp.res.SCC,
                            detection.rate.L=tmp.res.rate.L,
                            detection.rate.R=tmp.res.rate.R,
                            stringsAsFactors = F) %>%
  dplyr::filter(detection.rate.L > 0.05 & detection.rate.R > 0.05) %>%
  arrange(-PCC)

####------------3.3 identify eLR-------------

cutoff <- 0.1
tmp.var <- "PCC"
data.plot <- tmp.res.cor.1 %>%
  mutate(Rank=row_number(dplyr::desc(!!sym(tmp.var)))) %>%
  arrange(Rank) %>%
  mutate(group=ifelse((!!sym(tmp.var)) > cutoff,"p-eLR",ifelse((!!sym(tmp.var)) < -cutoff,"n-eLR","not sig"))) %>%
  mutate(group=factor(group,levels = c("p-eLR","not sig","n-eLR"))) %>%
  mutate(LRpairs=gsub(x = LRpairs,pattern = "_",replacement = "-"))

p <- ggplot(tmp.data.plot,aes(tmp.rank,gene.exp,group=group,col=group))+
  geom_point(alpha=0.3)+
  geom_smooth(se = F,method = "gam", formula = y ~ s(x, bs = "cs"),alpha=0.3)+
  xlab("rank")+
  ylab("log2(RPKM+1)")+
  ggtitle(label = paste0(tmp.Lgene,"-",tmp.Rgene,":","PCC=",signif(tmp.pcc,digits = 4)))+
  scale_color_manual(values = mypalette[1:2])+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  theme_cowplot(font_size = 28)+
  theme(legend.position = "top",legend.justification = "center",
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        plot.title = element_text(hjust = 0.5))

p
####----------3.4 link to putative TF--------------
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
  
  
  p <- ggplot(tmp.data.plot,aes(tmp.rank,gene.exp,group=group,col=group))+
    geom_point(alpha=0.3)+
    geom_smooth(se = F,method = "gam", formula = y ~ s(x, bs = "cs"),alpha=0.3)+
    xlab("rank by pseudotime")+
    ylab("log2(RPKM+1)")+
    ggtitle(label = paste0(tmp.Lgene,"-",tmp.Rgene,",",tmp.TF))+
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

####Gata6 positively regulate Fgfr2

tmp.TF.list <- c("Gata6","Sox2","Pou5f1","Nanog")
lapply(tmp.TF.list,FUN = function(ii){
  tmp.Lgene <- "Fgf4"
  tmp.Rgene <- "Fgfr1"
  tmp.TF <- ii
  p <- tmp.TF.IS.plot(tmp.Lgene = tmp.Lgene,tmp.Rgene = tmp.Rgene,tmp.TF = tmp.TF)
  myggsave(p = p,
           prefix = paste0("res/fig/Fig2_show_case_",tmp.Lgene,"_",tmp.Rgene,"_",tmp.TF),
           suffix = ".png",
           width = 8,
           height = 8,
           dpi=350)
})

tmp.TF.list <- c("Gata6","Sox2","Pou5f1","Nanog")
lapply(tmp.TF.list,FUN = function(ii){
  tmp.Lgene <- "Fgf4"
  tmp.Rgene <- "Fgfr2"
  tmp.TF <- ii
  p <- tmp.TF.IS.plot(tmp.Lgene = tmp.Lgene,tmp.Rgene = tmp.Rgene,tmp.TF = tmp.TF)
  myggsave(p = p,
           prefix = paste0("res/fig/Fig2_show_case_",tmp.Lgene,"_",tmp.Rgene,"_",tmp.TF),
           suffix = ".png",
           width = 8,
           height = 8,
           dpi=350)
})




tmp.Lgene <- "Fgf4"
tmp.Rgene <- "Fgfr1"
tmp.TF <- "Gata6"
p <- tmp.TF.IS.plot(tmp.Lgene = tmp.Lgene,tmp.Rgene = tmp.Rgene,tmp.TF = tmp.TF)
p
myggsave(p = p,
         prefix = paste0("res/fig/Fig2_show_case",tmp.Lgene,"_",tmp.Rgene,"_",tmp.TF),
         suffix = ".png",
         width = 8,
         height = 8,
         dpi=350)


tmp.Lgene <- "Fgf4"
tmp.Rgene <- "Fgfr1"
tmp.TF <- "Sox2"
p <- tmp.TF.IS.plot(tmp.Lgene = tmp.Lgene,tmp.Rgene = tmp.Rgene,tmp.TF = tmp.TF)
p
myggsave(p = p,
         prefix = paste0("res/fig/Fig2_show_case",tmp.Lgene,"_",tmp.Rgene,"_",tmp.TF),
         suffix = ".png",
         width = 8,
         height = 8,
         dpi=350)

tmp.Lgene <- "Fgf4"
tmp.Rgene <- "Fgfr2"
tmp.TF <- "Sox2"
p <- tmp.TF.IS.plot(tmp.Lgene = tmp.Lgene,tmp.Rgene = tmp.Rgene,tmp.TF = tmp.TF)
p
myggsave(p = p,
         prefix = paste0("res/fig/Fig2_show_case",tmp.Lgene,"_",tmp.Rgene,"_",tmp.TF),
         suffix = ".png",
         width = 8,
         height = 8,
         dpi=350)



tmp.Lgene <- "Fgf4"
tmp.Rgene <- "Fgfr2"
tmp.TF <- "Gata6"

p <- tmp.TF.IS.plot(tmp.Lgene = tmp.Lgene,tmp.Rgene = tmp.Rgene,tmp.TF = tmp.TF)
p
myggsave(p = p,
         prefix = paste0("res/fig/Fig2_show_case",tmp.Lgene,"_",tmp.Rgene,"_",tmp.TF),
         suffix = ".png",
         width = 8,
         height = 8,
         dpi=350)

####Nanog inhibit Fgfr1
tmp.Lgene <- "Fgf4"
tmp.Rgene <- "Fgfr1"
tmp.TF <- "Nanog"
p <- tmp.TF.IS.plot(tmp.Lgene = tmp.Lgene,tmp.Rgene = tmp.Rgene,tmp.TF = tmp.TF)
myggsave(p = p,
         prefix = paste0("res/fig/Fig2_show_case",tmp.Lgene,"_",tmp.Rgene,"_",tmp.TF),
         suffix = ".png",
         width = 8,
         height = 8,
         dpi=350)


tmp.Lgene <- "Fgf4"
tmp.Rgene <- "Fgfr2"
tmp.TF <- "Nanog"
p <- tmp.TF.IS.plot(tmp.Lgene = tmp.Lgene,tmp.Rgene = tmp.Rgene,tmp.TF = tmp.TF)
p
myggsave(p = p,
         prefix = paste0("res/fig/Fig2_show_case",tmp.Lgene,"_",tmp.Rgene,"_",tmp.TF),
         suffix = ".png",
         width = 8,
         height = 8,
         dpi=350)




#####---------3.5 correlation with TF and regulon----------------

gene.TF <- read.delim(file = "database/Animal_TF_DB/Mus_musculus_TF.txt",stringsAsFactors = F) %>%
  pull(Symbol)
gene.TF.cofactor <- read.delim(file = "database/Animal_TF_DB/Mus_musculus_TF_cofactors.txt",stringsAsFactors = F) %>%
  pull(Symbol)


tmp.Lgene <- "Fgf4"
tmp.Rgene <- "Fgfr2"
tmp.TF <- "Gata6"

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


tmp.data.plot <- rbind(tmp.data.plot.3,tmp.data.plot.4)


tmp.data.plot$group <- factor(tmp.data.plot$group,
                              levels = c("TF","IS"))

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
  scale_color_manual(values = tmp.color[3:4],name=NULL)+
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
p_res

myggsave(p_res,prefix = "res/fig/fig2_Fgf4_Fgfr2_Gata6_causal",suffix = ".png",width = 8,height = 8,dpi = 400)
myggsave(p_res,prefix = "res/fig/fig2_Fgf4_Fgfr2_Gata6_causal",suffix = ".pdf",width = 8,height = 8,dpi = 400)



#####---------3.6 get regulon by RTN------------

#####---------3.6.1 load data---------
# seu <- readRDS(file = "res/R/early.scRNAseq.seurat_2022031015.rds")
# seu <- RunALRA(seu)

#### gene expression data
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")
beforeEPI.exp <- myRemoveNA(beforeEPI.exp)
beforeEPI.exp <- myNormalize(beforeEPI.exp)

beforeEPI.cds <- readRDS(file = "res/R/early.scRNAseq.moncole_2022031116.rds")
beforeEPI.meta <- pData(beforeEPI.cds)

tmp.cell.id <- rownames(beforeEPI.meta)
beforeEPI.exp <- beforeEPI.exp[,tmp.cell.id]

#### TF
gene.TF <- read.delim(file = "database/Animal_TF_DB/Mus_musculus_TF.txt",stringsAsFactors = F) %>%
  pull(Symbol)
gene.TF.cofactor <- read.delim(file = "database/Animal_TF_DB/Mus_musculus_TF_cofactors.txt",stringsAsFactors = F) %>%
  pull(Symbol)
gene.TF <- intersect(gene.TF,rownames(beforeEPI.exp))





####------------3.6.2 Construct TNI object---------------

rtni <- tni.constructor(expData = as.matrix(beforeEPI.exp),
                        regulatoryElements = gene.TF)


# Please set nPermutations >= 1000
options(cluster=snow::makeCluster(spec=10, "SOCK"))
rtni <- tni.permutation(rtni,nPermutations = 1000,pValueCutoff = 1e-7)
rtni <- tni.bootstrap(rtni)
stopCluster(getOption("cluster"))

# Compute the DPI-filtered regulatory network
rtni <- tni.dpi.filter(rtni, eps = NA)
tni.regulon.summary(rtni)



saveRDS(rtni,file = "res/R/early.rtni_20220316.rds")


###------3.6.3 TRA reference----------


###load data
rtni <- readRDS(file = "res/R/early.rtni_20220316.rds")
tmp.cds <- readRDS(file = "res/R/early.scRNAseq.moncole_2022031116.rds")
head(fData(tmp.cds))
tmp.fData <- fData(tmp.cds) 


tmp.vst.variance <- tmp.fData$vst.variance
names(tmp.vst.variance) <- rownames(tmp.fData)
head(tmp.vst.variance)
tmp.vst.gene <- tmp.fData %>%
  filter(vst.variable==T) %>%
  rownames()

###

rtna <- tni2tna.preprocess(rtni, 
                           phenotype=tmp.vst.variance ,
                           hits=tmp.vst.gene)

?tna.mra
#run MRA analysis pipeline
rtna <- tna.mra(rtna)

# options(cluster=snow::makeCluster(10, "SOCK"))
# rtna <- tna.gsea1(rtna, nPermutations=1000)
# stopCluster(getOption("cluster"))

# options(cluster=snow::makeCluster(10, "SOCK"))
# rtna <- tna.gsea2(rtna, nPermutations=1000)
# stopCluster(getOption("cluster"))

tmp.master.regulon <- tna.get(rtna,what = "mra")

saveRDS(tmp.master.regulon,
        file = myFileName(prefix = "res/R/putative_master_regulon",suffix = ".rds"))

saveRDS(rtna,
        file = "res/R/early.rtna_20220316.rds")

#### Illustrate cases

which(rownames(tmp.master.regulon) %in% 
        c("Pou5f1","Sox2","Klf4","Cdx2","Tead4"))

tmp.master.regulon[ c("Pou5f1","Zscan4f","Zscan4d","Sox2","Klf4","Cdx2","Tead4"),]




####---------3.7 TimeTalk prototype version-------------------

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


####---------3.7.1 filter to get the master regulator of the result-------------

tmp.master.regulon <- readRDS(file = "res/R/putative_master_regulon_2022032109.rds")
head(tmp.master.regulon)
### add mean and sd to data frame
tmp.TF.mean <- rowMeans(beforeEPI.exp[tmp.master.regulon$Regulon,])
tmp.TF.sd <- apply(beforeEPI.exp[tmp.master.regulon$Regulon,],MARGIN = 1,sd)
tmp.master.regulon$mean.exp <- tmp.TF.mean
tmp.master.regulon$sd <- tmp.TF.sd
tmp.df <- tmp.master.regulon %>%
  filter(mean.exp > 1)

####As a by product we need to analysis master regulator


####---------3.7.2 plot the master regulator heatmap -----------

####---------3.7.2.1 prepare data------------
tmp.min <- -3
tmp.max <- 3
tmp.mat <- beforeEPI.exp[tmp.df$Regulon,tmp.cell.id]
tmp.mat <- pheatmap:::scale_rows(tmp.mat)
tmp.mat[tmp.mat > tmp.max] <- tmp.max
tmp.mat[tmp.mat < tmp.min] <- -tmp.min

####-------3.7.2.2 draw draft pheatmap----------

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
tmp.cluster.color <- monet.sunumbrella.women(num_clusters)
names(tmp.cluster.color) <- 1:num_clusters

colnames(annotation_row) <- "tTFs"

ann_colors = list(
  Pseudotime = col.spectral(100),
  Stage = tmp.Stage.color,
  tTFs = tmp.cluster.color
)

#####------3.7.2.3 re annotation levels-----------
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


annotation_row$tTFs <- plyr::mapvalues(annotation_row$tTFs,
                                          from = c(3,1,6,4,2,5),
                                          to = 1:num_clusters)
annotation_row$tTFs <- factor(annotation_row$tTFs,
                                 levels = 1:num_clusters)
annotation_row.df <- annotation_row %>%
  rownames_to_column("gene")
annotation_row.df <- annotation_row.df[tmp.order,] %>%
  group_by(tTFs) %>%
  mutate(tmp.rank = row_number()) %>%
  ungroup() %>%
  arrange(tTFs,tmp.rank) %>%
  dplyr::select(gene,tTFs) %>%
  column_to_rownames("gene")
tmp.mat.plot <- tmp.mat[rownames(annotation_row.df),]

early_embryo_TF <- annotation_row.df %>%
  rownames_to_column("gene")

early_embryo_TF %>%
  filter(tTFs == 4)


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

saveRDS(annotation_row.df,file = myFileName(prefix = "res/R/early_embryo_core_TF_annotation_row",suffix = ".rds"))

####-----------3.7.2.4 show cluster TF temporal pattern -----------------


annotation_row.df <- readRDS(file = "res/R/early_embryo_core_TF_annotation_row_2022032521.rds")
tmp.TF.meta <- annotation_row.df %>%
  rownames_to_column("gene")
num_clusters <- 6
tmp.cluster.color <- pal_npg()(num_clusters)

head(tmp.meta.df.col)

data.plot <- tmp.meta.df.col %>%
  mutate(tmp.rank = row_number(),
         tmp.rect = 1)

head(data.plot)
p_bottom <- ggplot(data.plot,aes(x = tmp.rank,
                     y = tmp.rect,fill=Stage))+
  geom_tile()+
  scale_fill_manual(values = tmp.Stage.color)+
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0))+
  ylab(NULL)+
  xlab("rank by pseudotime")+
  theme_cowplot(font_size = 28,rel_small = 8/14)+
  NoLegend()+
  theme(axis.ticks = element_line(size = 2),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 22),
        axis.line.y = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.y = element_blank())
p_bottom



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
                color = tmp.cluster.color[ii])+
    xlab("rank by pesudotime")+
    ylab("zscore")+
    scale_x_continuous(expand = c(0,0))+
    ggtitle(paste0("C",ii))+
    theme_cowplot(font_size = 28,rel_small = 8/14)+
    theme(axis.line = element_line(size = 2),
          legend.position = "top",
          legend.justification = "center",
          axis.ticks = element_line(size = 2),
          plot.title = element_text(hjust = 0.5))
  p <- p1+xlab(NULL)+theme(axis.ticks.x = element_blank(),
                           axis.text.x = element_blank(),
                      axis.line.x = element_blank()) 
  p_res <- wrap_plots(p,p_bottom,ncol = 1,heights = c(9,1)) &
    theme(plot.margin = margin(t = 0,  
                               r = 0,  
                               b = 0,  
                               l = 0))

  p_res
  return(p_res)
})
p <- wrap_plots(plot.list,nrow = 2)
myggsave(p,prefix = "res/fig/fig2_TF_cluster",suffix = ".png",width = 12,height = 8,dpi=350)
myggsave(p,prefix = "res/fig/fig2_TF_cluster",suffix = ".pdf",width = 12,height = 8,dpi=350)

lapply(1:6,FUN = function(ii){
  p <- plot.list[[ii]]
  myggsave(p,prefix = paste0("res/fig/fig2_TF_cluster",ii),suffix = ".png",width = 6,height = 6,dpi=350)
  myggsave(p,prefix = paste0("res/fig/fig2_TF_cluster",ii),suffix = ".pdf",width = 6,height = 6,dpi=350)
  
})

p <- plot.list[[5]]
myggsave(p,prefix = "res/fig/fig2_TF_cluster5",suffix = ".png",width = 6,height = 6,dpi=350)
myggsave(p,prefix = "res/fig/fig2_TF_cluster5",suffix = ".pdf",width = 6,height = 6,dpi=350)

p <- plot.list[[1]]
myggsave(p,prefix = "res/fig/fig2_TF_cluster1",suffix = ".png",width = 6,height = 6,dpi=350)
myggsave(p,prefix = "res/fig/fig2_TF_cluster1",suffix = ".pdf",width = 6,height = 6,dpi=350)

dev.off()

#####---------3.8 perform TimeTalk to get eLR ----------------------- 


####----------3.8.1 load data------------------

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


#####-------3.8.2 test cluster  causal test------------

ii <- 5
tmp.TF.gene <- tmp.TF.meta %>%
  filter(Cluster == ii) %>%
  pull(gene)
tmp.TF.exp.score <- colMeans(beforeEPI.exp[tmp.TF.gene,tmp.cell.id])


tmp.Lgene <- "Fgf4"
tmp.Rgene <- "Fgfr2"
tmp.TF <- "Gata6"

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
  mutate(group = "C5")

tmp.data.plot.3$gene.exp <- tmp.TF.exp.score

tmp.data.plot.4 <- tmp.data.plot.2 %>%
  mutate(gene.exp = sqrt(tmp.data.plot.1$gene.exp * tmp.data.plot.2$gene.exp)) %>%
  mutate(group = "IS")


tmp.data.plot <- rbind(tmp.data.plot.3,tmp.data.plot.4)


tmp.data.plot$group <- factor(tmp.data.plot$group,
                              levels = c("C5","IS"))

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
  ggtitle(label = paste0(tmp.Lgene,"-",tmp.Rgene,",","master TFs in C5","\n",
                         "p_TFtoLR=",signif(tmp.pvalue.TFtoLR,digits = 4),"\n",
                         "p_LRtoTF=",signif(tmp.pvalue.LRtoTF,digits = 4)))+
  scale_color_manual(values = c("#f39b7f","#631879FF"),name=NULL)+
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
p_res

myggsave(p_res,prefix = "res/fig/fig2_Fgf4_Fgfr2_C5_causal",suffix = ".png",width = 8,height = 8,dpi = 400)
myggsave(p_res,prefix = "res/fig/fig2_Fgf4_Fgfr2_C5_causal",suffix = ".pdf",width = 8,height = 8,dpi = 400)


lapply(1:6, FUN = function(ii){
  
  cat(ii,sep = "\n")
  tmp.TF.gene <- tmp.TF.meta %>%
    filter(Cluster == ii) %>%
    pull(gene)
  tmp.TF.exp.score <- colMeans(beforeEPI.exp[tmp.TF.gene,tmp.cell.id])
  
  tmp.name <- paste0("C",ii)
    
  tmp.Lgene <- "Fgf4"
  tmp.Rgene <- "Fgfr2"
  
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
    mutate(group = tmp.name)
  
  tmp.data.plot.3$gene.exp <- tmp.TF.exp.score
  
  tmp.data.plot.4 <- tmp.data.plot.2 %>%
    mutate(gene.exp = sqrt(tmp.data.plot.1$gene.exp * tmp.data.plot.2$gene.exp)) %>%
    mutate(group = "IS")
  
  
  tmp.data.plot <- rbind(tmp.data.plot.3,tmp.data.plot.4)
  
  
  tmp.data.plot$group <- factor(tmp.data.plot$group,
                                levels = c(tmp.name,"IS"))
  
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
    ggtitle(label = paste0(tmp.Lgene,"-",tmp.Rgene,",","master TFs in ",tmp.name,"\n",
                           "p_TFtoLR=",signif(tmp.pvalue.TFtoLR,digits = 4),"\n",
                           "p_LRtoTF=",signif(tmp.pvalue.LRtoTF,digits = 4)))+
    scale_color_manual(values = c(tmp.cluster.color[ii],"#631879FF"),name=NULL)+
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

  myggsave(p_res,prefix = paste0("res/fig/fig2_Fgf4_Fgfr2_",tmp.name,"_causal"),suffix = ".png",width = 8,height = 8,dpi = 400)
  myggsave(p_res,prefix = paste0("res/fig/fig2_Fgf4_Fgfr2_",tmp.name,"_causal"),suffix = ".pdf",width = 8,height = 8,dpi = 400)
  
})

lapply(1:6, FUN = function(ii){
  
  cat(ii,sep = "\n")
  tmp.TF.gene <- tmp.TF.meta %>%
    filter(Cluster == ii) %>%
    pull(gene)
  tmp.TF.exp.score <- colMeans(beforeEPI.exp[tmp.TF.gene,tmp.cell.id])
  
  tmp.name <- paste0("C",ii)
  
  tmp.Lgene <- "Bmp4"
  tmp.Rgene <- "Bmpr2"
  
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
    mutate(group = tmp.name)
  
  tmp.data.plot.3$gene.exp <- tmp.TF.exp.score
  
  tmp.data.plot.4 <- tmp.data.plot.2 %>%
    mutate(gene.exp = sqrt(tmp.data.plot.1$gene.exp * tmp.data.plot.2$gene.exp)) %>%
    mutate(group = "IS")
  
  
  tmp.data.plot <- rbind(tmp.data.plot.3,tmp.data.plot.4)
  
  
  tmp.data.plot$group <- factor(tmp.data.plot$group,
                                levels = c(tmp.name,"IS"))
  
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
    ggtitle(label = paste0(tmp.Lgene,"-",tmp.Rgene,",","master TFs in ",tmp.name,"\n",
                           "p_TFtoLR=",signif(tmp.pvalue.TFtoLR,digits = 4),"\n",
                           "p_LRtoTF=",signif(tmp.pvalue.LRtoTF,digits = 4)))+
    scale_color_manual(values = c(tmp.cluster.color[ii],"#631879FF"),name=NULL)+
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
  
  myggsave(p_res,prefix = paste0("res/fig/fig2_",tmp.Lgene,"_",tmp.Rgene,"_",tmp.name,"_causal"),suffix = ".png",width = 8,height = 8,dpi = 400)
  myggsave(p_res,prefix = paste0("res/fig/fig2_",tmp.Lgene,"_",tmp.Rgene,"_",tmp.name,"_causal"),suffix = ".pdf",width = 8,height = 8,dpi = 400)
  
})




####-------2.2 perform granger causal test -----------

#### Actualy the following code is the TimeTalk version 0.0.1


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

tmp.res.eLR.df <- readRDS(file = "res/R/early_embryo_eLR_TF_interaction_information_2022060620.rds")



#####---------4. plot for illustrate -------------


#####---------4.1 show how get the high correlated eLR-----------


### This function were adapted from ROSE
### https://bitbucket.org/young_computation/rose/src/master/ROSE_callSuper.R

calculate_cutoff <- function(inputVector, drawPlot=TRUE,...){
  
  inputVector <- sort(inputVector)
  inputVector[inputVector<0]<-0 #set those regions with more control than ranking equal to zero
  slope <- (max(inputVector)-min(inputVector))/length(inputVector) #This is the slope of the line we want to slide. This is the diagonal.
  # slope <- (max(inputVector)-inputVector)/length(inputVector)
  # tmp.idx <- which(abs(inputVector - (max(inputVector)-0.1)) < 0.0001)
  # inputVector[length(inputVector) - tmp.idx]
  # inputVector[1]
  # slope <- max(inputVector)
  
  xPt <- floor(optimize(numPts_below_line,lower=1,upper=length(inputVector),myVector= inputVector,slope=slope)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  y_cutoff <- inputVector[xPt] #The y-value at this x point. This is our cutoff.
  
  # ## second slope add by wlt to get stringent cutoff
  # slope_final <- (max(inputVector) - y_cutoff)/(length(inputVector)-xPt)
  # xPt_final <- floor(optimize(numPts_below_line,
  #                             lower=xPt,
  #                             upper=length(inputVector),
  #                             myVector= inputVector,
  #                             slope=slope_final)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  # y_cutoff_final <- inputVector[xPt_final] #The y-value at this x point. This is our cutoff.
  
  # ## second slope add by wlt to get soft cutoff
  # slope_2 <- (y_cutoff - min(inputVector))/xPt
  # xPt_2 <- floor(optimize(numPts_below_line,
  #                             lower=1,
  #                             upper=xPt,
  #                             myVector= inputVector,
  #                             slope=slope_2)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  # y_cutoff_2 <- inputVector[xPt_2] #The y-value at this x point. This is our cutoff.
  
  # ## third slope add by wlt to get soft cutoff
  # slope_final <- (y_cutoff_2 - min(inputVector))/xPt_2
  # xPt_final <- floor(optimize(numPts_below_line,
  #                             lower=1,
  #                             upper=xPt_2,
  #                             myVector= inputVector,
  #                             slope=slope_final)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  # y_cutoff_final <- inputVector[xPt_final] #The y-value at this x point. This is our cutoff.
  
  ## second slope add by wlt to get soft cutoff
  slope_final <- (y_cutoff - min(inputVector))/xPt
  xPt_final <- floor(optimize(numPts_below_line,
                          lower=1,
                          upper=xPt,
                          myVector= inputVector,
                          slope=slope_final)$minimum) #Find the x-axis point where a line passing through that point has the minimum number of points below it. (ie. tangent)
  y_cutoff_final <- inputVector[xPt_final] #The y-value at this x point. This is our cutoff.
  
  
  
  if(drawPlot){  #if TRUE, draw the plot
    
    tmp.data.plot <- data.frame(rank = 1:length(inputVector),
                                PCC = inputVector,
                                stringsAsFactors = F)
    
    b <- y_cutoff-(slope* xPt)
    b_final <- y_cutoff_final-(slope_final * xPt_final)
    
    p <- ggplot(tmp.data.plot,aes(rank,PCC))+
      geom_line(size = 1)+
      geom_point(x = 1,y = inputVector[1],color = "skyblue",size = 3)+
      geom_point(x = length(inputVector),y = inputVector[length(inputVector)],color = "skyblue",size = 3 )+
      geom_path(data = data.frame(x = c(1,length(inputVector)),
                                  y = c(inputVector[1],inputVector[length(inputVector)]),
                                  stringsAsFactors = F),
                mapping = aes(x,y),
                size = 1,color = "skyblue")+
      geom_path(data = data.frame(x = c(1,xPt),
                                  y = c(inputVector[1],inputVector[xPt]),
                                  stringsAsFactors = F),
                mapping = aes(x,y),
                size = 1,color = "skyblue")+
      geom_abline(intercept = b,slope = slope,color = "skyblue",size = 1)+
      geom_abline(intercept = b_final,slope = slope_final,color = "skyblue",size = 1)+
      #geom_hline(yintercept = y_cutoff_final,color = "skyblue",size = 1,linetype = "dashed")+
      geom_point(x = xPt,y = y_cutoff,color = "skyblue",size = 3)+
      geom_point(x = xPt_final,y = y_cutoff_final,color = "skyblue",size = 3)+
      geom_text(x = 1+30,
                y = inputVector[1]+0.03,
                label = "A",
                size = 8)+
      geom_text(x = length(inputVector),
                y = inputVector[length(inputVector)]+0.02,
                label = "B",
                size = 8)+
      geom_text(x = xPt,
                y = y_cutoff-0.03,
                label = "C",
                size = 8)+
      geom_text(x = xPt_final+100,
                y = y_cutoff_final-0.03,
                label = paste0("D","(","x=",xPt_final,",","y=",round(y_cutoff_final,digits = 4),")"),
                size = 8)+
      ggtitle(label = paste0("cutoff is ",signif(y_cutoff_final,digits = 3)))+
      # scale_x_continuous(limits = c(0,700),expand = c(0,0))+
      # scale_y_continuous(limits = c(0,1),expand = c(0,0))+
      theme_cowplot(font_size = 20)+
      theme(plot.title = element_text(hjust = 0.5))

    # plot(1:length(inputVector), inputVector,type="l",...)
    # 
    # b <- y_cutoff-(slope* xPt)
    # abline(v= xPt,h= y_cutoff,lty=2,col=8)
    # points(xPt,y_cutoff,pch=16,cex=0.9,col=2)
    # abline(coef=c(b,slope),col=2)
    # 
    # b_final <- y_cutoff_final-(slope_final * xPt_final)
    # abline(v= xPt_final,h= y_cutoff_final,lty=2,col=8)
    # points(xPt_final,y_cutoff_final,pch=16,cex=0.9,col=2)
    # abline(coef=c(b_final,slope_final),col=2)
    # 
    # # title(paste("x=",xPt,"\ny=",signif(y_cutoff,3),"\nFold over Median=",
    # #             signif(y_cutoff/median(inputVector),3),"x\nFold over Mean=",
    # #             signif(y_cutoff/mean(inputVector),3),"x",sep=""))
    # 
    # title(paste("x=",xPt,"\ny=",signif(y_cutoff_final,3),"\nFold over Median=",
    #             signif(y_cutoff_final/median(inputVector),3),"x\nFold over Mean=",
    #             signif(y_cutoff_final/mean(inputVector),3),"x",sep=""))
    
    #axis(1,sum(inputVector==0),sum(inputVector==0),col.axis="pink",col="pink") #Number of regions with zero signal
  }
  # return(list(absolute=y_cutoff,
  #             overMedian=y_cutoff/median(inputVector),
  #             overMean=y_cutoff/mean(inputVector)))
  return(list(absolute=y_cutoff_final,
              overMedian=y_cutoff_final/median(inputVector),
              overMean=y_cutoff_final/mean(inputVector),
              plot = p))
}


numPts_below_line <- function(myVector,slope,x){
  yPt <- myVector[x]
  b <- yPt-(slope*x)
  xPts <- 1:length(myVector)
  return(sum(myVector<=(xPts*slope+b)))
}


cutoff = 0.1

data.plot <- readRDS(file = "res/R/putative_eLR_pairs_1094_2022032109.rds")
data.plot <- data.plot %>%
  mutate(group="not sig") %>%
  mutate(group=ifelse(abs(PCC)> cutoff,"sig",group)) %>%
  mutate(group=factor(group,levels = c("sig","not sig"))) %>%
  arrange(-PCC) %>%
  mutate(Rank=row_number())
tmp.var = "PCC"


### This step to explain why we choose to use 0.1 as cutoff
### negative cutoff
tmp.vec <- data.plot %>%
  dplyr::filter(PCC > 0) %>%
  pull(PCC)
tmp.res.1 <- calculate_cutoff(abs(tmp.vec))
tmp.res.1$plot
ggsave(filename = myFileName(prefix = "res/fig/fig2_PCC_arrange_PCC_for_revision_p_eLR",
                             suffix = ".jpg"),
       width = 8,height = 6,dpi = 350)
ggsave(filename = myFileName(prefix = "res/fig/fig2_PCC_arrange_PCC_for_revision_p_eLR",
                             suffix = ".pdf"),
       width = 8,height = 6)

### negative cutoff
tmp.vec <- data.plot %>%
  dplyr::filter(PCC < 0) %>%
  pull(PCC)
tmp.res.2 <- calculate_cutoff(abs(tmp.vec))
tmp.res.2$plot+
  ylab("-PCC")
ggsave(filename = myFileName(prefix = "res/fig/fig2_PCC_arrange_PCC_for_revision_n_eLR",
                             suffix = ".jpg"),
       width = 8,height = 6,dpi = 350)
ggsave(filename = myFileName(prefix = "res/fig/fig2_PCC_arrange_PCC_for_revision_n_eLR",
                             suffix = ".pdf"),
       width = 8,height = 6)


# the result is 0.1
round(max(c(tmp.res.1$absolute,tmp.res.2$absolute)),digits = 2)




p <- ggplot(data.plot,aes(Rank,(!!sym(tmp.var)),
                          color=group))+
  geom_point(alpha=0.8,shape=16,size=6)+
  scale_color_manual(values = c("#c00000","grey"),name=NULL,label=c("co-vary","not"))+
  scale_x_continuous(breaks = seq(0,1200,100),labels = seq(0,1200,100))+
  scale_y_continuous(breaks = seq(-1,1,0.2),labels = seq(-1,1,0.2),limits = c(-1,1))+
  geom_hline(yintercept = cutoff,linetype="dotdash")+
  geom_hline(yintercept = -cutoff,linetype="dotdash")+
  theme_cowplot(font_size = 30)+
  theme(axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        axis.text.x = element_text(angle = 45,vjust = 0.5),
        legend.position = "top",
        legend.justification = "center")
p

head(data.plot)
nrow(data.plot)



# ggsave(p,filename = myFileName(prefix = "res/fig/Fig3_eLR_show_get_high_cor_eLR",suffix = ".png"),
#        width = 6,height = 6,dpi = 350)
# ggsave(p,filename = myFileName(prefix = "res/fig/Fig3_eLR_show_get_high_cor_eLR",suffix = ".pdf"),
#        width = 6,height = 6,dpi = 350)
myggsave(p = p,prefix = "res/fig/Fig3_eLR_show_get_high_cor_eLR",suffix = ".png",width = 6,height = 6,dpi=350)
myggsave(p = p,prefix = "res/fig/Fig3_eLR_show_get_high_cor_eLR",suffix = ".pdf",width = 6,height = 6,dpi=350)
#####---------4.2  show how get regulate TF-----------

tmp.color <- c("brown3","navy","#f2be58","grey")
names(tmp.color) <-  c("forward","backward","feedback","background")


data.plot <- readRDS(file = "res/R/putative_eLR_pairs_1094_2022032109.rds")
data.plot <- data.plot %>%
  filter(abs(PCC) > 0.1) %>%
  #filter(LRtoTF < 0.01 | TFtoLR < 0.01) %>%
  mutate(group = "background") %>%
  mutate(group = ifelse(LRtoTF < 0.01 & TFtoLR > 0.01,"forward",group)) %>%
  mutate(group = ifelse(LRtoTF > 0.01 & TFtoLR < 0.01,"backward",group)) %>%
  mutate(group = ifelse(LRtoTF < 0.01 & TFtoLR < 0.01,"feedback",group)) %>%
  mutate(group = factor(group,levels = c("forward","backward","feedback","background"))) %>%
  mutate(LRtoTF=-log10(LRtoTF)) %>%
  mutate(TFtoLR=-log10(TFtoLR))


p <- ggplot(data.plot,aes(LRtoTF,TFtoLR,color=group))+
  geom_point(size=3,shape=16,alpha=0.8)+
  geom_hline(yintercept = 2,linetype="dotdash")+
  geom_vline(xintercept = 2,linetype="dotdash")+
  scale_color_manual(values = tmp.color,
                     name=NULL,
                     guide=guide_legend(nrow = 2))+
  scale_x_continuous(breaks = seq(0,16,by=2))+
  scale_y_continuous(breaks = seq(0,16,by=2))+
  xlab(expression(-log[10](P[LRtoTF])))+
  ylab(expression(-log[10](P[TFtoLR])))+
  theme_cowplot(font_size = 28)+
  theme(axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2),
        legend.position = "top",
        legend.justification = "center")


#?geom_point
p
# ggsave(p,filename = myFileName(prefix = "res/fig/Fig3_eLR_show_get_eLR",suffix = ".png"),
#        width = 6,height = 6,dpi = 350)
# ggsave(p,filename = myFileName(prefix = "res/fig/Fig3_eLR_show_get_eLR",suffix = ".pdf"),
#        width = 6,height = 6,dpi = 350)
myggsave(p = p,prefix = "res/fig/Fig3_eLR_show_get_eLR",suffix = ".png",width = 6,height = 6,dpi=350)
myggsave(p = p,prefix = "res/fig/Fig3_eLR_show_get_eLR",suffix = ".pdf",width = 6,height = 6,dpi=350)









