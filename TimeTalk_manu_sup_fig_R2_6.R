###--------0. Library-------------
library(tidyverse)
library(ggplot2)
library(GenomicFeatures)
library(GenomeInfoDb)
source("code/myUtils.R")
library(ggsignif)
###-------1. load DDD database and perform some analysis----------------

###--------1.1 investigate DDD database----------

head(tmp.data)
table(tmp.data$organ.specificity.list)


tmp.file <- "database/DDG2P_18_5_2023.csv.gz"
DDD_database <- read.csv(file = tmp.file,stringsAsFactors = F)

tmp.organ <- Reduce(f = c,
                    strsplit(x = tmp.data$organ.specificity.list,split = ";"))
table(tmp.organ)

table(tmp.data$disease.name)


###-------1.2 convert the human symbol list------------
tmp.path <- "database/HOM_MouseHumanSequence.rpt"
HOM_database <- read.delim(file = tmp.path,stringsAsFactors = F)

tmp.file <- "database/DDG2P_18_5_2023.csv.gz"
DDD_database <- read.csv(file = tmp.file,stringsAsFactors = F)
tmp.gene <- DDD_database$gene.symbol

tmp.gene <- intersect(tmp.gene,
                      unique(HOM_database$Symbol))

tmp.res.list <- lapply(1:length(tmp.gene),
                       FUN = function(ii){
                         cat(ii,sep = "\n")
                         x <- tmp.gene[ii]
                         
                         tmp.id <- HOM_database %>%
                           dplyr::filter(Symbol == x) %>%
                           pull(HomoloGene.ID)
                         
                         y <- HOM_database %>%
                           dplyr::filter(NCBI.Taxon.ID == 10090 & 
                                           HomoloGene.ID == tmp.id) %>%
                           pull(Symbol)
                         
                         if (is_empty(y)) {
                           y <- NA
                         }
                         
                         tmp.res <- data.frame(H_symbol = x,
                                               M_symbol = y,
                                               stringsAsFactors = F)
                         return(tmp.res)
                         
                       })

tmp.df <- Reduce(rbind,tmp.res.list)

tmp_DDD_use <- merge(DDD_database,tmp.df,
                     by.x = "gene.symbol",
                     by.y = "H_symbol") %>%
  dplyr::filter(is.na(M_symbol) == F)


###-------1.3 Perform Enrichment----------

mm9KG_txdb <- makeTxDbFromGFF(file = "/mnt/e/project/reference/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf")
mm9.gene <- GenomicFeatures::genes(mm9KG_txdb)
mm9.gene.list <- mm9.gene$gene_id
### load DDD gene 
?genes

DDD.gene <- tmp_DDD_use$M_symbol
DDD.gene <- intersect(DDD.gene,mm9.gene.list)

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
M <- length(DDD.gene)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,DDD.gene))
n <- length(LRgene)
x <- M/N
y <- k/n
p <- phyper(k-1,M,N-M,n,lower.tail = F)
tmp.df.1 <- data.frame(value = c(x,y),group=c("all","LR"),stringsAsFactors = F)
pvalue.list <- p

####-------1.4 eLR------------------
eLR.df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds") %>%
  rownames_to_column("LRpairs") %>%
  mutate(Lgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ ii[2]})))


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

setdiff(data.plot$LRpairs,eLR.df$LRpairs)


Lgenelist <- eLR.df$Lgene 
Rgenelist <- eLR.df$Rgene
Lgene <- unique(Lgenelist)
Rgene <- unique(Rgenelist)
LRgene <- c(Lgene,Rgene)
LRgene <- intersect(LRgene,mm9.gene.list)
eLRgene <- LRgene
length(unique(eLRgene))
#export.bed(test,myFileName(prefix = "res/tmp/test_mm9_overlap",suffix = ".bed"))
M <- length(DDD.gene)
N <- length(mm9.gene.list)
k <- length(intersect(LRgene,DDD.gene))
n <- length(LRgene)
x <- M/N
y <- k/n


pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(LRgene.all,DDD.gene)),
                                                  N = length(LRgene.all),
                                                  k = length(intersect(eLRgene,DDD.gene)),
                                                  n = length(eLRgene)))

pvalue.list <- c(pvalue.list,my_enrichment.pvalue(M = length(intersect(mm9.gene.list,DDD.gene)),
                                                  N = length(mm9.gene.list),
                                                  k = length(intersect(eLRgene,DDD.gene)),
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
  ggtitle("DDD gene enrichment")+
  scale_fill_manual(values = c("darkgrey",blues(3)[2:3]))+
  theme_cowplot(font_size = 40)+
  theme(legend.position = "none",
        plot.title = element_text(size = 36,hjust = 0.5),
        axis.line = element_line(size = 2),
        axis.ticks = element_line(size = 2))

ggsave(filename = myFileName(prefix = "res/fig/FigS_DDD_gene_enrichment",
                             suffix = ".jpg"),
       width = 8,height = 8,dpi = 350)

ggsave(filename = myFileName(prefix = "res/fig/FigS_DDD_gene_enrichment",
                             suffix = ".pdf"),
       width = 8,height = 8,dpi = 350)

#####------1.5 euler plot---------
tmp.list <- list(all = DDD.gene,
                 LR = LRgene.all,
                 eLR = eLRgene)
names(tmp.list)[1] <- "DDD gene"
jpeg(filename = myFileName(prefix = "res/fig/FigS_DDD_gene_euler",suffix = ".jpg"),
    width = 8,height = 8,units = "in",bg = "grey",res = 350)
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="DDD gene overlap",fontsize = 24),
     bg="grey")
dev.off()

pdf(file = myFileName(prefix = "res/fig/FigS_DDD_gene_euler",suffix = ".pdf"),
    width = 8,height = 8,bg = "grey")
####beacuse venn digaram involved random
set.seed(666)
plot(eulerr::euler(tmp.list,shape = "ellipse"),
     quantities = list(fontsize=24),
     label=list(fontsize=24),
     fill=c("white",blues(4)[-1]),
     col="black",
     main=list(label="DDD gene overlap",fontsize = 24),
     bg="grey")
dev.off()


### load HK gene
tmp.dir <- "database/HRT_housekeeping_gene_database/Housekeeping_GenesMouse.csv"
tmp.df <- read.delim(file = tmp.dir,sep = ";",stringsAsFactors = F,header = T)
HK.gene <- tmp.df %>%
  pull(Gene) %>%
  unique()
HK.gene <- intersect(HK.gene,mm9.gene.list)



