### TimeTalk supplementary analysis
### Author: wlt
### Date: 20230422
### Finish Date: 20230425

###Co-evolution analysis

###------0.load package---------------
library(rentrez)
library(msa)
library(tidyverse)
library(Biostrings)
library(MSA2dist)
library(rtracklayer)
library(genbankr)
library(orthologr)
library(testthat)
library(ggridges)
library(cowplot)
library(RColorBrewer)
library(future.apply)
library(ggbeeswarm)
library(gghalves)
library(ggsignif)
source("./code/myUtils.R")

myGetCDSFromGBK <- function(rawgbk, verbose = FALSE) {
  bf = proc.time()["elapsed"]
  feats = rawgbk$FEATURES
  sq = rawgbk$ORIGIN
  typs = sapply(feats, function(x) if (length(x) > 0) 
    x$type[1]
    else NA_character_)
  empty = is.na(typs)
  feats = feats[!empty]
  typs = typs[!empty]
  featspl = split(feats, typs)
  srcs = genbankr:::fill_stack_df(featspl$source)
  circ = rep(grepl("circular", rawgbk$LOCUS), times = length(srcs))
  genom = genbankr:::gn_from_vers(rawgbk$VERSION)
  sqinfo = Seqinfo(seqlevels(srcs), width(srcs), circ, genom)
  if (verbose) 
    message(Sys.time(), " Starting creation of gene GRanges")
  gns = genbankr:::make_genegr(featspl$gene, sqinfo)
  if (verbose) 
    message(Sys.time(), " Starting creation of CDS GRanges")
  if (!is.null(featspl$CDS)) 
    cdss = genbankr:::make_cdsgr(featspl$CDS, gns, sqinfo)
  else cdss = GRanges(seqinfo = sqinfo)
  
  cdsseq = rawgbk$ORIGIN[cdss]
  tmp.meta <- elementMetadata(cdss)
  tmp.meta$SOURCE <- rawgbk$SOURCE$source
  names(cdsseq) <- paste0(tmp.meta$gene,"_",tmp.meta$SOURCE)
  names(cdsseq) <- gsub(pattern = " ",replacement = "_",names(cdsseq))
  elementMetadata(cdsseq) <- tmp.meta
  af = proc.time()["elapsed"]
  if (verbose) 
    message(Sys.time(), " - Done creating GenBankRecord object [ ", 
            af - bf, " seconds ]")
  return(cdsseq)
}

myGetOrthologoSeq <- function(tmp.symbol = tmp.symbol,
                              tmp.db.type = "homologene",
                              myret_obj = TRUE,
                              tmp.path = "res/fasta/CDS/LR_gene/",
                              tmp.taxid = c("9031","10116","9606","10090","7955","9913","8364")){
  #tmp.symbol <- "Oprl1"
  cat("search homologene\n")
  r_search <- entrez_search(db = tmp.db.type, term=paste0(tmp.symbol,"[GENE] AND mouse[ORGN]"))
  r_search$ids ### 1522
  
  tmp.links <- entrez_link(dbfrom = tmp.db.type, 
                           id = r_search$ids, 
                           db = "gene")
  
  tmp.id <- tmp.links$links$homologene_gene
  cat("search entrez gene\n")
  #tmp.res$links$gene_nuccore_refseqrna
  upload <- entrez_post(db = "gene",
                        id = tmp.id)
  tmp.res <- entrez_link(dbfrom="gene", 
                         db="nuccore", 
                         cmd="neighbor_history",
                         web_history=upload)
  cat("get summary\n")
  #tmp.res$links$gene_nuccore_refseqrna
  tmp_summary <- entrez_summary(db="nuccore",
                                web_history = tmp.res$web_histories$gene_nuccore_refseqrna)
  
  tmp.df <- lapply(tmp_summary, FUN = function(x){
    ttt <- data.frame(title = x$title,
                      caption = x$caption,
                      length = x$statistics$count[1],
                      createdate = x$createdate,
                      organism = x$organism,
                      taxid = x$taxid,
                      stringsAsFactors = F)
    return(ttt)
  })
  tmp.df <- Reduce(rbind,tmp.df)
  #tmp.taxid <- c("9031","10116","9606","10090","7955","9913","8364")
  tmp.df <- tmp.df %>%
    dplyr::filter(taxid %in% tmp.taxid) %>%
    dplyr::filter(grepl('NM', caption)) %>%
    unique() %>%
    mutate(caption_id = gsub(pattern = "NM_",replacement = "",caption)) %>%
    mutate(caption_id = as.integer(caption_id)) %>%
    group_by(taxid) %>%
    dplyr::filter(length == max(length)) %>%
    #dplyr::filter(createdate == min(createdate)) %>%
    dplyr::filter(caption_id == min(caption_id)) %>%
    ungroup()

  tmp.list <- lapply(tmp.df$caption,FUN = function(id){
    tmp.res <- entrez_fetch(db = "nuccore",
                            id = id,
                            rettype = "gbk")
    tmp.obj <- genbankr::parseGenBank(text = tmp.res)
    tmp.seq <- myGetCDSFromGBK(rawgbk = tmp.obj,verbose = T)
    return(tmp.seq)
  })
  tmp.res <- do.call(c,tmp.list)
  Biostrings::writeXStringSet(tmp.res,filepath = paste0(tmp.path,tmp.symbol,"_cds_n.fasta"))
  if(myret_obj == TRUE){
    return(tmp.res)
  }
}

myGetCodonAlign <- function(tmp.symbol,
                            tmp.path = "res/fasta"){
  
  tmp.obj <- myGetOrthologoSeq(tmp.symbol = tmp.symbol)
  tmp.aln <- msa(cds2aa(tmp.obj),order = "input", 
                 method="ClustalW",
                 output = "clustal")
  tmp.aln.1 <- as(tmp.aln,"AAStringSet")
  tmp.res <- cds2aa(tmp.obj)
  Biostrings::writeXStringSet(tmp.aln.1,
                              filepath = paste0(tmp.path,tmp.symbol,"_cds_p_aln.fasta"))
  Biostrings::writeXStringSet(tmp.obj,
                              filepath = paste0(tmp.path,tmp.symbol,"_cds_n.fasta"))
  
  
  codon_aln.1 <- codon_aln(file_aln = paste0(tmp.path,tmp.symbol,"_cds_p_aln.fasta"),
                           file_nuc = paste0(tmp.path,tmp.symbol,"_cds_n.fasta"),
                           format   = "clustal",
                           tool     = "pal2nal",
                           get_aln  = TRUE)
  
  codon_aln.2 <- codon_aln.1 |> aln2dnastring()
  
  return(codon_aln.2)
}

###------1. Getting sequence------------

### all LR pairs ### these are expressed pairs
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_2022031117.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 

LRpairs.df <- data.frame(LRpairs = paste0(Lgenelist,"-",Rgenelist),
                         Lgene = Lgenelist,
                         Rgene = Rgenelist,
                         stringsAsFactors = F)

tmp.gene.list <- unique(c(Lgenelist,Rgenelist))



tmp.res.df <- as.data.frame(tmp.gene.list)


###-------1.1 batch query LR------------------------
lapply(tmp.gene.list,function(x){
  cat(x,sep = "\n")
  cat(x,file = "res/fasta/CDS/LR_gene/log.txt",sep = "\n",append = T)
  tryCatch(
    expr = {
      try_again(10,myGetOrthologoSeq(x,tmp.path = "res/fasta/CDS/LR_gene/"))
  },
  error = function(e){
    cat(paste0(x,"don't have homology"),file = "res/fasta/CDS/LR_gene/log.txt",sep = "\n",append = T)
    message('skip')
  })
})

length(tmp.gene.use)


###----1.2 download append--------
my_appgene <- setdiff(tmp.gene.list,tmp.gene.use)
lapply(my_appgene,FUN = function(x){
  cat(x,sep = "\n")
  cat(x,file = "res/fasta/CDS/LR_gene/log_append.txt",sep = "\n",append = T)
  tryCatch(
    expr = {
      try_again(10,myGetOrthologoSeq(x,tmp.path = "res/fasta/CDS/LR_gene/"))
    },
    error = function(e){
      cat(paste0(x,"don't have homology"),file = "res/fasta/CDS/LR_gene/log.txt",sep = "\n",append = T)
      message('skip')
    })
})

tmp.append.function <- function(tmp.symbol = NULL,
                                tmp.id = NULL,
                                tmp.path = "res/fasta/CDS/LR_gene/"){
  
  
  upload <- entrez_post(db = "gene",
                        id = tmp.id)
  tmp.res <- entrez_link(dbfrom="gene", 
                         db="nuccore",
                         cmd="neighbor_history",
                         web_history=upload)
  cat("get summary\n")
  #tmp.res$links$gene_nuccore_refseqrna
  tmp_summary <- entrez_summary(db="nuccore",
                                version = "1.0",
                                #rettype = "json",
                                web_history = tmp.res$web_histories$gene_nuccore_refseqrna)
  
  tmp.df <- lapply(tmp_summary, FUN = function(x){
    #x <- tmp_summary$`2462534297`
    ttt <- data.frame(title = x$Title,
                      caption = x$Caption,
                      length = x$Length,
                      createdate = x$CreateDate,
                      #organism = x$organism,
                      taxid = x$TaxId,
                      stringsAsFactors = F)
    return(ttt)
  })
  tmp.df <- Reduce(rbind,tmp.df)
  tmp.taxid <- c("9031","10116","9606","10090","7955","9913","8364")
  tmp.df <- tmp.df %>%
    dplyr::filter(taxid %in% tmp.taxid) %>%
    dplyr::filter(grepl('NM', caption)) %>%
    unique() %>%
    group_by(taxid) %>%
    filter(length == max(length))
  
  tmp.list <- lapply(tmp.df$caption,FUN = function(id){
    tmp.res <- entrez_fetch(db = "nuccore",
                            id = id,
                            rettype = "gbk")
    tmp.obj <- genbankr::parseGenBank(text = tmp.res)
    tmp.seq <- myGetCDSFromGBK(rawgbk = tmp.obj,verbose = T)
    return(tmp.seq)
  })
  tmp.res <- do.call(c,tmp.list)
  Biostrings::writeXStringSet(tmp.res,filepath = paste0(tmp.path,tmp.symbol,"_cds_n.fasta"))
  
}

tmp.symbol <- "Nrcam"
tmp.id <- c("1834","666279")

cat("search homologene\n")
r_search <- entrez_search(db = tmp.db.type, term=paste0(tmp.symbol,"[GENE] AND mouse[ORGN]"))
r_search$ids ### 1522

tmp.links <- entrez_link(dbfrom = tmp.db.type, 
                         id = r_search$ids, 
                         db = "gene")

tmp.id <- tmp.links$links$homologene_gene

tmp.res <- entrez_link(dbfrom = "gene",
                       id = tmp.id,
                       db = "nuccore")
length(tmp.res$links$gene_nuccore_refseqrna)

# for( seq_start in seq(1,length(tmp.res$links$gene_nuccore_refseqrna),50)){
#   seq_start <- 1
#   tmp_summary <- entrez_summary(db="nuccore",
#                                 retmax = 50,
#                                 retstart = seq_start,
#                                 web_history = tmp.res$web_histories$gene_nuccore_refseqrna)
#   
#   cat(recs, file="snail_coi.fasta", append=TRUE)
#   cat(seq_start+49, "sequences downloaded\r")
# }

upload <- entrez_post(db = "gene",
                      id = tmp.id)
tmp.res <- entrez_link(dbfrom="gene", 
                       db="nuccore",
                       cmd="neighbor_history",
                       web_history=upload)
cat("get summary\n")
#tmp.res$links$gene_nuccore_refseqrna
tmp_summary <- entrez_summary(db="nuccore",
                              version = "1.0",
                              #rettype = "json",
                              web_history = tmp.res$web_histories$gene_nuccore_refseqrna)

tmp.df <- lapply(tmp_summary, FUN = function(x){
  #x <- tmp_summary$`2462534297`
  ttt <- data.frame(title = x$Title,
                    caption = x$Caption,
                    length = x$Length,
                    createdate = x$CreateDate,
                    #organism = x$organism,
                    taxid = x$TaxId,
                    stringsAsFactors = F)
  return(ttt)
})
tmp.df <- Reduce(rbind,tmp.df)
#tmp.taxid <- c("9031","10116","9606","10090","7955","9913","8364")
tmp.df <- tmp.df %>%
  dplyr::filter(taxid %in% tmp.taxid) %>%
  dplyr::filter(grepl('NM', caption)) %>%
  unique() %>%
  group_by(taxid) %>%
  filter(length == max(length))

tmp.list <- lapply(tmp.df$caption,FUN = function(id){
  tmp.res <- entrez_fetch(db = "nuccore",
                          id = id,
                          rettype = "gbk")
  tmp.obj <- genbankr::parseGenBank(text = tmp.res)
  tmp.seq <- myGetCDSFromGBK(rawgbk = tmp.obj,verbose = T)
  return(tmp.seq)
})
tmp.res <- do.call(c,tmp.list)
Biostrings::writeXStringSet(tmp.res,filepath = paste0(tmp.path,tmp.symbol,"_cds_n.fasta"))

x <- "Cacna1c"
try_again(10,myGetOrthologoSeq(x,tmp.path = "res/fasta/CDS/LR_gene/"))


### HomologGene data
tmp.path <- "database/homologene.data.txt"
Homology_database <- read.delim(file = tmp.path,stringsAsFactors = F,header = F)

tmp.taxid <- c("9031","10116","9606","10090","7955","9913","8364")

tmp.id.use <- Homology_database %>%
  dplyr::filter(V2 %in% tmp.taxid) %>%
  group_by(V1) %>%
  summarise( n = n()) %>%
  dplyr::filter(n == 7) %>%
  pull(V1)

# length(tmp.id.use)
# [1] 7859

###------1.3. read results and append sequence--------

tmp.path <- "res/fasta/CDS/LR_gene/"

tmp.file.use <- list.files(path = tmp.path,
                           pattern = "_cds_n.fasta")

tmp.gene.use <- gsub(pattern = "_cds_n.fasta",
                     replacement = "",
                     tmp.file.use)

setdiff(tmp.gene.list,tmp.gene.use)

myspecies_use <- c("Homo_sapiens_(human)",
                   "Rattus_norvegicus_(Norway_rat)",
                   "Mus_musculus_(house_mouse)",
                   "Bos_taurus_(cattle)",
                   "Danio_rerio_(zebrafish)",
                   "Xenopus_tropicalis_(tropical_clawed_frog)")

myspecies_use <- c("Homo_sapiens_(human)","Mus_musculus_(house_mouse)")

tmp.res.list <- lapply(tmp.file.use,FUN = function(x){
  #x <- grep(pattern = "Fgf4",x = tmp.file.use,value = T)
  cat(x,sep = "\n")
  tmp.path <- "res/fasta/CDS/LR_gene/"
  tmp.res <- readDNAStringSet(filepath = paste0(tmp.path,x))
  tmp.species <- str_split(names(tmp.res),pattern = "_",n = 2)
  tmp.species <- unlist(lapply(tmp.species,
                               FUN = function(x) x[2]))
  if(all(myspecies_use %in% tmp.species)){
    return(tmp.res)
  }
})

# x <- grep(pattern = "Cd44",x = tmp.file.use,value = T)
# 
# cat(x,sep = "\n")
# tmp.path <- "res/fasta/CDS/LR_gene/"
# tmp.res <- readDNAStringSet(filepath = paste0(tmp.path,x))
# tmp.species <- str_split(names(tmp.res),pattern = "_",n = 2)
# tmp.species <- unlist(lapply(tmp.species,
#                              FUN = function(x) x[2]))

x <- tmp.res.list
names(x) <- tmp.gene.use
x[sapply(x,is.null)] <- NULL
grep("Tgf",names(x),value = T)

tmp.homology.gene.list <- names(x)
tmp.query <- setdiff(tmp.gene.list,names(x))



tmp.symbol <- "Bmp8b"

# x <- "Bmp8b"
# try_again(10,myGetOrthologoSeq(x,tmp.path = "res/fasta/CDS/LR_gene/"))
entrez_db_searchable(db = "gene")
.temp_append <- function(tmp.symbol){
  cat(tmp.symbol,sep = "\n")
  #tmp.symbol <- "Sct"
  tmp.query <- paste0(tmp.symbol,"[PREF]")
  #entrez_db_searchable("gene")
  tmp.gene <- entrez_search(db = "gene",term = tmp.query,
                            retmax = 30)
  tmp.append.function(tmp.symbol = tmp.symbol,
                      tmp.path = "res/fasta/CDS/LR_gene/",
                      tmp.id = tmp.gene$ids)
}

temp_gene <- setdiff(tmp.gene.list,names(x))
.temp_append(tmp.symbol = "Clec2d")
lapply(temp_gene,FUN = function(x){
  try_again(times = 10,.temp_append(x))
})
tmp.query <- paste0(tmp.symbol,"[GENE]")
#entrez_db_searchable("gene")
tmp.gene <- entrez_search(db = "gene",term = tmp.query,
                          retmax = 30)
tmp.append.function(tmp.symbol = tmp.symbol,
                    tmp.path = "res/fasta/CDS/LR_gene/",
                    tmp.id = tmp.gene$ids)

LRpairs.have5homolog <- LRpairs.df %>%
  dplyr::filter((Lgene %in% names(x)) & (Rgene %in% names(x)))


#### load eLR data
eLRpairs.df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds") %>%
  rownames_to_column("LRpairs") %>%
  mutate(Lgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ii[2]})))

tmp.df <- LRpairs.have5homolog %>%
  dplyr::filter(LRpairs %in% eLRpairs.df$LRpairs)

setdiff(eLRpairs.df$LRpairs,tmp.df$LRpairs)


####-------2. read and perform analysis------------


tmp.path <- "res/fasta/CDS/LR_gene/"
tmp.file.use <- list.files(path = tmp.path,
                           pattern = "_cds_n.fasta")
tmp.gene.use <- gsub(pattern = "_cds_n.fasta",
                     replacement = "",
                     tmp.file.use)
length(tmp.LRgenelist)
length(tmp.gene.use)
setdiff(tmp.gene.list,tmp.gene.use)

# myspecies_use <- c("Homo_sapiens_(human)",
#                    "Rattus_norvegicus_(Norway_rat)",
#                    "Mus_musculus_(house_mouse)",
#                    "Bos_taurus_(cattle)",
#                    "Danio_rerio_(zebrafish)",
#                    "Xenopus_tropicalis_(tropical_clawed_frog)")
myspecies_use <- c("Homo_sapiens_(human)",
                   "Rattus_norvegicus_(Norway_rat)",
                   "Mus_musculus_(house_mouse)",
                   "Bos_taurus_(cattle)")

#myspecies_use <- c("Homo_sapiens_(human)","Mus_musculus_(house_mouse)")

tmp.res.list <- lapply(tmp.file.use,FUN = function(x){
  #x <- grep(pattern = "Fgf4",x = tmp.file.use,value = T)
  cat(x,sep = "\n")
  tmp.path <- "res/fasta/CDS/LR_gene/"
  tmp.res <- readDNAStringSet(filepath = paste0(tmp.path,x))
  tmp.species <- str_split(names(tmp.res),pattern = "_",n = 2)
  tmp.species <- unlist(lapply(tmp.species,
                               FUN = function(x) x[2]))
  if(all(myspecies_use %in% tmp.species)){
    return(tmp.res)
  }
})

# x <- grep(pattern = "Cd44",x = tmp.file.use,value = T)
# 
# cat(x,sep = "\n")
# tmp.path <- "res/fasta/CDS/LR_gene/"
# tmp.res <- readDNAStringSet(filepath = paste0(tmp.path,x))
# tmp.species <- str_split(names(tmp.res),pattern = "_",n = 2)
# tmp.species <- unlist(lapply(tmp.species,
#                              FUN = function(x) x[2]))

x <- tmp.res.list
names(x) <- tmp.gene.use
x[sapply(x,is.null)] <- NULL
grep("Tgf",names(x),value = T)

tmp.query <- setdiff(tmp.gene.list,names(x))


#### load eLR data
eLRpairs.df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds") %>%
  rownames_to_column("LRpairs") %>%
  mutate(Lgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ii[2]})))
LRpairs.have5homolog <- LRpairs.df %>%
  dplyr::filter((Lgene %in% names(x)) & (Rgene %in% names(x)))


tmp.df.eLR <- LRpairs.have5homolog %>%
  dplyr::filter(LRpairs %in% eLRpairs.df$LRpairs)

tmp.df.non.eLR <- LRpairs.have5homolog %>%
  dplyr::filter(!(LRpairs %in% eLRpairs.df$LRpairs))


tmp.idx <- 1
tmp.path <- "res/fasta/CDS/LR_gene/"
x <- tmp.df.eLR$Lgene[tmp.idx]
y <- tmp.df.eLR$Rgene[tmp.idx]

x
y
.temp_cds_aln <- function(tmp.symbol,
                          myspecies_use =  c("Homo_sapiens_(human)",
                                             "Rattus_norvegicus_(Norway_rat)",
                                             "Mus_musculus_(house_mouse)",
                                             "Bos_taurus_(cattle)")){
  #tmp.symbol <- "Oprl1"
  tmp.path <- "res/fasta/CDS/LR_gene/"
  tmp.file.path <- paste0(tmp.path,tmp.symbol,
                          "_cds_n.fasta")
  tmp.obj <- readDNAStringSet(filepath = tmp.file.path)
  tmp.name.list <- str_split(names(tmp.obj),
                             pattern = "_",n=2)
  tmp.name.df <- lapply(tmp.name.list,FUN = function(idx){
    res.df <- data.frame(gene_symbol = idx[1],
                         organism = idx[2],stringsAsFactors = F)
    return(res.df)
  })
  tmp.name.df <- Reduce(rbind,tmp.name.df)
  
  
  tmp.idx <- tmp.name.df %>%
    mutate(myseqname = paste0(gene_symbol,'_',organism)) %>%
    dplyr::filter(organism %in% myspecies_use) %>%
    dplyr::arrange(match(organism,myspecies_use)) %>%
    pull(myseqname)
  
  tmp.obj <- tmp.obj[tmp.idx]
  tmp.obj
  tmp.aln <- msa(cds2aa(tmp.obj),
                 order = "input", 
                 method="ClustalW",
                 output = "clustal")
  tmp.aln.1 <- as(tmp.aln,"AAStringSet")
  #tmp.res <- cds2aa(tmp.obj)
  
  tmp.path <- "res/fasta/CDS/codon_aln/"
  Biostrings::writeXStringSet(tmp.aln.1,
                              filepath = paste0(tmp.path,tmp.symbol,"_cds_p_aln.fasta"))
  Biostrings::writeXStringSet(tmp.obj,
                              filepath = paste0(tmp.path,tmp.symbol,"_cds_n.fasta"))
  
  codon_aln.1 <- codon_aln(file_aln = paste0(tmp.path,tmp.symbol,"_cds_p_aln.fasta"),
                           file_nuc = paste0(tmp.path,tmp.symbol,"_cds_n.fasta"),
                           format   = "clustal",
                           tool     = "pal2nal",
                          get_aln  = TRUE)
  
  codon_aln.2 <- codon_aln.1 |> aln2dnastring() 
  
  return(codon_aln.2)
}

tmp.pair.df.use <- LRpairs.have5homolog
tmp.coevole_metric.list <- lapply(1:nrow(tmp.pair.df.use),
                                FUN = function(tmp.idx){
                                  tryCatch(expr = {
                                    cat(tmp.pair.df.use$LRpairs[tmp.idx],sep = "\n")
                                    
                                    # x <- "Calm3"
                                    # y <- "Grm7"
                                    x <- tmp.pair.df.use$Lgene[tmp.idx]
                                          
                                    y <- tmp.pair.df.use$Rgene[tmp.idx]
                                           
                                           
                                    x.aln <- .temp_cds_aln(tmp.symbol = x,myspecies_use = myspecies_use)
                                           
                                    #x.res <- x.aln |> cds2aa() |> aastring2dist(score = granthamMatrix())
                                    
                                    x.res <- x.aln |> dnastring2dist(model = "K80")
                                    
                                    x.res <- x.res$distSTRING
                                           
                                    y.aln <- .temp_cds_aln(tmp.symbol = y,myspecies_use = myspecies_use)
                                           
                                    #y.res <- y.aln |> cds2aa() |> aastring2dist(score = granthamMatrix())
                                    
                                    y.res <- y.aln |> dnastring2dist(model = "K80")       
                                    y.res <- y.res$distSTRING
                                           
                                    tmp.res <- allCorrelations(x.res,y.res,ncomp1 = 3,ncomp2 = 3)
                                           
                                    tmp.res.df <- data.frame(t(tmp.res))
                                           
                                    tmp.res.df$LR <- paste0(x,"-",y)
                                           
                                    return(tmp.res.df)
                                           
                                    },
                                    
                                    error = function(e){
                                      message("skip")       
                                           
                                      })
                                  
            
                                })
                              
tmp.coevole_metric.df <- Reduce(rbind,tmp.coevole_metric.list)


tmp.res.df <- tmp.coevole_metric.df %>%
  dplyr::mutate(group = ifelse(LR %in% eLRpairs.df$LRpairs,
                               "eLR","not_eLR"))

x <- tmp.res.df %>%
  dplyr::filter(group == "eLR") %>%
  pull(RV2)
y <- tmp.res.df %>%
  dplyr::filter(group == "not_eLR") %>%
  pull(RV2)

tmp.p.res <- wilcox.test(x,y,
            alternative = "two.sided")


p <- ggplot(tmp.res.df,aes(x = RV2,y = group,fill=group))+
  geom_density_ridges2()+
  scale_x_continuous(limits = c(0.8,1),expand = c(0,0))+
  theme_cowplot(font_size = 28)+
  scale_fill_manual(values = c("#8FD3F7","#E8E0D1"))+
  ylab(NULL)+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  ggtitle(label = paste0("p=",signif(tmp.p.res$p.value,4)))
p  
ggsave(filename = paste0(myFileName(prefix = "res/fig/FigS_co_evovle_RV2_ridge_plot",suffix = ".jpg")),
       width = 8,height = 6,dpi = 350)  
ggsave(filename = paste0(myFileName(prefix = "res/fig/FigS_co_evovle_RV2_ridge_plot",suffix = ".pdf")),
       width = 8,height = 6)  


tmp.pair.df.use <- LRpairs.have5homolog
head(tmp.pair.df.use)
nrow(tmp.pair.df.use)
ncol(tmp.pair.df.use)

tmp.idx <- 1

tmp.coevole_kaks.list <- lapply(1:nrow(tmp.pair.df.use),
                                  FUN = function(tmp.idx){
                                    tryCatch(expr = {
                                      
                                      cat(tmp.pair.df.use$LRpairs[tmp.idx],sep = "\n")
                                      
                                      # x <- "Calm3"
                                      # y <- "Grm7"
                                      x <- tmp.pair.df.use$Lgene[tmp.idx]
                                      
                                      y <- tmp.pair.df.use$Rgene[tmp.idx]
                                      
                                      x.aln <- .temp_cds_aln(tmp.symbol = x)
                                      
                                      #x.res <- x.aln |> cds2aa() |> aastring2dist(score = granthamMatrix())
                                      
                                      x.res <- x.aln |> dnastring2kaks() %>%
                                        mutate(kaks_ratio = ka/ks) %>%
                                        pull(kaks_ratio) 

                                      y.aln <- .temp_cds_aln(tmp.symbol = y)
                                      
                                      #y.res <- y.aln |> cds2aa() |> aastring2dist(score = granthamMatrix())
                                      
                                      y.res <- y.aln |> dnastring2kaks() %>%
                                        mutate(kaks_ratio = ka/ks) %>%
                                        pull(kaks_ratio) 
                                      
                                      tmp.p.value <- cor.test(x.res,y.res)
                                      tmp.pcc <- cor(x.res,y.res)
                                      tmp.pcc[is.na(tmp.pcc)] <- 0
                                      
                                      tmp.res.df <- data.frame(PCC = tmp.pcc,
                                                               PCC.tes = tmp.p.value$p.value,
                                                               stringsAsFactors = F)
                                      
                                      tmp.res.df$LR <- paste0(x,"-",y)
                                      
                                      return(tmp.res.df)
                                      
                                    },
                                    
                                    error = function(e){
                                      message("skip")       
                                      
                                    })
                                    
                                    
                                  })

tmp.coevole_kaks.df <- Reduce(rbind,tmp.coevole_kaks.list)
head(tmp.coevole_kaks.df)

plot(tmp.coevole_kaks.df$PCC)

tmp.coevole_kaks.df <- tmp.coevole_kaks.df %>%
  dplyr::mutate(group = ifelse(LR %in% eLRpairs.df$LRpairs,
                               "eLR","not_eLR")) %>%
  mutate(Lgene=unlist(lapply(str_split(LR,pattern = "-"),FUN = function(ii){ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LR,pattern = "-"),FUN = function(ii){ii[2]}))) %>%
  arrange(-PCC)




###--------2.1 plot distribution-----------------

tmp.data.plot <- tmp.coevole_kaks.df 

x <- tmp.data.plot %>%
  filter(group == "eLR") %>%
  pull(PCC)

y <- tmp.data.plot %>%
  filter(group == "not_eLR") %>%
  pull(PCC)

tmp.res.p <- wilcox.test(x,y)
tmp.res.p$p.value
head(tmp.data.plot)


ggplot(tmp.data.plot,aes(PCC,group,fill = group))+
  geom_density_ridges2()+
  scale_x_continuous(limits = c(-1,1))+
  theme_cowplot(font_size = 28)+
  scale_fill_manual(values = c("#8FD3F7","#E8E0D1"))+
  ylab(NULL)+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")+
  ggtitle(label = paste0("p=",signif(tmp.p.res$p.value,4)))
ggsave(filename = paste0(myFileName(prefix = "res/fig/FigS_co_evovle_PCC_ridge_plot",suffix = ".jpg")),
       width = 8,height = 6,dpi = 350)  
ggsave(filename = paste0(myFileName(prefix = "res/fig/FigS_co_evovle_PCC_ridge_plot",suffix = ".pdf")),
       width = 8,height = 6)  
  

###---------2.2 EDA------------------

PCC.data.plot <- readRDS(file = "res/R/putative_eLR_pairs_1094_2022032109.rds") %>%
  mutate(LR = paste0(Lgene,"-",Rgene),PCC_covary = PCC) %>%
  select(LR,PCC_covary)

tmp.pair.df.use <- tmp.coevole_kaks.df %>%
  mutate(Lgene=unlist(lapply(str_split(LR,pattern = "-"),FUN = function(ii){ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LR,pattern = "-"),FUN = function(ii){ii[2]}))) %>%
  arrange(-PCC)

length(tmp.pair.df.use$LR)
length(PCC.data.plot$LR)
setdiff(PCC.data.plot$LR,tmp.pair.df.use$LR)

setdiff(tmp.pair.df.use$LR,PCC.data.plot$LR)
intersect(tmp.pair.df.use$LR,PCC.data.plot$LR)


tmp.pair.df.use <- merge(tmp.pair.df.use,PCC.data.plot,by = "LR")

head(tmp.pair.df.use)


head(tmp.pair.df.use)

tmp.idx <- 1
.temp_coevovle_plot <- function(x,y){
  x.aln <- .temp_cds_aln(tmp.symbol = x)
  #x.res <- x.aln |> cds2aa() |> aastring2dist(score = granthamMatrix())
  x.res <- x.aln |> dnastring2kaks() %>%
    mutate(kaks_ratio = ka/ks) %>%
    pull(kaks_ratio) 
  
  y.aln <- .temp_cds_aln(tmp.symbol = y)
  #y.res <- y.aln |> cds2aa() |> aastring2dist(score = granthamMatrix())
  y.res <- y.aln |> dnastring2kaks() %>%
    mutate(kaks_ratio = ka/ks) %>%
    pull(kaks_ratio) 
  
  tmp.pcc <- cor(x.res,y.res)
  tmp.pcc[is.na(tmp.pcc)] <- 0
  
  tmp.data.plot <- data.frame(kaks_L = x.res,
                              kaks_R = y.res,
                              PCC = cor(x.res,y.res),
                              LRpairs = paste0(x,"-",y),
                              stringsAsFactors = F)
  
  p <- ggplot(tmp.data.plot,aes(kaks_L,kaks_R))+
    geom_point(size = 3)+
    geom_smooth(method = "lm",se = F,color = "blue")+
    ggtitle(label = paste0(tmp.data.plot$LRpairs,"\n",
                           "PCC=",
                           signif(tmp.data.plot$PCC,digits = 5)))+
    xlab("Ligand Ka/Ks")+
    ylab("Receptor Ka/Ks")+
    theme_cowplot(font_size = 28)+
    theme(plot.title = element_text(hjust = 0.5))
  p
  ggsave(filename = myFileName(prefix = paste0("res/fig/",x,"-",y,"_coevolve_","KaKs_ratio"),
                               suffix = ".jpg"),
         width = 8,height = 6,dpi = 350)
  ggsave(filename = myFileName(prefix = paste0("res/fig/",x,"-",y,"_coevolve_","KaKs_ratio"),
                               suffix = ".pdf"),
         width = 8,height = 6,dpi = 350)
  return(p)
}

head(tmp.pair.df.use)

.temp_coevovle_plot(x = "Ccl21b",
                    y = "Cxcr3")

.temp_coevovle_plot(x = "Ccl21b",
                    y = "Cxcr3")

.temp_coevovle_plot(x = "Fgf10",
                    y = "Fgfr2")

#Fgf8-Fgfr1

plot(tmp.pair.df.use$PCC)

.temp_coevovle_plot(x = "Fgf8",
                    y = "Fgfr1")




tmp.pair.df.use <- tmp.coevole_kaks.df %>%
  dplyr::mutate(group = ifelse(LR %in% eLRpairs.df$LRpairs,
                               "eLR","not_eLR"))

myspecies_use


###-------3. gene overlap and hyper geometric test------------

###-------3.1 read data---------------

tmp.path <- "res/fasta/CDS/LR_gene/"
tmp.file.use <- list.files(path = tmp.path,
                           pattern = "_cds_n.fasta")

tmp.gene.use <- gsub(pattern = "_cds_n.fasta",
                     replacement = "",
                     tmp.file.use)
# myspecies_use <- c("Homo_sapiens_(human)",
#                    "Rattus_norvegicus_(Norway_rat)",
#                    "Mus_musculus_(house_mouse)",
#                    "Bos_taurus_(cattle)",
#                    "Danio_rerio_(zebrafish)",
#                    "Xenopus_tropicalis_(tropical_clawed_frog)")
myspecies_use <- c("Homo_sapiens_(human)",
                   "Rattus_norvegicus_(Norway_rat)",
                   "Mus_musculus_(house_mouse)",
                   "Bos_taurus_(cattle)")

#myspecies_use <- c("Homo_sapiens_(human)","Mus_musculus_(house_mouse)")

tmp.res.list <- lapply(tmp.file.use,FUN = function(x){
  #x <- grep(pattern = "Fgf4",x = tmp.file.use,value = T)
  cat(x,sep = "\n")
  tmp.path <- "res/fasta/CDS/LR_gene/"
  tmp.res <- readDNAStringSet(filepath = paste0(tmp.path,x))
  tmp.species <- str_split(names(tmp.res),pattern = "_",n = 2)
  tmp.species <- unlist(lapply(tmp.species,
                               FUN = function(x) x[2]))
  if(all(myspecies_use %in% tmp.species)){
    return(tmp.res)
  }
})

x <- tmp.res.list
names(x) <- tmp.gene.use
x[sapply(x,is.null)] <- NULL

homology_gene_use <- names(x)

####--------3.2 load LR dataframe -----------------------

#### load LR data.frame
LRpairs <- readRDS(file = "res/R/cytotalkdb_LRpairs_2022031117.rds")
Lgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[1])) 
Rgenelist <- unlist(lapply(strsplit(LRpairs,split = "_"),FUN = function(x) x[2])) 
LRpairs.df <- data.frame(LRpairs = paste0(Lgenelist,"-",Rgenelist),
                         Lgene = Lgenelist,
                         Rgene = Rgenelist,
                         stringsAsFactors = F)
#### load eLR data.frame
eLRpairs.df <- readRDS(file = "res/R/eLR_row_annotaion_df_2022032700.rds") %>%
  rownames_to_column("LRpairs") %>%
  mutate(Lgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LRpairs,pattern = "-"),FUN = function(ii){ii[2]})))

### 922 eLR composed by homology LRpairs data.frame
LRpairs.have5homolog <- LRpairs.df %>%
  dplyr::filter((Lgene %in% homology_gene_use) & (Rgene %in% homology_gene_use))



####--------3.3 calculate the KaKs metric---------------

tmp.idx <- 1
tmp.path <- "res/fasta/CDS/LR_gene/"

.temp_cds_aln <- function(tmp.symbol,
                          myspecies_use =  c("Homo_sapiens_(human)",
                                             "Rattus_norvegicus_(Norway_rat)",
                                             "Mus_musculus_(house_mouse)",
                                             "Bos_taurus_(cattle)")){
  #tmp.symbol <- "Oprl1"
  tmp.path <- "res/fasta/CDS/LR_gene/"
  tmp.file.path <- paste0(tmp.path,tmp.symbol,
                          "_cds_n.fasta")
  tmp.obj <- readDNAStringSet(filepath = tmp.file.path)
  tmp.name.list <- str_split(names(tmp.obj),
                             pattern = "_",n=2)
  tmp.name.df <- lapply(tmp.name.list,FUN = function(idx){
    res.df <- data.frame(gene_symbol = idx[1],
                         organism = idx[2],stringsAsFactors = F)
    return(res.df)
  })
  tmp.name.df <- Reduce(rbind,tmp.name.df)
  
  
  tmp.idx <- tmp.name.df %>%
    mutate(myseqname = paste0(gene_symbol,'_',organism)) %>%
    dplyr::filter(organism %in% myspecies_use) %>%
    dplyr::arrange(match(organism,myspecies_use)) %>%
    pull(myseqname)
  
  tmp.obj <- tmp.obj[tmp.idx]
  tmp.obj
  tmp.aln <- msa(cds2aa(tmp.obj),
                 order = "input", 
                 method="ClustalW",
                 output = "clustal")
  tmp.aln.1 <- as(tmp.aln,"AAStringSet")
  #tmp.res <- cds2aa(tmp.obj)
  
  tmp.path <- "res/fasta/CDS/codon_aln/"
  Biostrings::writeXStringSet(tmp.aln.1,
                              filepath = paste0(tmp.path,tmp.symbol,"_cds_p_aln.fasta"))
  Biostrings::writeXStringSet(tmp.obj,
                              filepath = paste0(tmp.path,tmp.symbol,"_cds_n.fasta"))
  
  codon_aln.1 <- codon_aln(file_aln = paste0(tmp.path,tmp.symbol,"_cds_p_aln.fasta"),
                           file_nuc = paste0(tmp.path,tmp.symbol,"_cds_n.fasta"),
                           format   = "clustal",
                           tool     = "pal2nal",
                           get_aln  = TRUE)
  
  codon_aln.2 <- codon_aln.1 |> aln2dnastring() 
  
  return(codon_aln.2)
}


tmp.pair.df.use <- LRpairs.have5homolog


tmp.coevole_kaks.list <- lapply(1:nrow(tmp.pair.df.use),
                                FUN = function(tmp.idx){
                                  tryCatch(expr = {
                                    
                                    cat(tmp.pair.df.use$LRpairs[tmp.idx],sep = "\n")
                                    
                                    # x <- "Calm3"
                                    # y <- "Grm7"
                                    x <- tmp.pair.df.use$Lgene[tmp.idx]
                                    
                                    y <- tmp.pair.df.use$Rgene[tmp.idx]
                                    
                                    x.aln <- .temp_cds_aln(tmp.symbol = x)
                                    
                                    #x.res <- x.aln |> cds2aa() |> aastring2dist(score = granthamMatrix())
                                    
                                    x.res <- x.aln |> dnastring2kaks() %>%
                                      mutate(kaks_ratio = ka/ks) %>%
                                      pull(kaks_ratio) 
                                    
                                    y.aln <- .temp_cds_aln(tmp.symbol = y)
                                    
                                    #y.res <- y.aln |> cds2aa() |> aastring2dist(score = granthamMatrix())
                                    
                                    y.res <- y.aln |> dnastring2kaks() %>%
                                      mutate(kaks_ratio = ka/ks) %>%
                                      pull(kaks_ratio) 
                                    
                                    tmp.p.value <- cor.test(x.res,y.res)
                                    tmp.pcc <- cor(x.res,y.res)
                                    tmp.pcc[is.na(tmp.pcc)] <- 0
                                    
                                    tmp.res.df <- data.frame(PCC = tmp.pcc,
                                                             PCC.tes = tmp.p.value$p.value,
                                                             stringsAsFactors = F)
                                    
                                    tmp.res.df$LR <- paste0(x,"-",y)
                                    
                                    return(tmp.res.df)
                                    
                                  },
                                  
                                  error = function(e){
                                    message("skip")       
                                    
                                  })
                                  
                                  
                                })

tmp.coevole_kaks.df <- Reduce(rbind,tmp.coevole_kaks.list)
head(tmp.coevole_kaks.df)

plot(tmp.coevole_kaks.df$PCC)

tmp.coevole_kaks.df <- tmp.coevole_kaks.df %>%
  dplyr::mutate(group = ifelse(LR %in% eLRpairs.df$LRpairs,
                               "eLR","not_eLR")) %>%
  mutate(Lgene=unlist(lapply(str_split(LR,pattern = "-"),FUN = function(ii){ii[1]}))) %>%
  mutate(Rgene=unlist(lapply(str_split(LR,pattern = "-"),FUN = function(ii){ii[2]}))) %>%
  arrange(-PCC)

plot(tmp.coevole_kaks.df$PCC)

saveRDS(tmp.coevole_kaks.df,file = "res/R/LR_co-evolve_20230430.rds")

####-----3.4 visulization-----------------

.temp_coevovle_plot <- function(x,y){
  x.aln <- .temp_cds_aln(tmp.symbol = x)
  #x.res <- x.aln |> cds2aa() |> aastring2dist(score = granthamMatrix())
  x.res <- x.aln |> dnastring2kaks() %>%
    mutate(kaks_ratio = ka/ks) %>%
    pull(kaks_ratio) 
  
  y.aln <- .temp_cds_aln(tmp.symbol = y)
  #y.res <- y.aln |> cds2aa() |> aastring2dist(score = granthamMatrix())
  y.res <- y.aln |> dnastring2kaks() %>%
    mutate(kaks_ratio = ka/ks) %>%
    pull(kaks_ratio) 
  
  tmp.pcc <- cor(x.res,y.res)
  tmp.pcc[is.na(tmp.pcc)] <- 0
  
  tmp.data.plot <- data.frame(kaks_L = x.res,
                              kaks_R = y.res,
                              PCC = cor(x.res,y.res),
                              LRpairs = paste0(x,"-",y),
                              stringsAsFactors = F)
  
  p <- ggplot(tmp.data.plot,aes(kaks_L,kaks_R))+
    geom_point(size = 3)+
    geom_smooth(method = "lm",se = F,color = "blue")+
    ggtitle(label = paste0(tmp.data.plot$LRpairs,"\n",
                           "PCC=",
                           signif(tmp.data.plot$PCC,digits = 5)))+
    xlab("Ligand Ka/Ks")+
    ylab("Receptor Ka/Ks")+
    theme_cowplot(font_size = 28)+
    theme(axis.line = element_line(size = 1),
          axis.ticks = element_line(size = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none")
  p
  ggsave(filename = myFileName(prefix = paste0("res/fig/",x,"-",y,"_coevolve_","KaKs_ratio"),
                               suffix = ".jpg"),
         width = 8,height = 6,dpi = 350)
  ggsave(filename = myFileName(prefix = paste0("res/fig/",x,"-",y,"_coevolve_","KaKs_ratio"),
                               suffix = ".pdf"),
         width = 8,height = 6,dpi = 350)
  return(p)
}

####Ensure LR pairs is accurate
####setdiff(tmp.coevole_kaks.df$LR,LRpairs.df$LRpairs)


####Plot examples
.temp_coevovle_plot(x = "Ccl21b",
                    y = "Cxcr3")

.temp_coevovle_plot(x = "Ccl21b",
                    y = "Cxcr3")

.temp_coevovle_plot(x = "Fgf10",
                    y = "Fgfr2")

saveRDS(tmp.coevole_kaks.df,
        file = myFileName(prefix = "res/R/Homologogy_LR_pairs_coevolve_KaKs",
                          suffix = ".rds"))

### 922 LR pairs
tmp.coevole_kaks.df <- readRDS(file = "res/R/Homologogy_LR_pairs_coevolve_KaKs_2023050110.rds")
table(tmp.coevole_kaks.df$group)

####Plot overall plot
tmp.data.plot <- tmp.coevole_kaks.df 
x <- tmp.data.plot %>%
  filter(group == "eLR") %>%
  pull(PCC)
y <- tmp.data.plot %>%
  filter(group == "not_eLR") %>%
  pull(PCC)

tmp.res.p <- wilcox.test(x,y)
tmp.res.p$p.value
head(tmp.data.plot)


ggplot(tmp.data.plot,
       aes(group,PCC,fill = group))+
  geom_boxplot()
 

ggplot(tmp.data.plot,aes(group,PCC,fill = group))+
  geom_boxplot()+
  geom_signif(comparisons = list(c("eLR","not_eLR")),
              annotations = paste0("p=",signif(tmp.res.p$p.value,4)),
              y_position = 1.2,
              size = 1,
              tip_length = 0.05,
              textsize = 10)+
  scale_y_continuous(limits = c(-1,1.5))+
  theme_cowplot(font_size = 28)+
  scale_fill_manual(values = c("#8FD3F7","#E8E0D1"))+
  xlab(NULL)+
  ylab("KaKs PCC")+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none")
ggsave(filename = paste0(myFileName(prefix = "res/fig/FigS_co_evovle_PCC_ridge_plot",suffix = ".jpg")),
       width = 8,height = 6,dpi = 350)  
ggsave(filename = paste0(myFileName(prefix = "res/fig/FigS_co_evovle_PCC_ridge_plot",suffix = ".pdf")),
       width = 8,height = 6)  







