
myCellChatToSeu <- function(object,slot.use="data"){
  meta <- object@meta
  if (is.list(object@idents)) {
    meta$group.cellchat <- object@idents$joint
  }
  else {
    meta$group.cellchat <- object@idents
  }
  w10x <- Seurat::CreateSeuratObject(counts = slot(object = object,name = slot.use),
                                     meta.data = meta)
  Seurat::Idents(w10x) <- "group.cellchat"
  return(w10x)
}

myscAreaPlot <- function(tmp.seu = tmp.seu, 
                         features = "Col1a1",
                         tmp.color= tmp.color){
  data <- Seurat::FetchData(object = tmp.seu, vars = features, slot = "data")
  data$label <- Idents(tmp.seu)
  data <- data %>%
    gather(key = "features",value = "gene.exp",-label)
  
  data <- data %>%
    group_by(features,label) %>%
    mutate(feature = (gene.exp-min(gene.exp))/(max(gene.exp)-min(gene.exp))) %>%
    arrange(-feature) %>%
    mutate(index = row_number()) %>%
    mutate(index = (index - min(index))/(max(index)-min(index))) %>%
    ungroup()
  
  data$features <- factor(data$features,levels = features)
  
  p <- ggplot(data = data,aes(x=-index,y= feature ,fill=label))+
    geom_area(color=NA)+
    facet_grid(features~label)+
    scale_fill_manual(values = tmp.color)+
    xlab(NULL)+
    ylab("Relative gene expression level per cluster")+
    scale_y_continuous(expand = c(0,0))+
    theme_few()+
    theme(axis.line  = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          strip.background.x = element_rect(color = "black",fill="white"),
          strip.background.y = element_rect(color = "black",fill="white"),
          strip.text.y = element_text(angle = 0,face = "italic"),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          panel.spacing = unit(0,"lines"),
          legend.position = "none")
  p
  return(p)
}



myPlotCellChatFeature <- function (object, features = NULL, signaling = NULL, enriched.only = TRUE, 
                                   using.slot="data",
                                   type = c("violin", "dot","scArea"), color.use = NULL, 
                                   group.by = NULL, ...) 
{
  type <- match.arg(type)
  meta <- object@meta
  
  if (is.list(object@idents)) {
    meta$group.cellchat <- object@idents$joint
  }
  else {
    meta$group.cellchat <- object@idents
  }
  w10x <- Seurat::CreateSeuratObject(counts = slot(object = object,name = using.slot), 
                                     meta.data = meta)
  if (is.null(group.by)) {
    group.by <- "group.cellchat"
  }
  Seurat::Idents(w10x) <- group.by
  
  if (is.null(color.use)) {
    numCluster <- length(levels(Seurat::Idents(w10x)))
    color.use <- scPalette(numCluster)
  }
  if (!is.null(features) & !is.null(signaling)) {
    warning("`features` will be used when inputing both `features` and `signaling`!")
  }
  if (!is.null(features)) {
    feature.use <- features
  }
  else if (!is.null(signaling)) {
    res <- extractEnrichedLR(object, signaling = signaling, 
                             geneLR.return = TRUE, enriched.only = enriched.only)
    feature.use <- res$geneLR
  }
  if (type == "violin") {
    gg <- StackedVlnPlot(w10x, features = feature.use, color.use = color.use, 
                         ...)
  }
  if (type == "dot") {
    gg <- dotPlot(w10x, features = feature.use, ...)
  }
  if (type == "scArea"){
    gg <- myscAreaPlot(w10x,features = feature.use,tmp.color = color.use)
  }
  return(gg)
}



####show color bar

my_showcolorbar <- function(tmp.colorbar,tmp.name,n=100){
  image(
    1:n, 1, as.matrix(1:n),
    col = tmp.colorbar(n),
    xlab = tmp.name, 
    ylab = "", 
    xaxt = "n", 
    yaxt = "n", 
    bty = "n"
  )
}


####Define function


my_getLR <- function(tmp.L.exp = data.plot[Lgene.list,tmp.id],
         tmp.R.exp = data.plot[Rgene.list,tmp.id],
         tmp.p.cutoff=0.05,tmp.cell.cutoff=3){
  
  cell_number <- ncol(tmp.L.exp)
  
  tmp.res.PCC <- sapply(1:nrow(tmp.L.exp), function(i) {
    cor(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),method = "pearson")
  })
  tmp.res.PCC.p <- sapply(1:nrow(tmp.L.exp), function(i) {
    tmp <- cor.test(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),method = "pearson")
    return(tmp$p.value)
  })
  
  tmp.res.SCC <- sapply(1:nrow(tmp.L.exp), function(i) {
    cor(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),method = "spearman")
  })
  tmp.res.SCC.p <- sapply(1:nrow(tmp.L.exp), function(i) {
    tmp <- cor.test(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),method = "spearman")
    return(tmp$p.value)
  })
  
  tmp.res.KCC <- sapply(1:nrow(tmp.L.exp), function(i) {
    cor(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),method = "kendall")
  })
  tmp.res.KCC.p <- sapply(1:nrow(tmp.L.exp), function(i) {
    tmp <- cor.test(as.numeric(tmp.L.exp[i,]), as.numeric(tmp.R.exp[i,]),method = "kendall")
    return(tmp$p.value)
  })
  
  tmp.res.rate.L <- sapply(1:nrow(tmp.L.exp), function(i) sum(as.numeric(tmp.L.exp[i,])>0)/length(as.numeric(tmp.L.exp[i,])>0))
  tmp.res.rate.R <- sapply(1:nrow(tmp.R.exp), function(i) sum(as.numeric(tmp.R.exp[i,])>0)/length(as.numeric(tmp.R.exp[i,])>0))
  tmp.res.cor.1 <- data.frame(Lgene=Lgene.list,
                              Rgene=Rgene.list,
                              LRpairs=paste0(Lgene.list,"_",Rgene.list),
                              PCC=tmp.res.PCC,
                              PCC.p = tmp.res.PCC.p,
                              SCC = tmp.res.SCC,
                              SCC.p = tmp.res.SCC.p,
                              KCC = tmp.res.KCC,
                              KCC.p = tmp.res.KCC.p,
                              detection.rate.L=tmp.res.rate.L,
                              detection.rate.R=tmp.res.rate.R,
                              stringsAsFactors = F) %>%
    dplyr::filter(detection.rate.L > tmp.cell.cutoff/cell_number & detection.rate.R > tmp.cell.cutoff/cell_number)
  
  tmp.res.cor.1$PCC[is.na(tmp.res.cor.1$PCC)] <- 0 
  tmp.res.cor.1 <- tmp.res.cor.1 %>%
    dplyr::filter(PCC.p < tmp.p.cutoff) %>%
    arrange(-PCC) 
  
  return(tmp.res.cor.1)
  
}







my_LRScatterPlot <- function(mat = data.plot,
                          mat.meta = tmp.meta,
                          tmp.Lgene = tmp.res.cor.1$Lgene[tmp.idx],
                          tmp.Rgene = tmp.res.cor.1$Rgene[tmp.idx],
                          PCC = tmp.res.cor.1$PCC[tmp.idx]){
  mat.L <- mat[tmp.Lgene,]
  mat.R <- mat[tmp.Rgene,]
  
  res.df <- data.frame(Lgene = mat.L,
                       Rgene = mat.R,
                       stringsAsFactors = F)
  colnames(res.df) <- c(tmp.Lgene,tmp.Rgene)
  res.df <- res.df %>%
    rownames_to_column("cell_id") %>%
    left_join(mat.meta,by = "cell_id") %>%
    gather(key = "gene",value = "gene.exp",-cell_id,-Time) %>%
    mutate(gene=factor(gene,levels = c(tmp.Lgene,tmp.Rgene)))
  
  
  
  p <- ggplot(data = res.df,aes_string(x="Time",y="gene.exp",colour="gene"))+
    geom_point(alpha=0.1,size=1)+
    geom_smooth(aes_string(group="gene"),se = F)+
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



#### Define function
myDavid_GOplot <- function(df=tmp,
                           use.FDR=F,
                           FDR.cut=0.01,
                           p.cut=0.01,
                           show_number=20,
                           title="LR_cor_ge_0.1_Lgene",
                           fill.color="#E64B35FF",
                           plot.ylab="BP",
                           term.pos.adjust=0.8,
                           font.size = 12){
  
  data.plot <- df
  
  data.plot$logp <- -log10(data.plot$PValue)
  data.plot$logFDR <- -log10(data.plot$FDR)
  data.plot$logp.adjust <- -log10(data.plot$Benjamini)
  
  data.plot$Term <- gsub(pattern = ".*~",replacement = "",data.plot$Term)
  
  if(use.FDR){
    data.plot <- data.plot %>%
      filter(FDR < FDR.cut) %>%
      top_n(n = show_number,wt = logFDR)
  }else{
    data.plot <- data.plot %>%
      filter(PValue < p.cut) %>%
      top_n(n = show_number,wt = logp)
  }
  
  data.plot$Term <- factor(data.plot$Term,levels = rev(data.plot$Term))
  data.term <- data.frame(x = 0,
                          y = row_number(data.plot$Term),
                          label = data.plot$Term,
                          stringsAsFactors = F)
  p1 <- ggplot() + 
    geom_bar(data = data.plot, aes_string(x = "Term", 
                                          y = ifelse(use.FDR,"logFDR","logp")),
             stat = "identity",
             fill=fill.color)+
    coord_flip() +
    ggtitle(title)+
    xlab(plot.ylab)+
    ylab(ifelse(use.FDR,"-log10(FDR)","-log10(Pvalue)"))+
    geom_text(data = data.term, 
              aes(x=y,y=x,label = label), 
              color = "black",
              hjust = term.pos.adjust,
              size = 5.8)+
    theme_cowplot(font_size = font.size) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust=0.5))
  
  return(p1)
}

#### myenrichrGo_plot
myenrichr_GOplot <- function(df=tmp,
                           p.adjust.cut=0.01,
                           show_number=20,
                           title="LR_cor_ge_0.1_Lgene",
                           fill.color="#E64B35FF",
                           plot.ylab="BP",
                           term.pos.adjust=0.8,
                           font.size = 12){
  data.plot <- df
  data.plot$logp.adjust <- -log10(data.plot$p.adjust)
  data.plot$Term <- data.plot$Description
  data.plot <- data.plot %>%
    dplyr::filter(p.adjust < p.adjust.cut) %>%
    slice_max(order_by = logp.adjust,n=show_number,with_ties=F)
    
    
  data.plot$Term <- factor(data.plot$Term,levels = rev(data.plot$Term))
  data.term <- data.frame(x = 0,
                          y = row_number(data.plot$Term),
                          label = data.plot$Term,
                          stringsAsFactors = F)
  p1 <- ggplot() + 
    geom_bar(data = data.plot, aes_string(x = "Term", 
                                          y = "logp.adjust"),
             stat = "identity",
             fill=fill.color)+
    coord_flip() +
    ggtitle(title)+
    xlab(plot.ylab)+
    ylab(expression(-log[10](p.just)))+
    geom_text(data = data.term, 
              aes(x=y,y=x,label = label), 
              color = "black",
              hjust = term.pos.adjust,
              size = 5.8)+
    theme_cowplot(font_size = font.size) +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(hjust=0.5))
  
  return(p1)
}





myCreateSeurat <- function (data=data, project = "SeuratProject", assay = "RNA", min.cells = 0, 
                            min.features = 0, names.field = 1, names.delim = "_", meta.data = NULL) {
  if (!is.null(x = meta.data)) {
    if (is.null(x = rownames(x = meta.data))) {
      stop("Row names not set in metadata. Please ensure that rownames of metadata match column names of data matrix")
    }
    if (length(x = setdiff(x = rownames(x = meta.data), y = colnames(x = data)))) {
      warning("Some cells in meta.data not present in provided data matrix.")
      meta.data <- meta.data[intersect(x = rownames(x = meta.data), 
                                       y = colnames(x = data)), ]
    }
    if (is.data.frame(x = meta.data)) {
      new.meta.data <- data.frame(row.names = colnames(x = data))
      for (ii in 1:ncol(x = meta.data)) {
        new.meta.data[rownames(x = meta.data), colnames(x = meta.data)[ii]] <- meta.data[, 
                                                                                         ii, drop = FALSE]
      }
      meta.data <- new.meta.data
    }
  }
  assay.data <- CreateAssayObject(data=as.matrix(data))
  Key(object = assay.data) <- paste0(tolower(x = assay), "_")
  assay.list <- list(assay.data)
  names(x = assay.list) <- assay
  init.meta.data <- data.frame(row.names = colnames(x = assay.list[[assay]]))
  idents <- factor(x = unlist(x = lapply(X = colnames(x = assay.data), 
                                         FUN = Seurat:::ExtractField, field = names.field, delim = names.delim)))
  if (any(is.na(x = idents))) {
    warning("Input parameters result in NA values for initial cell identities. Setting all initial idents to the project name")
  }
  ident.levels <- length(x = unique(x = idents))
  if (ident.levels > 100 || ident.levels == 0 || ident.levels == 
      length(x = idents)) {
    idents <- rep.int(x = factor(x = project), times = ncol(x = assay.data))
  }
  names(x = idents) <- colnames(x = assay.data)
  object <- new(Class = "Seurat", assays = assay.list, meta.data = init.meta.data, 
                active.assay = assay, active.ident = idents, project.name = project, 
                version = packageVersion(pkg = "Seurat"))
  object[["orig.ident"]] <- idents
  n.calc <- Seurat:::CalcN(object = assay.data)
  if (!is.null(x = n.calc)) {
    names(x = n.calc) <- paste(names(x = n.calc), assay, 
                               sep = "_")
    object[[names(x = n.calc)]] <- n.calc
  }
  if (!is.null(x = meta.data)) {
    object <- AddMetaData(object = object, metadata = meta.data)
  }
  return(object)
}

myRemoveNA <- function(df){
  df[is.na(df)] <- 0
  return(df)
}


myNormalize <- function(data.exp){
  data.exp <- log2(data.exp+1)
  return(data.exp)
}

myMatrixFilter <- function(data.exp,filter=0.95){
  data.exp <- apply(data.exp,2,FUN = function(x) ifelse(x >= quantile(x,filter),quantile(x,filter),x))
  return(data.exp)
}

myBoxPlot <- function(data.plot,x,y,xlab,ylab,...){
  p <- ggplot(data.plot,aes_string(x,y,...))+
    geom_boxplot()+
    theme_cowplot(font_size = 15)+
    xlab(xlab)+
    ylab(ylab)+
    theme(legend.position = "right",
          axis.line.x = element_line( size = 1 ),
          axis.line.y = element_line( size = 1 ),
          axis.text.x = element_text(hjust = 0.8,
                                     angle = 30))
  return(p)
}

### file generate

myFileName <- function(prefix,suffix){
  res <- paste0(prefix,"_",format(Sys.time(),"%Y%m%d%H"),suffix)
  return(res)
}

### Capitlized

myCapitalized <- function(tmp.string){
  tmp.res <- Hmisc::capitalize(tolower(tmp.string))
}





myrunDAVID <- function (gene=eg$ENTREZID, 
                        idType = "ENTREZ_GENE_ID", 
                        annotation = "GOTERM_BP_DIRECT",
                        david.user="longtengwang@pku.edu.cn") {
  
  david.pkg <- "RDAVIDWebService"
  pkgs <- installed.packages()[, 1]
  if (!david.pkg %in% pkgs) {
    stop("You should have RDAVIDWebService package installed before using enrichDAVID...")
  }
  require(david.pkg, character.only = TRUE)
  
  david <- DAVIDWebService$new(email = david.user, 
                               url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")
  idType <- match.arg(idType, getIdTypes(david))
  david.res <- addList(david, gene, idType = idType, listName = "clusterProfiler", 
                       listType = "Gene")
  if (david.res$inDavid == 0) {
    stop("All id can not be mapped. Please check 'idType' parameter...")
  }
  
  setAnnotationCategories(david, annotation)
  x <- getFunctionalAnnotationChart(david)
  if (length(x@.Data) == 0) {
    warning("No significant enrichment found...")
    return(NULL)
  }
  term.df <- x@.Data
  names(term.df) <- x@names
  term.df <- as.data.frame(term.df)
  return(term.df)
}

myPerformDavidGO <- function(gene.lits,output.table,output.fig,plot.title,annotation){
  eg = bitr(gene.lits, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)
  term.df = myrunDAVID(gene = eg$ENTREZID, 
                       idType="ENTREZ_GENE_ID",
                       annotation=annotation,
                       david.user = "longtengwang@pku.edu.cn")
  term.df$Genes <- as.character(term.df$Genes)
  
  test <- lapply(term.df$Genes,function(x){
    eg <- bitr(unlist(strsplit(x, split = ", ")),
               fromType = "ENTREZID",
               toType = "SYMBOL",
               OrgDb = org.Mm.eg.db)
    res <- paste(eg$SYMBOL,collapse = ", ")
    return(res)
  })
  term.df$Genes <- unlist(test)
  write.table(term.df,file = output.table,sep = "\t",row.names = F,quote = F)
  
  p <- myDavid_GOplot(df = term.df,
                      use.FDR = F,
                      p.cut = 0.01,
                      plot.ylab = annotation,
                      title = plot.title,
                      show_number = 20)
  ggsave(p,filename = output.fig,width = 8,height = 8,dpi = 300)
}


myMatToZeroOne <- function(mat,cutoff){
  mat <- ifelse(mat > cutoff,1,0)
  return(mat)
}

mutate_cond <- function(.data, condition, ..., envir = parent.frame()) {
  condition <- eval(substitute(condition), .data, envir)
  .data[condition, ] <- .data[condition, ] %>% mutate(...)
  .data
}

myMatCutoff <- function(df,cutoff){
  res <- df[apply(df,1,function(x) all(x>cutoff)),]
  return(res)
}

median_center <- function(x){
  res <- x - apply(x, 2, median)
  return(res)
}


myPerformDavidGO <- function(gene.lits,output.table,output.fig,plot.title,annotation){
  eg = bitr(gene.lits, fromType="SYMBOL", toType="ENTREZID", OrgDb= org.Mm.eg.db)
  term.df = myrunDAVID(gene = eg$ENTREZID, 
                       idType="ENTREZ_GENE_ID",
                       annotation=annotation,
                       david.user = "longtengwang@pku.edu.cn")
  term.df$Genes <- as.character(term.df$Genes)
  
  test <- lapply(term.df$Genes,function(x){
    eg <- bitr(unlist(strsplit(x, split = ", ")),
               fromType = "ENTREZID",
               toType = "SYMBOL",
               OrgDb = org.Mm.eg.db)
    res <- paste(eg$SYMBOL,collapse = ", ")
    return(res)
  })
  term.df$Genes <- unlist(test)
  write.table(term.df,file = output.table,sep = "\t",row.names = F,quote = F)
  
  p <- myDavid_GOplot(df = term.df,
                      use.FDR = F,
                      p.cut = 0.01,
                      plot.ylab = annotation,
                      title = plot.title,
                      show_number = 20)
  ggsave(p,filename = output.fig,width = 8,height = 8,dpi = 300)
}

myBoxplot_advanced <- function(data = df,
                               x = "dose",
                               y = "len",
                               fill = "supp",
                               color = "supp",
                               mycolor = mypalette,
                               boxwidth = 0.6,
                               errorbarwidth = 0.25,
                               dodgingwidth = 0.8,
                               ylimts = c(0,40),
                               ybreaks = seq(0,40,10),
                               fontsize = 18,
                               linesize = 1,
                               medianlinesize = 1,
                               axislinesize = 1,
                               outlier.shape = 19,
                               outlier.alpha = 0.1){
  bxp <- ggplot(data = data,aes_string(x = x,
                                       y = y,
                                       fill = fill,
                                       color= color))+
    geom_boxplot(linetype = "dashed",
                 position = position_dodge(dodgingwidth),
                 width = boxwidth, 
                 size = linesize,
                 outlier.shape = outlier.shape,
                 outlier.alpha = outlier.alpha) +
    stat_boxplot(geom = "errorbar", 
                 aes(ymin = ..ymax..),
                 width = errorbarwidth , 
                 size = linesize,
                 position = position_dodge(dodgingwidth)) +
    stat_boxplot(geom = "errorbar", 
                 aes(ymax = ..ymin..),
                 width = errorbarwidth ,
                 position = position_dodge(dodgingwidth)) +
    stat_boxplot(aes(ymin = ..lower.., 
                     ymax = ..upper..),
                 outlier.alpha = outlier.alpha,
                 outlier.shape = outlier.shape,
                 width = boxwidth,
                 size = linesize,
                 position = position_dodge(dodgingwidth)) +
    ### The bang-bang operator !! forces a single object. One common case for  !! is to substitute an environment-variable (created with <-) with a data-variable (inside a data frame).
    stat_boxplot(data = data,geom = "errorbar", 
                 aes(x = !!sym(x),
                     y = !!sym(y),
                     ymin = ..middle..,
                     ymax = ..middle..),
                 size = medianlinesize,
                 col="black",
                 width = boxwidth ,
                 position = position_dodge(dodgingwidth))
  if(is.null(ylimts)){
    bxp <- bxp + scale_fill_manual(values = mycolor)+
      scale_color_manual(values =mycolor)+
      theme_cowplot(font_size = fontsize)+
      theme(legend.position = "right",
            axis.line.x = element_line( size = axislinesize ),
            axis.line.y = element_line( size = axislinesize ),
            axis.text.x = element_text(hjust = 0.8,
                                       angle = 30))
  }else{
    bxp <- bxp + scale_y_continuous(limits = ylimts,breaks = ybreaks,expand = expansion(mult = c(0,0.1)))+
      scale_fill_manual(values = mycolor)+
      scale_color_manual(values = mycolor)+
      theme_cowplot(font_size = fontsize)+
      theme(legend.position = "right",
            axis.line = element_line( size = axislinesize ),
            axis.ticks = element_line( size = axislinesize ),
            axis.text.x = element_text(hjust = 0.8,
                                       angle = 30))
  }
  return(bxp)
}

#### about color bar

# code from
# https://stackoverflow.com/questions/9314658/colorbar-from-custom-colorramppalette
# edited
# Function to plot color bar
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab=NA, yaxt='n', ylab=NA, main=title)
  #axis(side = 4, ticks, las=1)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}

# code from
# https://stackoverflow.com/questions/15006211/how-do-i-generate-a-mapping-from-numbers-to-colors-in-r 
map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}

#### myshowcolorbar

my_showcolorbar <- function(tmp.colorbar,tmp.name,n=100){
  image(
    1:n, 1, as.matrix(1:n),
    col = tmp.colorbar(n),
    xlab = tmp.name, 
    ylab = "", 
    xaxt = "n", 
    yaxt = "n", 
    bty = "n"
  )
}

#####started to edit this Utils function
#####202108
#####20210927
#####20211023

#####---------1. Operation function--------

myMinMaxScale <- function(x){
  (x-min(x))/(max(x)-min(x))
}



myGetGRangesFromsTxDb <- function(
  txdb = mm9KG_txdb,
  standard.chromosomes = TRUE,
  verbose = TRUE
) {
  if (!requireNamespace("biovizBase", quietly = TRUE)) {
    stop("Please install biovizBase\n",
         "https://www.bioconductor.org/packages/biovizBase/")
  }
  # convert seqinfo to granges
  whole.genome <-  as(object = seqinfo(x = txdb), Class = "GRanges")
  if (standard.chromosomes) {
    whole.genome <- keepStandardChromosomes(whole.genome, pruning.mode = "coarse")
  }
  # extract genes from each chromosome
  if (verbose) {
    ####crunch featching GRanges from various data source
    tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
      cat(seqlevels(whole.genome)[x],sep = "\n")
      biovizBase::crunch(
        obj = txdb,
        which = whole.genome[x],
        columns = c("tx_id", "tx_name", "gene_id"))
    })
  } else {
    tx <- sapply(X = seq_along(whole.genome), FUN = function(x){
      suppressMessages(expr = biovizBase::crunch(
        obj = txdb,
        which = whole.genome[x],
        columns = c("tx_id", "tx_name","gene_id")))
    })
  }
  # combine
  tx <- do.call(what = c, args = tx)
  return(tx)
}


p2star <- function(p){
  symnum(p,cutpoints = c(0,0.001,0.01,0.05,1),
         symbols = c('***','**','*',''),na = NA)
}

my_enrichment.pvalue <- function(M,N,k,n,lower.tail=F){
  p <- phyper(k-1,M,N-M,n,lower.tail = lower.tail)
  return(p)
}
myTrunclog <- function(x){
  res <- ifelse(x>1,log10(x),0)
}

####extract and build Cellchat object
mySeuratToCellChat <- function(tmp.seu = NULL,tmp.group.by="cell.type"){
  tmp.mat <- GetAssayData(tmp.seu,slot = "data")
  tmp.meta <- tmp.seu[[]]
  ####Create CellChat object
  ####require exp.mat to be matrix;
  cellchat <- createCellChat(object = as.matrix(tmp.mat),
                             meta = tmp.meta,
                             group.by = tmp.group.by)
  return(cellchat)
}




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


#####---------2. IO--------

myBegraphToBigwig <- function(tmp.input,tmp.output,region=NULL,tmp.seqlength=NULL,tmp.levels=NULL){
  cat("load bedgraph",sep = "\n")
  if(!is.null(region)){
    region_data <- rtracklayer::import.bedGraph(
      con = tmp.input,
      which = as(region,"GRanges")
    )
  }else{
    region_data <- rtracklayer::import.bedGraph(
      con = tmp.input,
    )
  }
  #rtracklayer::close(con = tmp.input)
  seqlevels(region_data) <- tmp.levels
  seqlengths(region_data) <- tmp.seqlength
  cat("export bedgraph",sep = "\n")
  rtracklayer::export.bw(object = region_data,con = tmp.output)
  #rtracklayer::close(con = tmp.output)
}

myGet_ggsize <- function(plot,tmp.unit="in") {
  gtab <- patchwork:::patchworkGrob(plot)
  
  has_fixed_dimensions <- 
    !gtab$widths %>% map(~is.null(attr(.x, "unit"))) %>% unlist() %>% any() |
    !gtab$heights %>% map(~is.null(attr(.x, "unit"))) %>% unlist() %>% any()
  
  if (has_fixed_dimensions) {
    width <- grid::convertWidth(sum(gtab$widths) + unit(1, "mm"), unitTo = tmp.unit, valueOnly = TRUE)
    height <- grid::convertHeight(sum(gtab$heights) + unit(1, "mm"), unitTo = tmp.unit, valueOnly = TRUE)
    c(width = width, height = height)
  } else {
    c(width = NA, height = NA)
  }
}

get_pheatmap_dims <- function(heat_map){
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}

calc_ht_size = function(ht, unit = "inch") {
  pdf(NULL)
  ht = draw(ht)
  w = ComplexHeatmap:::width(ht)
  w = convertX(w, unit, valueOnly = TRUE)
  h = ComplexHeatmap:::height(ht)
  h = convertY(h, unit, valueOnly = TRUE)
  dev.off()
  c(w, h)
}

### default to save to jpg
myggsave <- function(p,prefix,suffix,width=8,height=8,...){
  if(suffix==".pdf"){
    ggsave(p,filename = myFileName(prefix = prefix,suffix = suffix),width = width,height = height,useDingbats=F,...)
  }
  else{
    ggsave(p,filename = myFileName(prefix = prefix,suffix = suffix),width = width,height = height,...)
  }
}



#####---------3. visualization--------

#####process colors
mycolor.bar <- function(my.colors,
                        min,
                        max=-min,
                        breaks=seq(min,max,
                                   length.out=length(my.colors)+1),
                        vertical=T,
                        showticks=T,
                        tmp.title=NULL){
  
  if(vertical){
    z=matrix(breaks,nrow=1)
    x=1
    y=breaks # supposing 3 and 2345 are the range of your data
    image(x,y,z,col = my.colors,
          breaks = breaks,
          axes = FALSE,xlab="",ylab="")
    axis(4,tick = showticks,lwd.ticks = 1)
    title(main = tmp.title)
  }
  else{
    x = breaks
    y = 1
    z=matrix(breaks,ncol=1)
    image(x,y,z,col=my.colors,
          breaks = breaks,
          axes=FALSE,xlab="",ylab="")
    axis(1,tick = showticks,lwd.ticks = 1)
    title(main = tmp.title)
  }
}


myShowColors <- function (colours, labels = TRUE, borders = NULL, cex_label = 1, 
                          ncol = NULL,tmp.title = NULL) {
  n <- length(colours)
  ncol <- ncol %||% ceiling(sqrt(length(colours)))
  nrow <- ceiling(n/ncol)
  colours <- c(colours, rep(NA, nrow * ncol - length(colours)))
  colours <- matrix(colours, ncol = ncol, byrow = TRUE)
  old <- par(pty = "s", mar = c(0, 0, 1, 0))
  on.exit(par(old))
  size <- max(dim(colours))
  plot(c(0, size), c(0, -size), type = "n", xlab = "", 
       ylab = "", axes = FALSE, main = tmp.title)
  rect(col(colours) - 1, -row(colours) + 1, col(colours), -row(colours), 
       col = colours, border = borders)
  if (labels) {
    hcl <- farver::decode_colour(colours, "rgb", "hcl")
    label_col <- ifelse(hcl[, "l"] > 50, "black", 
                        "white")
    text(col(colours) - 0.5, -row(colours) + 0.5, colours, 
         cex = cex_label, col = label_col)
  }
}

myBedTrack <- function(bed.input = tmp.files,
                       tmp.region,
                       track.color,
                       font_size = 18,
                       tmp.label = "2cell"){
  tmp.df <- read.delim(file = bed.input,header = F)
  tmp.region <- as(tmp.region,"GRanges")
  names(tmp.df)[1:3] <- c("chr", "xmin", "xmax")
  tmp.df.plot <- tmp.df %>%
    dplyr::filter(chr == as.character(chrom(tmp.region)))
  
  p <- ggplot()+
    geom_rect(ggplot2::aes(xmin = xmin, ymin = 0, 
                           xmax = xmax, ymax = 1),
              fill = track.color,
              data = tmp.df.plot)+
    ylab(tmp.label)+
    scale_x_continuous(limits = c(start(tmp.region),end(tmp.region)),
                       expand = c(0,0))+
    scale_y_continuous(expand = c(0,0))+
    theme_cowplot(font_size = font_size)+
    theme(axis.line.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_line(size = 0.5),
          axis.title.y = element_text(angle = 0,vjust = 0.5),
          axis.ticks.y = element_blank())
  return(p)
}


myBigwigTrack <- function(
  region = tmp.region,
  bigwig = K4me3,
  smooth = 200,
  lognorm = F,
  type = "coverage",
  y_label = "Score",
  fontsize=18,
  track.color="blue",
  tmp.ylimits=c(0,16),
  max.downsample = 3000,
  downsample.rate = 0.1,
  tmp.seed=42
) {
  possible_types <- c("line", "heatmap", "coverage","bar")
  if (!(type %in% possible_types)) {
    stop(
      "Invalid type requested. Choose ",
      paste(possible_types, collapse = ", ")
    )
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/")
    return(NULL)
  }
  if (!inherits(x = region, what = "GRanges")) {
    region <- as(region,"GRanges")
  }
  ### to smooth track data
  region_data <- rtracklayer::import(
    con = bigwig,
    which = region,
    as = "NumericList"
  )[[1]]
  if (!is.null(x = smooth)) {
    region_data <- RcppRoll::roll_mean(x = region_data, n = smooth, fill = 0L)
  }
  region_data <- data.frame(
    position = start(x = region):end(x = region),
    score = region_data,
    stringsAsFactors = FALSE
  )
  window.size = width(x = region)
  sampling <- min(max.downsample, window.size * downsample.rate)
  set.seed(tmp.seed)
  coverages <- slice_sample(.data = region_data, n = sampling)
  if(lognorm){
    coverages$score <- myTrunclog(coverages$score)
  }
  p <- ggplot(
    data = coverages,
    mapping = aes_string(x = "position", y = "score")
  )
  if (type == "line") {
    p <- p + geom_line()
  } else if (type == "heatmap") {
    # different downsampling needed for heatmap
    # cut into n bins and average within each bin
    region_data$bin <- floor(x = region_data$position / smooth)
    region_data <- group_by(region_data, bin)
    region_data <- mutate(region_data, score = mean(x = score))
    region_data <- ungroup(region_data)
    region_data <- unique(x = region_data[, c("bin", "score")])
    p <- ggplot(
      data = region_data,
      mapping = aes_string(x = "bin", y = 1, fill = "score")
    ) + geom_tile() + scale_fill_viridis_c()
  } else if (type == "coverage") {
    p <- p + 
      geom_area(fill=track.color,color=track.color)
  } else if (type == "bar") {
    p <- p + 
      geom_bar(stat="identity",fill=track.color,color=track.color)
  } 
  chromosome <- as.character(x = seqnames(x = region))
  p <- p  +
    xlab(label = paste0(chromosome, " position (bp)")) +
    ylab(label = y_label) +
    scale_x_continuous(expand = c(0,0),limits = c(start(region),end(region)))+
    scale_y_continuous(expand = c(0,0),limits = tmp.ylimits,
                       breaks = c(tmp.ylimits[2]),
                       labels = paste0("[",0,"-",tmp.ylimits[2],"]"))+
    theme_cowplot(font_size = fontsize)+
    theme(axis.line = element_line(size = 0.5),
          axis.ticks = element_line(size = 0.5),
          axis.title.y = element_text(angle = 0),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank())
  return(p)
}


split_body <- function(df, width = 1000) {
  wd <- df$end - df$start
  nbreak <- wd / width
  if (nbreak > 1) {
    steps <- 0:(nbreak)
    starts <- (width * steps) + df$start
    starts[starts > df$end] <- NULL
  } else {
    starts <- df$end
  }
  breaks <- data.frame(
    seqnames = df$seqnames[[1]],
    start = starts,
    end = starts + 1,
    strand = df$strand[[1]],
    gene_name = df$gene_name[[1]],
    type = "arrow"
  )
  return(breaks)
}

record_overlapping <- function(annotation = gene_bodies, min.gapwidth = 1000) {
  # convert back to granges
  annotation.stash <- annotation
  annotation$strand <- "*"
  gr <- makeGRangesFromDataFrame(
    df = annotation[annotation$type == "body", ], keep.extra.columns = TRUE
  )
  # work out which ranges overlap
  collapsed <- GenomicRanges::reduce(
    x = gr, with.revmap = TRUE, min.gapwidth = min.gapwidth
  )$revmap
  idx <- seq_along(gr)
  for (i in seq_along(collapsed)) {
    mrg <- collapsed[[i]]
    for (j in seq_along(mrg)) {
      idx[[mrg[[j]]]] <- j
    }
  }
  names(x = idx) <- gr$gene_name
  return(idx)
}

reformat_annotations <- function(
  annotation = annotation.subset,
  start.pos = start.pos,
  end.pos = end.pos,
  arrow_sbreaks=2000
) {
  annotation <- annotation[annotation$type == "exon"]
  exons <- as.data.frame(x = annotation)
  annotation <- split(
    x = annotation,
    f = as.character(annotation$gene_id)
  )
  annotation <- lapply(X = annotation, FUN = as.data.frame)
  
  # add gene total start / end
  gene_bodies <- lapply(seq_along(annotation), function(i){
    df <- data.frame(
      seqnames = annotation[[i]]$seqnames[[1]],
      start = min(annotation[[i]]$start),
      end = max(annotation[[i]]$end),
      strand = annotation[[i]]$strand[[1]],
      gene_name = annotation[[i]]$gene_id[[1]],
      type = "body"
    )
    # trim any that extend beyond region
    df$start <- ifelse(
      test = df$start < start.pos,
      yes = start.pos,
      no = df$start
    )
    df$end <- ifelse(
      test = df$end > end.pos,
      yes = end.pos,
      no = df$end
    )
    breaks <- split_body(df = df,width = arrow_sbreaks)
    df <- rbind(df, breaks)
    return(df)
  })
  
  gene_bodies <- do.call(what = rbind, args = gene_bodies)
  
  # record if genes overlap
  overlap_idx <- record_overlapping(annotation = gene_bodies, min.gapwidth = 1000)
  gene_bodies$dodge <- overlap_idx[as.character(gene_bodies$gene_name)]
  exons$dodge <- overlap_idx[as.character(exons$gene_id)]
  
  label_df <- gene_bodies[gene_bodies$type == "body", ]
  label_df$width <- label_df$end - label_df$start
  label_df$position <- label_df$start + (label_df$width / 2)
  
  onplus <- gene_bodies[gene_bodies$strand %in% c("*", "+"), ]
  onminus <- gene_bodies[gene_bodies$strand == "-", ]
  
  return(
    list(
      "labels" = label_df,
      "exons" = exons,
      "plus" = onplus,
      "minus" = onminus
    )
  )
}



myGenePlot <- function(annotation=tx, 
                       region="chr7:152042291-152056148",showlabel=T,
                       arrow_sbreaks=2000,font_size=18,label_size=8) {
  
  if (is.null(x = annotation)) {
    return(NULL)
  }
  if (!inherits(x = region, what = "GRanges")) {
    region <- as(region,"GRanges")
  }
  start.pos <- start(x = region)
  end.pos <- end(x = region)
  chromosome <- seqnames(x = region)
  
  # get names of genes that overlap region, then subset to include only those
  # genes. This avoids truncating the gene if it runs outside the region
  annotation.subset <- subsetByOverlaps(x = annotation, ranges = region)
  genes.keep <- unique(x = annotation.subset$gene_id)
  annotation.subset <- annotation[
    fastmatch::fmatch(x = annotation$gene_id, table = genes.keep, nomatch = 0L) > 0L
  ]
  
  if (length(x = annotation.subset) == 0) {
    # make empty plot
    p <- ggplot(data = data.frame())
    y_limit <- c(0, 1)
  } else {
    annotation_df_list <- reformat_annotations(
      annotation = annotation.subset,
      start.pos = start.pos,
      end.pos = end.pos,
      arrow_sbreaks=arrow_sbreaks
    )
    p <- ggplot() +
      # exons
      geom_segment(
        data = annotation_df_list$exons,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$exons$dodge,
          xend = "end",
          yend = annotation_df_list$exons$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 5
      ) +
      # gene body
      geom_segment(
        data = annotation_df_list$labels,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$labels$dodge,
          xend = "end",
          yend = annotation_df_list$labels$dodge,
          color = "strand"
        ),
        show.legend = FALSE,
        size = 1/2
      )
    if (nrow(x = annotation_df_list$plus) > 0) {
      # forward strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$plus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$plus$dodge,
          xend = "end",
          yend = annotation_df_list$plus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "last",
          type = "open",
          angle = 45,
          length = unit(x = 0.05, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    if (nrow(x = annotation_df_list$minus) > 0) {
      # reverse strand arrows
      p <- p + geom_segment(
        data = annotation_df_list$minus,
        mapping = aes_string(
          x = "start",
          y = annotation_df_list$minus$dodge,
          xend = "end",
          yend = annotation_df_list$minus$dodge,
          color = "strand"
        ),
        arrow = arrow(
          ends = "first",
          type = "open",
          angle = 45,
          length = unit(x = 0.05, units = "inches")
        ),
        show.legend = FALSE,
        size = 1/2
      )
    }
    # label genes
    n_stack <- max(annotation_df_list$labels$dodge)
    annotation_df_list$labels$dodge <- annotation_df_list$labels$dodge + (n_stack * 0.2)
    if(showlabel){
      p <- p + geom_text(
        data = annotation_df_list$labels,
        mapping = aes_string(x = "position", y = "dodge", label = "gene_name"),
        size = label_size,fontface="bold"
      )
    }
    y_limit <- c(0.9, n_stack + (n_stack * 0.5))
  }
  p <- p +
    theme_cowplot(font_size = font_size) +
    ylab("Genes") +
    xlab(label = paste0(chromosome, " position (bp)")) +
    # scale_x_continuous(expand = c(0,0),
    #                    limits = c(start.pos,end.pos),
    #                    breaks = seq(start.pos,end.pos,length.out=5))+
    scale_x_continuous(expand = c(0,0),
                       limits = c(start.pos,end.pos))+
    ylim(y_limit) +
    theme(
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      axis.line.y = element_blank()
    ) +
    scale_color_manual(values = c("darkblue", "darkgreen"))
  return(p)
}


get_plot_dims <- function(heat_map){
  plot_height <- sum(sapply(heat_map$gtable$heights, grid::convertHeight, "in"))
  plot_width  <- sum(sapply(heat_map$gtable$widths, grid::convertWidth, "in"))
  return(list(height = plot_height, width = plot_width))
}

# myFunction <- function(p){
#   tmp.range.x <- round(range(p$data$UMAP_1))-1
#   tmp.range.y <- round(range(p$data$UMAP_2))-1
#   tmp.ticks.x <- seq(tmp.range.x[1],tmp.range.x[2],length.out=20)
#   tmp.ticks.y <- seq(tmp.range.y[1],tmp.range.y[2],length.out=20)
#   list(geom_path(data = data.frame(x=tmp.ticks.x[1:2],
#                                    y=tmp.ticks.y[c(1,1)]),
#                  size = 1,
#                  mapping = aes(x,y),
#                  linejoin = "bevel", 
#                  lineend = "round",
#                  arrow = arrow(angle = 45,type = "closed",length = unit(0.1, "inches"))),
#        geom_path(data = data.frame(x=tmp.ticks.x[c(1,1)],
#                                    y=tmp.ticks.y[1:2]),
#                  size = 1,
#                  mapping = aes(x,y),
#                  linejoin = "bevel", 
#                  lineend = "round",
#                  arrow = arrow(angle = 45,type = "closed",length = unit(0.1, "inches"))),
#        annotate("text",x = mean(tmp.ticks.x[1:2]),
#                 y = 2*tmp.ticks.y[1]-tmp.ticks.y[2],label="UMAP1",fontface="bold",size=4),
#        annotate("text",x = 2*tmp.ticks.x[1]-tmp.ticks.x[2],
#                 y = mean(tmp.ticks.y[1:2]),label="UMAP2",fontface="bold",size=4,angle=90),
#        NoAxes())
# }

myFeatureplot <- function(seu=pbmc3k.final,dim.use="tSNE",font.size=14,pt.size=3,
                          features=tmp.markers,tmp.color=coolwarm(100),tmp.ncol=2){
  plot.list <- lapply(features, function(ii){
    p <- FeaturePlot(seu,features = ii,pt.size = pt.size,raster = T) 
    tmp.data.plot <- p$data %>%
      arrange(!!sym(eval(ii)))
    p <- ggplot(data = tmp.data.plot)+
      scattermore::geom_scattermore(mapping = aes(!!sym(paste0(dim.use,"_1")),
                                                  !!sym(paste0(dim.use,"_2")),
                                                  color=!!sym(eval(ii))),
                                    pointsize = 3)+
      scale_color_gradientn(colors  = tmp.color,name=NULL)+
      ggtitle(ii)+
      theme_cowplot(font_size = font.size)
    return(p)
  })
  wrap_plots(plot.list,ncol = tmp.ncol)
}

myAddSeuratAxis <- function(p,reduction="tsne",tmp.ratio=0.3,tmp.fontsize=6){
  if (reduction == "tsne") {
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
         theme(plot.title = element_text(size = 24),
               legend.text = element_text(size = 24)))
  }
  
  if(reduction == "umap"){
    tmp.range.x <- round(range(p$data$UMAP_1))-1
    tmp.range.y <- round(range(p$data$UMAP_2))-1
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
                  y = tmp.ticks.y[1]-tmp.ratio*tmp.step.y,
                  label="UMAP_1",fontface="bold",size=tmp.fontsize),
         annotate("text",x = tmp.ticks.x[1]-tmp.ratio*tmp.step.x,
                  y = mean(tmp.ticks.y[1:2]),
                  label="UMAP_2",fontface="bold",size=tmp.fontsize,angle=90),
         NoAxes(),
         theme(plot.title = element_text(size = 24),
               legend.text = element_text(size = 24)))
  }
}



pheatmap_fixed <- function (mat, color = colorRampPalette(rev(brewer.pal(n = 7, 
                                                                     name = "RdYlBu")))(100), kmeans_k = NA, breaks = NA, 
                        border_color = ifelse(nrow(mat) < 100 & ncol(mat) < 100, 
                                              "grey60", NA), cellwidth = NA, cellheight = NA, 
                        scale = "none", cluster_rows = TRUE, cluster_cols = TRUE, 
                        clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", 
                        clustering_method = "complete", clustering_callback = NA, 
                        cutree_rows = NA, cutree_cols = NA, treeheight_row = ifelse(class(cluster_rows) == 
                                                                                      "hclust" || cluster_rows, 50, 0), treeheight_col = ifelse(class(cluster_cols) == 
                                                                                                                                                  "hclust" || cluster_cols, 50, 0), legend = TRUE, 
                        legend_breaks = NA, legend_labels = NA, annotation_row = NA, 
                        annotation_col = NA, annotation = NA, annotation_colors = NA, 
                        annotation_legend = TRUE, annotation_names_row = TRUE, annotation_names_col = TRUE, 
                        drop_levels = TRUE, show_rownames = TRUE, show_colnames = TRUE, 
                        main = NA, fontsize = 10, fontsize_row = fontsize, fontsize_col = fontsize, 
                        angle_col = c("270", "0", "45", "90", 
                                      "315"), display_numbers = FALSE, number_format = "%.2f", 
                        number_color = "grey30", fontsize_number = 0.8 * fontsize, 
                        gaps_row = NULL, gaps_col = NULL, labels_row = NULL, labels_col = NULL, 
                        filename = NA, width = NA, height = NA, silent = FALSE, na_col = "#DDDDDD", 
                        name = NULL, heatmap_legend_param = list(), ..., run_draw = FALSE) 
{
  if (is.data.frame(mat)) {
    ComplexHeatmap:::warning_wrap("The input is a data frame, convert it to the matrix.")
    mat = as.matrix(mat)
  }
  if (!identical(kmeans_k, NA)) {
    ComplexHeatmap:::warning_wrap("argument `kmeans_k` is not suggested to use in pheatmap -> Heatmap translation because it changes the input matrix. You might check `row_km` and `column_km` arguments in Heatmap().")
    km = kmeans(mat, centers = kmeans_k)
    mat = km$centers
    rownames(mat) = paste0("Cluster: ", seq_along(km$size), 
                           ", Size: ", km$size)
  }
  if ("row" %in% scale) {
    if (any(is.na(mat))) {
      mat = (mat - rowMeans(mat, na.rm = TRUE))/rowSds(mat, 
                                                       na.rm = TRUE)
    }
    else {
      mat = t(scale(t(mat)))
    }
  }
  else if ("column" %in% scale) {
    if (any(is.na(mat))) {
      mat = t((t(mat) - colMeans(mat, na.rm = TRUE))/colSds(mat, 
                                                            na.rm = TRUE))
    }
    else {
      mat = scale(mat)
    }
  }
  ht_param = list(matrix = mat)
  if (!identical(scale, "none") && !identical(breaks, 
                                              NA)) {
    ComplexHeatmap:::warning_wrap("It not suggested to both set `scale` and `breaks`. It makes the function confused.")
  }
  if (is.function(color)) {
    ht_param$col = color
    if (!identical(breaks, NA)) {
      ComplexHeatmap:::warning_wrap("`breaks` is ignored when `color` is set as a color mapping function.")
    }
  }
  else {
    if (identical(breaks, NA)) {
      n_col = length(color)
      if (identical(scale, "row") || identical(scale, 
                                               "column")) {
        lim = max(abs(mat), na.rm = TRUE)
        ht_param$col = circlize::colorRamp2(seq(-lim, lim, length = n_col), 
                                            color)
      }
      else {
        ht_param$col = circlize::colorRamp2(seq(min(mat, na.rm = TRUE), 
                                                max(mat, na.rm = TRUE), length = n_col), color)
      }
    }
    else {
      if (length(breaks) == length(color) + 1) {
        ht_param$col = local({
          breaks = breaks
          color = color
          fun = function(x) {
            n = length(color)
            df = data.frame(start = c(-Inf, breaks[seq_len(n)], 
                                      breaks[n + 1]), end = c(breaks[1], breaks[1 + 
                                                                                  seq_len(n)], Inf))
            ind = numeric(length(x))
            for (i in seq_along(x)) {
              ind[i] = which(df$start <= x[i] & df$end > 
                               x[i])
            }
            ind = ind - 1
            ind[ind < 1] = 1
            ind[ind > n] = n
            color[ind]
          }
          attr(fun, "breaks") = breaks
          fun
        })
      }
      else if (length(breaks) == length(color)) {
        ht_param$col = circlize::colorRamp2(breaks, color)
      }
      else {
        n_col = length(color)
        ht_param$col = circlize::colorRamp2(seq(min(breaks), max(breaks), 
                                                length = n_col), color)
        ComplexHeatmap:::warning_wrap("`breaks` does not have the same length as `color`. The colors are interpolated from the minimal to the maximal of `breaks`.")
      }
    }
  }
  if (!identical(filename, NA)) {
    ComplexHeatmap:::warning_wrap("argument `filename` is not supported in pheatmap -> Heatmap translation, skip it.")
  }
  if (!identical(width, NA)) {
    ComplexHeatmap:::warning_wrap("argument `width` is not supported in pheatmap -> Heatmap translation, skip it.")
  }
  if (!identical(height, NA)) {
    ComplexHeatmap:::warning_wrap("argument `height` is not supported in pheatmap -> Heatmap translation, skip it.")
  }
  if (!identical(silent, FALSE)) {
    ComplexHeatmap:::warning_wrap("argument `silent` is not supported in pheatmap -> Heatmap translation, skip it.")
  }
  ht_param$rect_gp = gpar(col = border_color)
  if (nrow(mat) > 1000 || ncol(mat) > 1000) {
    if (!is.na(border_color)) {
      ComplexHeatmap:::warning_wrap("border color is set for the matrix with large numbers of rows or columns. You might only be able to see the border colors in the plot. Set `border_color = NA` to get rid of it.")
    }
  }
  if (!identical(cellwidth, NA)) {
    ht_param$width = ncol(mat) * unit(cellwidth, "pt")
  }
  if (!identical(cellheight, NA)) {
    ht_param$height = nrow(mat) * unit(cellheight, "pt")
  }
  if (identical(clustering_distance_rows, "correlation")) 
    clustering_distance_rows = "pearson"
  if (identical(clustering_distance_cols, "correlation")) 
    clustering_distance_cols = "pearson"
  ht_param$cluster_rows = cluster_rows
  ht_param$cluster_columns = cluster_cols
  ht_param$clustering_distance_rows = clustering_distance_rows
  ht_param$clustering_distance_columns = clustering_distance_cols
  ht_param$clustering_method_rows = clustering_method
  ht_param$clustering_method_columns = clustering_method
  if (!is.na(cutree_rows)) {
    if (inherits(cluster_rows, c("logical", "hclust", 
                                 "dendrogram"))) {
      ht_param$row_split = cutree_rows
      ht_param$row_gap = unit(4, "bigpts")
      ht_param["row_title"] = list(NULL)
    }
  }
  if (!is.na(cutree_cols)) {
    if (inherits(cluster_cols, c("logical", "hclust", 
                                 "dendrogram"))) {
      ht_param$column_split = cutree_cols
      ht_param$column_gap = unit(4, "bigpts")
      ht_param["column_title"] = list(NULL)
    }
  }
  ht_param$row_dend_width = unit(treeheight_row, "pt")
  ht_param$column_dend_height = unit(treeheight_col, "pt")
  ht_param$show_heatmap_legend = legend
  if (identical(scale, "row") || identical(scale, "column")) {
    if (identical(legend_breaks, NA)) {
      lim = quantile(abs(mat), 0.975)
      le = pretty(c(-lim, lim), n = 3)
      if (length(le) == 7 && le[1] == -3) {
        le = c(-3, -1.5, 0, 1.5, 3)
      }
      else if (!0 %in% le) {
        le = c(le[1], le[1]/2, 0, le[length(le)]/2, le[length(le)])
      }
      legend_breaks = le
    }
  }
  if (!identical(legend_breaks, NA)) {
    heatmap_legend_param$at = legend_breaks
  }
  if (!identical(legend_labels, NA)) {
    heatmap_legend_param$labels = legend_labels
  }
  ht_param$heatmap_legend_param = list(title_gp=gpar(fontsize=fontsize,fontface="bold"),
                                       labels_gp = gpar(fontsize = fontsize*0.8))
  if (identical(annotation_colors, NA)) {
    annotation_colors = list()
  }
  if (!identical(annotation_col, NA)) {
    acn = rownames(annotation_col)
    mcn = colnames(mat)
    if (!is.null(acn)) {
      if (acn[1] %in% mcn) {
        if (length(union(acn, mcn)) == length(mcn)) {
          if (!identical(acn, mcn)) {
            ComplexHeatmap:::warning_wrap("Column annotation has different order from matrix columns. Adjust the column annotation based on column names of the matrix.")
          }
          annotation_col = annotation_col[mcn, , drop = FALSE]
        }
      }
    }
    for (nm in colnames(annotation_col)) {
      if (nm %in% names(annotation_colors)) {
        if (is.null(names(annotation_colors[[nm]])) && 
            is.numeric(annotation_col[, nm])) {
          foo_x = annotation_col[, nm]
          foo_n_col = length(annotation_colors[[nm]])
          annotation_colors[[nm]] = circlize::colorRamp2(seq(min(foo_x), 
                                                             max(foo_x), length = foo_n_col), annotation_colors[[nm]])
        }
      }
    }
    ht_param$top_annotation = HeatmapAnnotation(df = annotation_col[, 
                                                                    ncol(annotation_col):1, drop = FALSE], col = annotation_colors, 
                                                show_legend = annotation_legend, show_annotation_name = annotation_names_col, 
                                                annotation_legend_param = list(title_gp=gpar(fontsize=fontsize,fontface="bold"),
                                                                               labels_gp = gpar(fontsize = fontsize*0.8)),
                                                gp = gpar(col = border_color), annotation_name_gp = gpar(fontsize = fontsize, 
                                                                                                         fontface = "bold"), simple_anno_size = unit(10, 
                                                                                                                                                     "bigpts"), gap = unit(2, "bigpts"))
  }
  if (!identical(annotation_row, NA)) {
    arn = rownames(annotation_row)
    mrn = rownames(mat)
    if (!is.null(arn)) {
      if (arn[1] %in% mrn) {
        if (length(union(arn, mrn)) == length(mrn)) {
          if (!identical(arn, mrn)) {
            ComplexHeatmap:::warning_wrap("Row annotation has different order from matrix rows. Adjust the row annotation based on row names of the matrix.")
          }
          annotation_row = annotation_row[mrn, , drop = FALSE]
        }
      }
    }
    for (nm in colnames(annotation_row)) {
      if (nm %in% names(annotation_colors)) {
        if (is.null(names(annotation_colors[[nm]])) && 
            is.numeric(annotation_row[, nm])) {
          foo_x = annotation_row[, nm]
          foo_n_col = length(annotation_colors[[nm]])
          annotation_colors[[nm]] = circlize::colorRamp2(seq(min(foo_x), 
                                                             max(foo_x), length = foo_n_col), annotation_colors[[nm]])
        }
      }
    }
    ht_param$left_annotation = rowAnnotation(df = annotation_row[, 
                                                                 ncol(annotation_row):1, drop = FALSE], col = annotation_colors, 
                                             show_legend = annotation_legend, show_annotation_name = annotation_names_row, 
                                             annotation_legend_param = list(title_gp=gpar(fontsize=fontsize,fontface="bold"),
                                                                            labels_gp = gpar(fontsize = fontsize*0.8)),
                                             gp = gpar(col = border_color), annotation_name_gp = gpar(fontsize = fontsize, 
                                                                                                      fontface = "bold"), simple_anno_size = unit(10, 
                                                                                                                                                  "bigpts"), gap = unit(2, "bigpts"))
  }
  if (!identical(annotation, NA)) {
    ComplexHeatmap:::warning_wrap("argument `annotation` is not supported in pheatmap -> Heatmap translation, skip it.")
  }
  if (identical(drop_levels, FALSE)) {
    ComplexHeatmap:::warning_wrap("argument `drop_levels` is enfored to be TRUE, skip it.")
  }
  ht_param$show_row_names = show_rownames
  ht_param$show_column_names = show_colnames
  ht_param$row_names_gp = gpar(fontsize = fontsize_row)
  ht_param$column_names_gp = gpar(fontsize = fontsize_col)
  angle_col = match.arg(angle_col)[1]
  angle_col = switch(angle_col, `0` = 0, `45` = 45, 
                     `90` = 90, `270` = 90, `315` = -45)
  ht_param$column_names_rot = angle_col
  if (angle_col == 0) {
    ht_param$column_names_centered = TRUE
  }
  if (is.logical(display_numbers)) {
    if (display_numbers) {
      ht_param$layer_fun = local({
        number_format = number_format
        number_color = number_color
        fontsize_number = fontsize_number
        mat = mat
        function(j, i, x, y, w, h, fill) {
          grid.text(sprintf(number_format, pindex(mat, 
                                                  i, j)), x = x, y = y, gp = gpar(col = number_color, 
                                                                                  fontsize = fontsize_number))
        }
      })
    }
  }
  else if (is.matrix(display_numbers)) {
    if (!identical(dim(display_numbers), dim(mat))) {
      stop_wrap("dimension of `display_numbers` should be the same as the input matrix.")
    }
    ht_param$layer_fun = local({
      number_color = number_color
      fontsize_number = fontsize_number
      mat = display_numbers
      function(j, i, x, y, w, h, fill) {
        grid.text(pindex(mat, i, j), x = x, y = y, gp = gpar(col = number_color, 
                                                             fontsize = fontsize_number))
      }
    })
  }
  if (!is.null(labels_row)) {
    ht_param$row_labels = labels_row
  }
  if (!is.null(labels_col)) {
    ht_param$column_labels = labels_col
  }
  if (!is.null(gaps_row)) {
    if (inherits(cluster_rows, c("hclust", "dendrogram"))) {
      stop_wrap("`gaps_row` should not be set when `cluster_rows` is set as a clustering object.")
    }
    if (identical(cluster_rows, TRUE)) {
      stop_wrap("`gaps_row` should not be set when `cluster_rows` is set to TRUE.")
    }
    slices = diff(c(0, gaps_row, nrow(mat)))
    ht_param$row_split = rep(seq_along(slices), times = slices)
    ht_param$row_gap = unit(4, "bigpts")
    ht_param["row_title"] = list(NULL)
  }
  if (!is.null(gaps_col)) {
    if (inherits(cluster_cols, c("hclust", "dendrogram"))) {
      stop_wrap("`gaps_col` should not be set when `cluster_cols` is set as a clustering object.")
    }
    if (identical(cluster_cols, TRUE)) {
      stop_wrap("`gaps_col` should not be set when `cluster_cols` is set to TRUE.")
    }
    slices = diff(c(0, gaps_col, ncol(mat)))
    ht_param$column_split = rep(seq_along(slices), times = slices)
    ht_param$column_gap = unit(4, "bigpts")
    ht_param["column_title"] = list(NULL)
  }
  if (!identical(clustering_callback, NA)) {
    if (!identical(ht_param$cluster_rows, FALSE)) {
      row_hclust = hclust(get_dist(mat, ht_param$clustering_distance_rows), 
                          ht_param$clustering_method_rows)
      row_hclust = clustering_callback(row_hclust, ...)
      ht_param$cluster_rows = row_hclust
    }
    if (!identical(ht_param$cluster_columns, FALSE)) {
      column_hclust = hclust(get_dist(t(mat), ht_param$clustering_distance_columns), 
                             ht_param$clustering_method_columns)
      column_hclust = clustering_callback(column_hclust, 
                                          ...)
      ht_param$cluster_columns = column_hclust
    }
  }
  ht_param$name = name
  ht_param$row_dend_reorder = FALSE
  ht_param$column_dend_reorder = FALSE
  if (!identical(main, NA)) {
    ht_param$column_title = main
    ht_param$column_title_gp = gpar(fontface = "bold", 
                                    fontsize = 1.3 * fontsize)
  }
  ht_param = c(ht_param, list(...))
  ht = do.call(Heatmap, ht_param)
  attr(ht, "translate_from") = "pheatmap"
  if (run_draw) {
    draw(ht)
  }
  else {
    ht
  }
  #return(ht_param)
}



####define function
myenrichGO_plot <- function(data.plot,
                            p.adjust.cut=0.01,
                            show_number=20,
                            title="The Most Enriched GO Terms",
                            fill.color="#8DA1CB",
                            plot.ylab="BP",
                            term.pos.adjust=0,
                            font.size = 12,
                            text.size=5.8,
                            tmp.alpha=0.8){
  
  data.plot <- tmp.df
  data.plot$logp.adjust <- -log10(data.plot$p.adjust)
  data.plot <- data.plot%>%
    dplyr::filter(p.adjust < p.adjust.cut) %>%
    slice_max(order_by = logp.adjust,
              n=show_number,
              with_ties=F)
  data.plot$Description <- factor(data.plot$Description,
                                     levels = rev(data.plot$Description))
  
  data.term <- data.frame(x = 0,
                          y = row_number(data.plot$Description),
                          label = data.plot$Description,
                          stringsAsFactors = F)
  p1 <- ggplot() + 
    geom_bar(data = data.plot, aes_string(x = "Description", 
                                          y = "logp.adjust"),
             stat = "identity",
             fill=fill.color,alpha=tmp.alpha)+
    coord_flip() +
    ggtitle(title)+
    xlab(plot.ylab)+
    ylab(expression(-log[10](p.just)))+
    geom_text(data = data.term, 
              aes(x=y,y=x,label = label), 
              color = "black",
              hjust = term.pos.adjust,
              size = text.size)+
    theme_cowplot(font_size = font.size) +
    scale_y_continuous(expand = c(0,0))+
    theme(axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(size = 1),
          axis.ticks.x = element_line(size = 1),
          plot.title = element_text(hjust=0.5))
  return(p1)
}



####Define my theme

myGridTheme <- function (font.size = 14) 
{
  theme_bw() + theme(axis.text.x = element_text(colour = "black", 
                                                size = font.size*6/7, vjust = 1), 
                     axis.text.y = element_text(colour = "black", 
                                                size = font.size, hjust = 1), 
                     axis.title = element_text(margin = margin(10, 5, 0, 0), 
                                               color = "black", size = font.size), 
                     axis.title.y = element_text(angle = 90),
                     legend.title = element_text(size = font.size),
                     legend.text = element_text(size = font.size*2/3))
}


##### draw same in base plot
# t <- data.frame(P_Value=1:10,Description=stringi::stri_rand_strings(10, 10),stringsAsFactors = F)
# stringi::stri_rand_strings(10, 5)
# t$P.Value
# barplot(as.numeric(t$P_Value),
#         horiz=T,
#         xlim=c(0,max(t$P_Value)+0.5),
#         axes=T,
#         #col="lightblue",
#         col = "pink",
#         xlab ="-log10(p-value)",
#         cex.axis=1.3,
#         cex.lab=1,
#         border = NA) 
# for (i in 1:nrow(t)){
#   text(0,(1.2*i-0.6),t$Description[i],cex=1.2,pos=4)
# }



#####---------4. defien colors --------
rdbu <- colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu"))) 
test.color.2 <- colorRampPalette(c("Orange","lightgoldenrod","darkgreen"))
col.spectral <- colorRampPalette(brewer.pal(11,'Spectral')[-6])
test.color.3 <- colorRampPalette(c("#f86e11","#e9970a","#71a701","#62b474","#71c3ac","#9fc4ca"))
rdwhbu <- colorRampPalette(c("navy", "white", "brown3"))
skyblueYelow <- colorRampPalette(c("skyblue","black","yellow"))
skybluered <- colorRampPalette(c("skyblue","black","orange"))
solarExtra <- colorRampPalette(c("#3361A5","#248AF3","#14B3FF","#88CEEF","#C1D5DC","#EAD397","#FDB31A","#E42A2A","#A31D1D"))
blues <- colorRampPalette(colors = brewer.pal(9,"Blues"))
ylord <- colorRampPalette(colors = brewer.pal(9,"YlOrRd"))
hic.red <- colorRampPalette(c("white","red"))
hic.orage <- colorRampPalette(brewer.pal(9,"Oranges"))
cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026'))
mypalette <- c(rev(cold(21)), warm(20))
coldwarm <- colorRampPalette(colors = mypalette)
coolwarm <- colorRampPalette(colors = c("#3B4CC0", "#4F6AD9", "#6585EC", "#7B9FF9", "#93B5FF", 
                                        "#AAC7FD", "#C0D4F5", "#D4DBE6", "#E5D8D1", "#F2CBB7",
                                        "#F7B89C", "#F6A081", "#EE8568",
                                        "#E0654F", "#CC4039", "#B40426"))
ramp <- colorRampPalette(c("white","pink","red","black"))

bkcolor <- c(colorRampPalette(c(brewer.pal(9,"Blues" )[4:9],"#1a1919"))(50),
             colorRampPalette(c("#1a1919",rev(brewer.pal( 9,"YlOrBr" ))[1:6]))(50))

bkcolor <- colorRampPalette(colors = bkcolor)
hic.pca.red <- colorRampPalette(c("blue","gray1","red"))
hic.pca.redwhite <- colorRampPalette(c("#1d1856","navyblue","white","red4","#861617"))
hic.pca.orange <- colorRampPalette(c("#2f2583","black","#f9b232"))
hic.pca.skyblue <- colorRampPalette(c("skyblue","black","orange"))

okabe_ito <- colorRampPalette(colors = c("#e59f01","#56b4e8","#009f73","#f0e442","#0072b1","#d55e00","#cc79a7","#999999","#000000"))
OrBl_div <- colorRampPalette(colors = c("#9f3d22","#be4d21","#db6525","#ef8531","#f1ac73","#d8d4c9","#a1bccf","#6fa3cb","#5689b6","#4171a1","#2b5b8b"))



###gundam
strike_freedom_gundam_color <- colorRampPalette(colors = c("brown3","#eff3ff","navy","#373b35","#b1b0c2","#f2be58"))

npg.color <- colorRampPalette(colors = ggsci:::ggsci_db$npg$nrc) 

#mycolor.bar(hic.pca.red(100),min = -1,max = 1)
divergentcolor <- function (n) {
  colorSpace <- c("#E41A1C", "#377EB8", "#4DAF4A", 
                  "#984EA3", "#F29403", "#F781BF", "#BC9DCC", 
                  "#A65628", "#54B0E4", "#222F75", "#1B9E77", 
                  "#B2DF8A", "#E3BE00", "#FB9A99", "#E7298A", 
                  "#910241", "#00CDD1", "#A6CEE3", "#CE1261", 
                  "#5E4FA2", "#8CA77B", "#00441B", "#DEDC00", 
                  "#B3DE69", "#8DD3C7", "#999999")
  if (n <= length(colorSpace)) {
    colors <- colorSpace[1:n]
  }
  else {
    colors <- (grDevices::colorRampPalette(colorSpace))(n)
  }
  return(colors)
}

####some more divergent color
####from monet
monet.sunset <- colorRampPalette(colors = c("#272924","#323c48","#67879c",
                  "#7d9390","#d7695a","#ba9f84"))
monet.sunumbrella.women <- colorRampPalette(colors = c("#aed1e4","#bcdde2","#e8e7d5",
                                                       "#f1ede9","#e6cec2","#e49d40",
                                                       "#ddcb3e","#a39524","#89af7a",
                                                       "#744054","#aea3a9"))
monet.cliff <- colorRampPalette(colors = c("#242a26","#3c5653","#bacbc3","#3d4b6e",
                                           "#6f799d","#b3bedc","#f0f0f2","#898788",
                                           "#ddd0bd"))
####show pictures
pic.green <- colorRampPalette(colors = c("#344A23","#726342","#383812","#879979"))
pic.blue <- colorRampPalette(colors = c("#0B2C26","#004054","#95B7B6","#DFE5E1","#C26441"))
pic.green2 <- colorRampPalette(colors = c("#61776B","#A6AA8B","#B6C7D7","#EFEBEC","#E79A61"))
pic.beauty <- colorRampPalette(colors = c("#2F5365","#8A8E8F","#435D2E","#09200B","#CDBBA4","#D04A49"))
pic.beauty2 <- colorRampPalette(colors = c("#4B6032","#436A4D","#DFB28D","#BF4E31","#87988E"))
pic.beauty3 <- colorRampPalette(colors = c("#DB9371","#FCE6D8","#EB9B63","#DBD6C9","#666E20"))
pic.greenblue <- colorRampPalette(colors = c("#86AFA9","#E9A419","#FC8550","#F7C9BA","#8C5E51"))
pic.boy <- colorRampPalette(colors = c("#7C7F62","#956A3D","#625A58","#B1B1B1","#D8D4D5"))
pic.fellows <- colorRampPalette(colors = c("#193B14","#85A063","#8A4925","#C5947E","#ACA9A2"))
pic.girl <- colorRampPalette(colors = c("#1B3E64","#BE1705","#D77363","#84865D","#E2B88D"))
pic.girl2 <- colorRampPalette(colors = c("#BF593F","#DAAE7D","#AAC5DA","#DAC1BA","#141A22"))
molandi.color <- c("#e4d4c5","#c6b1ac","#764e56",
                   "#f9d9c0","#d19477","#93675a",
                   "#f0e9cd","#b9a783","#796656",
                   "#cdc1b1","#a2967e","#656356",
                   "#d8e7e4","#9eb2b1","#5a6873",
                   "#ccd8b0","#7e8563","#50463d")
molandi.color <- colorRampPalette(colors = molandi.color)
#scales::show_col(pic.girl2(10))

####-------5.misc function-----------

####-------5.1 ----------
MyMonocle2DiffGene <- function (cds, fullModelFormulaStr = "~CellType", 
                                reducedModelFormulaStr = "~1", 
                                relative_expr = TRUE, 
                                cores = 10, 
                                verbose = FALSE) 
{
  status <- NA
  if (class(cds)[1] != "CellDataSet") {
    stop("Error cds is not of type 'CellDataSet'")
  }
  all_vars <- c(all.vars(formula(fullModelFormulaStr)), 
                all.vars(formula(reducedModelFormulaStr)))
  pd <- pData(cds)
  for (i in all_vars) {
    x <- pd[, i]
    if (any((c(Inf, NaN, NA) %in% x))) {
      stop("Error: Inf, NaN, or NA values were located in pData of cds in columns mentioned in model terms")
    }
  }
  if (relative_expr && cds@expressionFamily@vfamily %in% c("negbinomial", 
                                                           "negbinomial.size")) {
    if (is.null(sizeFactors(cds)) || sum(is.na(sizeFactors(cds)))) {
      stop("Error: to call this function with relative_expr==TRUE, you must first call estimateSizeFactors() on the CellDataSet.")
    }
  }
  if (cores > 1) {
    diff_test_res <- mcesApply(cds, 1, 
                               monocle:::diff_test_helper, 
                               c("BiocGenerics", "VGAM", "Matrix"), 
                               cores = cores, 
                               fullModelFormulaStr = fullModelFormulaStr, 
                               reducedModelFormulaStr = reducedModelFormulaStr, 
                               expressionFamily = cds@expressionFamily, 
                               relative_expr = relative_expr, 
                               disp_func = cds@dispFitInfo[["blind"]]$disp_func, 
                               verbose = verbose)
    diff_test_res
  }
  else {
    diff_test_res <- smartEsApply(cds, 1, diff_test_helper, 
                                  convert_to_dense = TRUE, fullModelFormulaStr = fullModelFormulaStr, 
                                  reducedModelFormulaStr = reducedModelFormulaStr, 
                                  expressionFamily = cds@expressionFamily, relative_expr = relative_expr, 
                                  disp_func = cds@dispFitInfo[["blind"]]$disp_func, 
                                  verbose = verbose)
    diff_test_res
  }
  diff_test_res <- do.call(rbind.data.frame, diff_test_res)
  diff_test_res$qval <- 1
  diff_test_res$qval[which(diff_test_res$status == "OK")] <- p.adjust(subset(diff_test_res, 
                                                                             status == "OK")[, "pval"], method = "BH")
  diff_test_res <- merge(diff_test_res, fData(cds), by = "row.names")
  row.names(diff_test_res) <- diff_test_res[, 1]
  diff_test_res[, 1] <- NULL
  return(diff_test_res[row.names(cds), ])
}
####-------5.2 CellChat----------
####This function was used for CellChatCompare
myCompareInteractions <- function (object, measure = c("count", "weight"), 
                                   color.use = NULL, group = NULL, 
                                   group.levels = NULL, group.facet = NULL, 
                                   group.facet.levels = NULL, 
                                   n.row = 1, color.alpha = 1, legend.title = NULL, 
                                   width = 0.6, title.name = NULL, 
                                   digits = 3, xlabel = NULL, 
                                   ylabel = NULL, remove.xtick = FALSE, 
                                   show.legend = TRUE, 
                                   label.size = 3,
                                   x.lab.rot = FALSE, angle.x = 45, 
                                   vjust.x = NULL, hjust.x = 1, 
                                   size.text = 10) 
{
  measure <- match.arg(measure)
  if (measure == "count") {
    df <- as.data.frame(sapply(object@net, function(x) sum(x$count)))
    if (is.null(ylabel)) {
      ylabel = "Number of inferred interactions"
    }
  }
  else if (measure == "weight") {
    df <- as.data.frame(sapply(object@net, function(x) sum(x$weight)))
    df[, 1] <- round(df[, 1], digits)
    if (is.null(ylabel)) {
      ylabel = "Interaction strength"
    }
  }
  colnames(df) <- "count"
  df$dataset <- names(object@net)
  if (is.null(group)) {
    group <- 1
  }
  df$group <- group
  df$dataset <- factor(df$dataset, levels = names(object@net))
  if (is.null(group.levels)) {
    df$group <- factor(df$group)
  }
  else {
    df$group <- factor(df$group, levels = group.levels)
  }
  if (is.null(color.use)) {
    color.use <- ggPalette(length(unique(group)))
  }
  if (!is.null(group.facet)) {
    if (all(group.facet %in% colnames(df))) {
      gg <- ggplot(df, aes(x = dataset, y = count, fill = group)) + 
        geom_bar(stat = "identity", width = width, 
                 position = position_dodge())
      gg <- gg + facet_wrap(group.facet, nrow = n.row)
    }
    else {
      df$group.facet <- group.facet
      if (is.null(group.facet.levels)) {
        df$group.facet <- factor(df$group.facet)
      }
      else {
        df$group.facet <- factor(df$group.facet, levels = group.facet.levels)
      }
      gg <- ggplot(df, aes(x = dataset, y = count, fill = group)) + 
        geom_bar(stat = "identity", width = width, 
                 position = position_dodge())
      gg <- gg + facet_wrap(~group.facet, nrow = n.row)
    }
  }
  else {
    gg <- ggplot(df, aes(x = dataset, y = count, fill = group)) + 
      geom_bar(stat = "identity", width = width, 
               position = position_dodge())
  }
  gg <- gg + geom_text(aes(label = count), 
                       vjust = -0.3, 
                       size = label.size, 
                       position = position_dodge(0.9))
  gg <- gg + ylab(ylabel) + xlab(xlabel) + theme_classic() + 
    labs(title = title.name) + theme(plot.title = element_text(size = 10, 
                                                               face = "bold", hjust = 0.5)) + theme(text = element_text(size = size.text), 
                                                                                                    axis.text = element_text(colour = "black"))
  gg <- gg + scale_fill_manual(values = alpha(color.use, alpha = color.alpha), 
                               drop = FALSE)
  if (remove.xtick) {
    gg <- gg + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  }
  if (is.null(legend.title)) {
    gg <- gg + theme(legend.title = element_blank())
  }
  else {
    gg <- gg + guides(fill = guide_legend(legend.title))
  }
  if (!show.legend) {
    gg <- gg + theme(legend.position = "none")
  }
  if (x.lab.rot) {
    gg <- gg + theme(axis.text.x = element_text(angle = angle.x, 
                                                hjust = hjust.x, vjust = vjust.x, size = size.text))
  }
  gg
  return(gg)
}


myrankSimilarity <- function (object=cellchat, slot.name = "netP",
                              type = c("functional", "structural"), 
                              comparison1 = NULL, 
                              comparison2 = c(1, 2), x.rotation = 90, title = NULL, 
                              color.use = NULL, bar.w = NULL, 
                              font.size = 8) 
{
  type <- match.arg(type)
  if (is.null(comparison1)) {
    comparison1 <- 1:length(unique(object@meta$datasets))
  }
  comparison.name <- paste(comparison1, collapse = "-")
  cat("Compute the distance of signaling networks between datasets", 
      as.character(comparison1[comparison2]), "\n")
  comparison2.name <- names(methods::slot(object, slot.name))[comparison1[comparison2]]
  Y <- methods::slot(object, slot.name)$similarity[[type]]$dr[[comparison.name]]
  group <- sub(".*--", "", rownames(Y))
  data1 <- Y[group %in% comparison2.name[1], ]
  data2 <- Y[group %in% comparison2.name[2], ]
  rownames(data1) <- sub("--.*", "", rownames(data1))
  rownames(data2) <- sub("--.*", "", rownames(data2))
  pathway.show = as.character(intersect(rownames(data1), rownames(data2)))
  data1 <- data1[pathway.show, ]
  data2 <- data2[pathway.show, ]
  euc.dist <- function(x1, x2) sqrt(sum((x1 - x2)^2))
  dist <- NULL
  for (i in 1:nrow(data1)) dist[i] <- euc.dist(data1[i, ], data2[i, ])
  df <- data.frame(name = pathway.show, dist = dist, row.names = pathway.show)
  df <- df[order(df$dist), , drop = F]
  df$name <- factor(df$name, levels = as.character(df$name))
  gg <- ggplot(df, aes(x = name, y = dist)) + 
    geom_bar(stat = "identity",width = bar.w,fill = color.use) + 
    theme_classic() + 
    theme(text = element_text(size = font.size), 
          axis.title.y = element_text(size = font.size)) + 
    xlab("") + 
    ylab("Pathway distance") + 
    coord_flip()
  if (!is.null(title)) {
    gg <- gg + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  }
  return(gg)
}

