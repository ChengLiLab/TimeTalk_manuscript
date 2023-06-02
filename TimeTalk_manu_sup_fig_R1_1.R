#### Try to address the review 1's suggestion
library(kBET)
library(tidyverse)
library(ggplot2)
library(cowplot)
source("code/myUtils.R")
####------0.test data on simulation data------------

testdata <- create_testset_multibatch(n.genes=1000,
                                      n.batch=2, plattform='any')
pca.data <- prcomp(testdata$data, center=TRUE)
batch.silhouette.sim <- batch_sil(pca.data, testdata$batch)
batch.estimate.sim <- kBET(df = t(testdata$data),batch = testdata$batch,)

batch.estimate.sim$outsider
batch.estimate.sim$average.pval
batch.estimate.sim$summary


####------1.load data and kBET analysis-----------
beforeEPI.exp <- readRDS(file = "res/R/beforeEPI_exp_20210114.rds")
beforeEPI.meta <- readRDS(file = "res/R/beforeEPI_meta_20210114.rds")

data.plot <- beforeEPI.exp
tmp.select <- beforeEPI.meta %>%
  dplyr::filter(Stage!="fibroblast") %>%
  pull(cell_id)
data.plot <- data.plot[,tmp.select]
batch <- c(rep(1,3),rep(2,ncol(data.plot)-3))

table(batch)
batch.estimate <- kBET(data.plot,batch = batch)

batch.estimate$stats

batch.estimate$summary



batch.estimate <- kBET(data.plot, batch, plot=FALSE)
plot.data <- data.frame(class=rep(c('observed', 'expected'), 
                                  each=length(batch.estimate$stats$kBET.observed)), 
                        data =  c(batch.estimate$stats$kBET.observed,
                                  batch.estimate$stats$kBET.expected))
g <- ggplot(plot.data, aes(class, data)) + geom_boxplot() + 
  labs(x='Test', 
       y='Rejection rate',
       title='kBET test results') +
  theme_cowplot(font_size = 22) +  
  scale_y_continuous(limits=c(0,1))+
  theme(axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1),
        plot.title = element_text(hjust = 0.5))
g
batch.estimate$summary

ggsave(filename = myFileName(prefix = "res/fig/figS_revision_R1_1_kBET_analysis",
                             suffix = ".jpg"),
       width = 6,height = 6,
       dpi = 350)
ggsave(filename = myFileName(prefix = "res/fig/figS_revision_R1_1_kBET_analysis",
                             suffix = ".pdf"),
       width = 6,height = 6,
       dpi = 350)





###--------2. silhouette width and PCA based measure Prepare data-----

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


batch <- c(rep(1,3),rep(2,ncol(data.plot)-3))
batch.silhouette <- batch_sil(pca.res, batch)

batch.silhouette





data <- beforeEPI.exp
table(batch)

pca.data

batch.estimate <- kBET(data,batch = batch)


batch.estimate$results
batch.estimate$average.pval
batch.estimate$summary
?kBET

batch <- rep(seq_len(10),each=20)
data <- matrix(rpois(n = 50000, lambda = 10)*rbinom(50000,1,prob=0.5), nrow=200)

batch.estimate <- kBET(data,batch)
batch.estimate$summary
?kBET


