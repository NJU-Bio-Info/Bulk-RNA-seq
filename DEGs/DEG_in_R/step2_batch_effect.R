invisible(gc())

if(!dir.exists('Figure')) dir.create('Figure')

library(sva)
library(DESeq2)
library(dplyr)
library(ggplot2)

#the reads summary results from step1
cts <- counts$counts

#if you can not run featureCounts
#cts <- readRDS(file = 'featureCounts.rds')

colnames(cts)

#check the batch effect
coldata <- data.frame(condition = c(rep(paste0('0h_batch', c(1, 2)), time = 3), rep(c('20h', '9h'), each = 3)))
rownames(coldata) <- colnames(cts)

dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~condition)
vsd <- vst(object = dds, blind = F)

plotPCA(object = vsd, intgroup = 'condition') + 
  ggtitle(label = 'PCA plot for the data without removing batch effect')
#because the plotPCA() function use ggplot2 to produce figure, we can use ggsave() to save figure
if(F){
  ggsave(filename = 'Figure/Fig1.PCA_with_batch.pdf')
}


#here we can find that there is batch effect from the PCA results (oh).
#use combat_seq() (from sva package) to remove batch effect
#you must use raw counts rather than RPKM/TPM to perform such task
batch <- c(rep(c('batch1', 'batch2'), time = 3), rep('batch1', time = 6))
group <- c(rep('0h', time = 6), rep(c('20h', '9h'), each = 3))
#the matrix here: row~gene; col~sample
new_counts <- ComBat_seq(counts = cts,#raw counts table
                         batch = batch,#batch vector
                         group = group,#biological condition of interest
                         full_mod = TRUE#include condition of interest in model
                         )
#check batch effect again
new_coldata <- data.frame(condition = c(rep(paste0('0h_batch', c(1, 2)), time = 3), rep(c('20h', '9h'), each = 3)))
rownames(new_coldata) <- colnames(new_counts)

new_dds <- DESeqDataSetFromMatrix(countData = new_counts,
                              colData = new_coldata,
                              design = ~condition)
new_vsd <- vst(object = new_dds, blind = F)
plotPCA(object = new_vsd, intgroup = 'condition') + 
  ggtitle(label = 'PCA plot for the data with removing batch effect')
#because the plotPCA() function use ggplot2 to produce figure, we can use ggsave() to save figure
if(F){
  ggsave(filename = 'Figure/Fig2.PCA_without_batch.pdf')
}

#save counts data, you do not need to run
saveRDS(object = new_counts, file = 'Counts_without_batch.rds')