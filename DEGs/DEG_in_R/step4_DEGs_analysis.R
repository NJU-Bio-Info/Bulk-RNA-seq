invisible(gc())

library(dplyr)
library(ggplot2)
library(DESeq2)

#################################### function to produce volcano plot #########################################
plotVol <- function(data, foldchange, p.adj, colors = c('red', 'blue', 'gray')){
  plotData <- data %>% 
    mutate(type = case_when(
      padj < p.adj & log2FoldChange > log2(foldchange) ~ 'up',
      padj < p.adj & log2FoldChange < -log2(foldchange) ~ 'down',
      TRUE ~ 'none'
    ))
  plotData$type = factor(plotData$type, levels = c('up', 'down', 'none'))
  title = paste0('up-regulated:', nrow(filter(plotData, type == 'up')), '\n',
                 'down-regulated:', nrow(filter(plotData, type == 'down')), '\n',
                 'Cut-off for foldchange:', foldchange, '\n',
                 'Cut-off for pvalue:', p.adj)
  ggplot(data = plotData, aes(x = log2FoldChange, y = -log10(padj))) + 
    geom_point(aes(color = type, size = log2(baseMean)), alpha = .5) +
    scale_color_manual(values = colors) + 
    geom_hline(yintercept = -log10(p.adj), color = 'gray', lty = 4, lwd = .8) + 
    geom_vline(xintercept = c(log2(foldchange), -log2(foldchange)), color = 'gray', lty = 4, lwd = .8) + 
    ggtitle(label = title) + 
    xlab(label = 'log2FC(0h/9h)') +
    ylab(label = '-log10(adjust p value)') + 
    theme_test() + 
    theme(plot.title = element_text(hjust = 0.5), 
          legend.position = "right",
          axis.title = element_text(color = 'black', size = 15),
          axis.text = element_text(color = 'black', size = 10))
}

#here, we will try to identify the differentially expressed genes between 0h and 9h
#################################### With DESeq2 ########################################
#the input data is raw counts matrix, row ~ gene; col ~ sample
#the count matrix from step2
#if you can not run step2:
if(F){
  new_counts <- readRDS(file = 'Counts_without_batch.rds')
}
#matrix
ctx <- new_counts[ , c(1:6, 10:12)]
#condition
coldata <- data.frame(condition = c(rep('0h', time = 6), rep('9h', time = 3)))
rownames(coldata) <- colnames(ctx)
coldata
#notice that the columns of the matrix and the rows of the coldata are in the same order.
dds <- DESeqDataSetFromMatrix(countData = ctx,#count matrix
                              colData = coldata,#sample information
                              design = ~condition#biological factor of interest
                              )
dds #40838 genes, 9 samples
#filter lowly expressed genes
keep <- rowSums(counts(dds)) >= 9
dds <- dds[keep, ] 
dds #11524 genes, 9 samples

dds <- DESeq(dds)
#save as dataframe
res_deseq <- results(object = dds, contrast = c('condition', '0h', '9h')) %>% na.omit() %>% as.data.frame()

#cut off for foldchange
foldchange = 5
#cut off for adjust p value
p.adj = 1e-3

plotVol(data = res_deseq,#DESeq2 result
        foldchange = foldchange,#cut off for foldchange
        p.adj = p.adj,#cut off for adjust p value
        colors = c('red', 'blue', 'gray')#specify the colors, red ~ up-regulated, blue ~ down-regulated, gray ~ none
        )
ggsave(filename = 'Figure/Fig4.volcano.pdf')
