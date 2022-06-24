invisible(gc())

library(limma)
library(Mfuzz)
library(dplyr)


#step1: remove the effect of library size
#the counts matrix from step2
#if you can not run step2:
if(F){
  new_counts <- readRDS(file = 'Counts_without_batch.rds')
}
normalize_factor <- colSums(new_counts)
normalize_counts <- t(t(new_counts)/normalize_factor)*1e6
#step2: construct mfuzz object
data <- t(normalize_counts)
rownames(data) <- c(rep('0h', time = 6), rep(c('20h', '9h'), each = 3))
data <- avereps(data)

#filter the genes with 0 values across 3 time points
keep <- colSums(data)!=0
data <- data[ , keep] %>% t()
data <- data[, c(1,3,2)]

eset <- new('ExpressionSet', exprs = data)
#make genes comparable
eset <- standardise(eset = eset)

m <- mestimate(eset = eset)
#centers = 6, 6 clusters
set.seed(2022)
cl <- mfuzz(eset = eset, centers = 6, m = m)

pdf('Figure/Fig3.time_series_cluster.pdf')
mfuzz.plot2(eset = eset,
            cl = cl,
            mfrow = c(3, 2),
            x11 = F,
            time.labels = colnames(data))
dev.off()

#you can also extract the genes in each cluster
#the cluster1 as an example
which(cl$cluster == 1) %>% names()