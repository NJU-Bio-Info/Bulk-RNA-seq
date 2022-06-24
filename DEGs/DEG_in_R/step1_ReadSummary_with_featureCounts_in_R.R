library(Rsubread)
library(clusterProfiler)
library(dplyr)

#bam file
bam <- list.files(path = 'bam', pattern = 'bam$')
bam

files <- file.path('bam', bam)
files

#use featureCounts() to perform reads summary and store results in a list
counts <- featureCounts(files = files,
                        annot.inbuilt = 'mm39',#use in-built mm39 genome
                        GTF.featureType = 'exon',#only count the reads from exons
                        GTF.attrType = 'gene_id',#perform meta-feature reads summary at gene level
                        useMetaFeatures = TRUE,#perform reads summary at meta-feature level
                        strandSpecific = 0,#non-strand specific RNA-seq
                        isPairedEnd = TRUE,#indicate pair end sequencing
                        countReadPairs = TRUE,#count read pairs rather than reads
                        requireBothEndsMapped = TRUE,#only count the pairs with both reads mapped
                        nthreads = 4)#use 4 threads to perform reads summary
names(counts)
#counts: the reads count of each gene in each sample.
View(counts$counts)
#annotation: detailed information about every gene include exon coordinates, effective length etc.
View(counts$annotation)
#targets: the bam files used for reads summary.
counts$targets
#stat: reads summary report for every sample.
counts$stat %>% head()

#do not need run
if(F){
  saveRDS(object = counts$counts, file = 'featureCounts.rds')
}