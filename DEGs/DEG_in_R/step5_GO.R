invisible(gc())

library(clusterProfiler)
library(ggplot2)
library(org.Mm.eg.db)

#the DESeq2 result from step4
#here we take the genes downregulated in 9h sample as an example
#extract the genes

#cut off for foldchange
foldchange = 5
#cut off for adjust p value
p.adj = 1e-3
genes <- res_deseq %>% 
  filter(padj < p.adj & log2FoldChange > log2(foldchange)) %>% 
  rownames()
length(genes)

#the gene id type
keytypes(org.Mm.eg.db)

GO_result <- enrichGO(gene = genes, #genes vector
                      OrgDb = 'org.Mm.eg.db', #database, org.Mm.eg.db for mouse; org.Hs.eg.db for human
                      keyType = 'ENTREZID', #gene id type
                      ont = 'ALL',
                      readable = T
                      )
dotplot(GO_result)
ggsave(filename = 'Figure/Fig5.goenrich.pdf')
