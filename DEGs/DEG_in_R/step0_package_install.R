############## packages install ##############
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")
}

if(!require("Rsubread", quietly = T)) BiocManager::install("Rsubread")
if(!require("dplyr", quietly = T)) install.packages('dplyr')
if(!require('ggplot2', quietly = T)) install.packages('ggplot2')
if(!require('clusterProfiler', quietly = T)) BiocManager::install('clusterProfiler')
if(!require('Mfuzz', quietly = T)) BiocManager::install('Mfuzz')
if(!require('limma', quietly = T)) BiocManager::install('limma')
if(!require('DEseq2', quietly = T)) BiocManager::install("DESeq2")
if(!require('sva', quietly = T)) devtools::install_github('zhangyuqing/sva-devel')