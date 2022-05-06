BiocManager::install("edgeR")
library(edgeR)
#####载入数据#####
path = "~/case/Wilcoxon.test/data/"#此处可替换自己的counts数据
fileName = dir(path)
mydata = read.table(file = paste(path,fileName[1],sep=""),
                       header = T,row.names = 1,stringsAsFactors = F,check.names = F)
for(k in 1:length(fileName)){
  mydata[k] = read.table(file = paste(path,fileName[k],sep=""),
                         header = T,row.names = 1,stringsAsFactors = F,check.names = F)
  names(mydata)[k] <- fileName[k]
}
#####准备样品信息矩阵#####
conditions<-read.table(file="condition.txt", header = F)#"condition.txt"见Wilcoxon.test文件夹
conditions<-factor(t(conditions))
#####使用 edgeR 包进行计数矩阵预处理#####
y <- DGEList(counts=mydata,group=conditions)
#过滤，删除低表达基因
keep <- filterByExpr(y)
y <- y[keep,keep.lib.sizes=FALSE]
#TMM 归一化并转化为CPM（每百万计数）
y <- calcNormFactors(y,method="TMM")
count_norm=cpm(y)#因为每个基因长度不会改变，因此CPM就够了
count_norm<-as.data.frame(count_norm)
#####对每个基因进行Wilcoxon秩和检验#####
pvalues <- sapply(1:nrow(count_norm),function(i){
  data<-cbind.data.frame(gene=as.numeric(t(count_norm[i,])),conditions)
  p=wilcox.test(gene~conditions, data)$p.value#wilcox检验适用于成对的数据且无法保证其正态分布
  return(p)
})
#BH校正后的p值
fdr=p.adjust(pvalues,method = "fdr")
#####计算每个基因的fold-change#####
conditionsLevel<-levels(conditions)
dataCon1=count_norm[,c(which(conditions==conditionsLevel[1]))]
dataCon2=count_norm[,c(which(conditions==conditionsLevel[2]))]
foldChanges=log2(rowMeans(dataCon2)/rowMeans(dataCon1))
#####基于p-adjust的输出结果#####
outRst<-data.frame(log2foldChange=foldChanges, pValues=pvalues, padj=fdr)
rownames(outRst)=rownames(count_norm)
fdrThres=0.05
outRst_filter <- outRst[outRst$padj<fdrThres,]
#the code for volcano plot火山图
library(ggplot2)
library(dplyr)
library(readxl)

plotdata <- outRst

log2foldchange_thres = 1
log10pvalue_thres = -log10(0.05)
plotdata <- mutate(plotdata, log10_pvalue = -log10(padj))
plotdata <- mutate(plotdata, type = ifelse(abs(log2foldChange)>=log2foldchange_thres&log10_pvalue>log10pvalue_thres,
                                   'change', 'nochange'))

ggplot(data = plotdata, aes(x = log2foldChange, y = log10_pvalue)) + 
  geom_point(aes(color = type)) + 
  scale_color_manual(values = c('red', 'gray')) + 
  geom_hline(yintercept = log10pvalue_thres) + 
  geom_vline(xintercept = c(log2foldchange_thres, -log2foldchange_thres)) +
  xlab(label = "log2 fold change") + 
  ylab(label = '-log10 adjust p-value') +
  theme(legend.position = 'none',
        panel.background = element_blank(),
        panel.border = element_rect(color = 'black', fill = NA))

