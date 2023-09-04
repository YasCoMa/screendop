library(DESeq2) 
data=as.matrix(read.csv('immunoreactive_count_matrix.txt', sep=',', row.names='id', header=TRUE)) 
design=read.csv('immunoreactive_design_matrix.txt', sep=',', row.names=1, header=TRUE) 
data[is.na(data)] <- 0
dataset <- DESeqDataSetFromMatrix(countData = data, colData = design, design = ~treatment) 
dds <- DESeq(dataset) 
res <- results(dds) 
data_ = as.data.frame(res) 
data_ = data_[!is.na(data_$log2FoldChange),] 
data_ = data_[!is.na(data_$padj),] 
write.table(data_, file='immunoreactive_result.csv', sep=',', quote=FALSE) 
