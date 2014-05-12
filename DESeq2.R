source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")

library("DESeq2")

countsTable <- read.delim("rsem.shrimp.counts.txt",header=TRUE)
countsTable <- read.delim("express_shrimp_D.counts.txt",header=TRUE)

colData <- read.delim("metaData.rsem.txt", header=FALSE)
colData <- read.delim("metaData.eXpress.txt", header=FALSE)

rownames(countsTable) <- countsTable$gene_id
rownames(countsTable) <- countsTable$target_id

countsTable <- countsTable[,-1]
rownames(countsTable)
head(countsTable$target_id)
rownames(colData) <- colData$V1
colData$V1 <- colData$V2
colData2 <- colData[1]
colnames(colData2) <- c("treatment")

ddsFullCountTable <- DESeqDataSetFromMatrix(
  countData = countsTable,
  colData = colData2,
  design = ~ treatment)
ddsFullCountTable

DGE <- DESeq(ddsFullCountTable, betaPrior = FALSE)
DGE
results<- results(DGE, contrast = c("treatment", "H4", "N4"))
results
resSig <- results[ which(results$padj < 0.1 ), ]

resSig
write.table(results, file = "H4vN4.deseq.rsem.results.txt", sep = "\t")

mcols(results, use.names=TRUE)
results

