source("http://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)

raw.data <- read.table( file = "rsem.shrimp.counts.txt" , header = TRUE )
counts <- raw.data[,-1]

#for express
colnames(counts) <- colnames(raw.data)[-1]
rownames(counts) <- raw.data$target_id
#for rsem
rownames(counts) <- raw.data$gene_id

head(counts)
#RSEM counts grouping
group <- c(rep("H24", 6) , rep("H4", 6), rep("HH24", 6) , rep("HH4", 6), rep("N24", 6) , rep("N4", 6))
#Express counts grouping
group <- c(rep("N4", 6) , rep("N24", 6), rep("H4", 6) , rep("H24", 6), rep("HH4", 6) , rep("HH24", 6))

cds <- DGEList( counts , group = group )
cds <- cds[rowSums(1e+06 * cds$counts/expandAsMatrix(cds$samples$lib.size, dim(cds)) > 1) >= 3, ]
cds <- calcNormFactors( cds )
cds <- estimateCommonDisp( cds )
cds <- estimateTagwiseDisp( cds )

de.N4vH4 <- exactTest( cds,  pair = c( "N4" , "H4" ) )
de.N4vHH4 <- exactTest( cds,  pair = c( "N4" , "HH4" ) )
de.N24vH24 <- exactTest( cds,  pair = c( "N24" , "H24" ) )
de.N24vHH24 <- exactTest( cds,  pair = c( "N24" , "HH24" ) )



resultsTbl <- topTags( de.N24vHH24 , n = nrow(de.N24vHH24$table) )$table
de.genes <- rownames(resultsTbl )[ resultsTbl$FDR <= 0.1 ]
length(de.genes)
write.table(resultsTbl, "N24vHH24.edgeR.rsem.txt")
write(de.genes, "genes.N24vHH24.edgeR.rsem.txt")



