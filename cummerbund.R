source("http://bioconductor.org/biocLite.R")
biocLite("cummeRbund")
library(cummeRbund)
cuff <- readCufflinks(replicates =T)
cuff
genes <- genes(cuff)
ann <- annotation(genes)
write.table(featureNames(genes), "gene_ids.txt", sep="\t")
write.table(ann$locus, "locus.txt", sep="\t")

write.table(table(featureNames(genes), ann$locus), "ids_to_locus.txt", sep="\t")

featureNames(genes)
dendro <- csDendro(genes(cuff), replicates=T)
geneIDs <- getSig(cuff, alpha = 1)
length(geneIDs)

sigID <- getSig(cuff, 'H24', 'HH24', alpha = .05, level='genes')
diffGenes<-getGenes(cuff,sigID)
length(diffGenes)
ann <- annotation(diffGenes)
write(ann$locus, "H24vHH24.cuff.txt", sep="\t")
      

length(N4vH4.sigID)
N4vH4.sigID
sigIDs <- getSig(cuff, alpha = .05)
length(sigIDs)
?getSig
myGenes<-getGenes(cuff, geneIDs)
map <- csHeatmap(myGenes, cluster='both', replicates = T)
