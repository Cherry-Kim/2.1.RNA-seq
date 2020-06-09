##Install cummeRbund #R version 3.6.3 (2020-02-29)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
BiocManager::install("cummeRbund")
library(cummeRbund)

#Reading cuffdiff output
cuff <- readCufflinks( )

#Determining relationships between conditions
csBoxplot<-csBoxplot(genes(cuff))
pdf("csBoxplot.pdf")
csBoxplot
dev.off()

dens<-csDensity(genes(cuff))
pdf("dens.pdf")
dens
dev.off()

pdf(file="scatterplot.pdf")
csScatter(genes(cuff), 'group1', 'group2',smooth=T)
dev.off()

pdf(file="volcanoplot.pdf")
csVolcano(genes(cuff), 'group1', 'group2', alpha= 0.05, showSignificant=T)
dev.off()

gene.matrix<-fpkmMatrix(genes(cuff))
gene.count.matrix<-countMatrix(genes(cuff))

##Analysis of differential expression
#Identifying differentially expressed genes
myGeneIds<-getSig(cuff,alpha=0.05,level="genes")
#head(sigGeneIds)
#length(sigGeneIds)

#Creating significant gene sets
myGenes<-getGenes(cuff,myGeneIds)
#head(fpkm(myGenes))

#Geneset level plots
csHeatmap<-csHeatmap(myGenes,cluster='both')
pdf("csHeatmap.pdf")
csHeatmap
dev.off()

expressionBarplot<-expressionBarplot(myGenes)
pdf("expressionBarplot.pdf")
expressionBarplot
dev.off()

##Individual Genes
myGeneId<-"Actn3"
#myGeneId<-"XLOC_019880"
myGene <- getGene(cuff,myGeneId)
#head(fpkm(myGene))

#Gene-level plots
gl<-expressionPlot(myGene)
pdf(file="expressionPlot.pdf")
gl
dev.off()

gb<-expressionBarplot(myGene)
pdf(file="expressionBarplot.pdf")
gb
dev.off()

gp<-csPie(myGene,level="isoforms")
pdf(file="csPie.pdf")
gp
dev.off()
