library(DESeq2)

rawCounts <- read.table("edgeR_inpit.txt", sep="\t", header=T, stringsAsFactor=F)

A <- rawCounts$GENE
B <-rawCounts[,2:453]
C <- cbind(A,B)

counts <- as.matrix(C[,-1])

rownames(counts) <- C[,1]
head(counts[1:3,1:3])

T=c()
for (i in colnames(B)[1:226]){
  s = unlist(strsplit(i,split=".T.genes.results"))
  T = c(T,s)
}
N=c()
for (i in colnames(B)[1:226]){
  s = unlist(strsplit(i,split=".T.genes.results"))
  N = c(N,s)
}

Sample=factor(c(T,N))
Group=factor(rep(c("T","N"),c(226,226)))
Group <- relevel(Group, ref = "N")
sampleData <- data.frame(Sample=Sample,Group)

dds <- DESeqDataSetFromMatrix(countData = round(counts), colData = sampleData, design = ~ Sample + Group)
dds <- DESeq(dds)

dds_results <- results(dds, contrast = c("Group","T", "N"))
dds_results_test<- results(dds)
dds_results <- dds_results[order(dds_results$padj),]
write.csv(dds_results, file="DESeq2_results.csv")
############################################################

par(mfrow = c(1,1))
plotMA(dds_results)

# Top 25 differentially expressed genes
par(mfrow = c(5,5))
for(i in 1:25){
  plotCounts(dds, gene = rownames(dds_results)[i], intgroup = "Group", pch = 19)
}

# Volcano plot
par(mfrow = c(1,1))
with(dds_results, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot"))

# Add colored points: blue if padj<0.01, red if log2FC>2 and padj<0.01)
with(subset(dds_results, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch=20, col="blue"))
with(subset(dds_results, padj<.01 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(pvalue), pch=20, col="red"))

# Get differentially expressed gene matrix
dds_significant <- dds_results[!is.na(dds_results$padj) &
                                 dds_results$padj<0.10 &
                                 abs(dds_results$log2FoldChange)>=1,]

# Sort in Log2FoldChange order
dds_significant<-dds_significant[order(dds_significant$log2FoldChange, decreasing = T),]
significant_genes_ID<-rownames(dds_significant)
# Convert ensembl ID to gene symbol
for(i in 1:nrow(dds_significant)){
  rownames(dds_significant)[i]<-genes$Gene.Name[which(genes$Gene.ID==rownames(dds_significant)[i])]
}

significant_genes_symbol<-rownames(dds_significant)
write.csv(dds_significant, file="DESeq2_significant_results.csv")

# Making heatmap
install.packages("RColorBrewer")
library(RColorBrewer)

BiocManager::install("gplots")
library(gplots)

# colors of the heat map
hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)

# Heatmap
significant_genes_count<-log2(counts(dds, normalized=T)[significant_genes_ID,] + 1)
rownames(significant_genes_count)<-significant_genes_symbol
heatmap.2(significant_genes_count,
          col = hmcol, scale="row", Rowv = TRUE, Colv = FALSE,
          dendrogram="row", trace="none", margin=c(4,6), cexRow=0.5, cexCol=1, keysize=1)

