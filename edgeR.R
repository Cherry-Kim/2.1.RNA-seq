#!/usr/bin/env Rscript

################################################
### STEP1. Construct input file to edgeR
###############################################

BiocManager::install("tximportData")
library(tximport)
path.case = "/home/hykim/RNA-seq/Quantified_Gene_Expression/"

files_N <- list.files(path = path.case, pattern = "N.genes.results$", full.names = TRUE)
files_T <- list.files(path = path.case, pattern = "T.genes.results$", full.names = TRUE)
files <- c(files_N, files_T)
print(length(files))

N.list<-dir(path=path.case, pattern = "N.genes.results")
sample.N.list <- list()
for( i in 1:length(N.list)){
	sample.N.list[i] <- unlist(strsplit(N.list[i],split=".genes.results"))[1]
	}	
T.list<-dir(path=path.case, pattern = "T.genes.results")
sample.T.list <- list()
for( i in 1:length(T.list)){
	sample.T.list[i] <- unlist(strsplit(T.list[i],split=".genes.results"))[1]
	}	
sample.list <- c(sample.N.list, sample.T.list)
names(files) <- sample.list

txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
a <- round(txi.rsem$counts[,], digit=0)
write.table(a, "edgeR_input.txt",col.names=NA, row.names=T, quote=F,sep='\t')
##write.table(txi.rsem$abundance, "edgeR_input.txt",col.names=NA, row.names=T, quote=F,sep='\t')

####################################################
### STEP2. edgeR
###################################################
library(edgeR)

rawdata <- read.delim("edgeR_input.txt", check.names=FALSE, stringsAsFactors=FALSE)
group <- rep(c("N","T"), each=226)
##group <- rep(c("N","T"),c(226,226))
y <- DGEList(counts=rawdata[,2:453], genes=rawdata[,1], group=group)

o <- order(rowSums(y$counts), decreasing=TRUE)
y <- y[o,]

d <- duplicated(y$genes$genes)
print(d)
y <- y[!d,]


y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)

png("MDS.png")
plotMDS(y)
dev.off()

y <- estimateDisp(y)

Group = factor(rep(c("N","T"),c(226,226)))
data.frame(Sample=colnames(y),Group)

design <- model.matrix(~Group)
rownames(design) <- colnames(y)

y <- estimateDisp(y, design, robust=FALSE)

et <- exactTest(y, pair=c("N","T"))
et_fdr <- topTags(et, adjust.method = "BH", sort.by = "PValue", p.value = 1, n=nrow(et$table))
write.table(et_fdr, "results_edgeR.txt", col.names=NA, row.names=T, quote=F,sep='\t')

a <- et_fdr$table
b <- which(a$FDR<=0.05)
sig.b <- a[b,]
write.table(sig.b, "results_edgeR_sig.txt", col.names=NA, row.names=T, quote=F,sep='\t')



