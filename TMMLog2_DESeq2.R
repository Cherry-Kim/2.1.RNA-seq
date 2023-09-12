#!/usr/bin/env Rscript
STEP1_tximport <- function(dir
           ){
        library(tximport)
        samples <- read.table(file.path('/BIO/', "mapfile.txt"), header = TRUE)
        head(samples,2)
#         sample  TREATMENT
# s1    N untreated
# s2    T   treated

        files_C <- paste0(  samples[which(samples$TREATMENT == "untreated"),1], ".W.genes.results")
        sample.C.list <- list()
        for( i in 1:length(files_C)){
                sample.C.list[i] <- unlist(strsplit(files_C[i], split=".genes.results"))[1]
        }

        files_T <- paste0(  samples[which(samples$TREATMENT == "treated"),1], ".B.genes.results")
        sample.T.list <- list()
        for( i in 1:length(files_T)){
                sample.T.list[i] <- unlist(strsplit(files_T[i], split=".genes.results"))[1]
        }

        files <- c(files_T, files_C)

        sample.list <- c(sample.T.list, sample.C.list)
        names(files) <- sample.list

        txi.rsem <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE)
        write.table(txi.rsem$counts, "txi.rsem.counts.txt", col.names=NA, row.names=T, quote=F,sep='\t')
        write.table(txi.rsem$abundance, "txi.rsem.abundance.TPM.txt", col.names=NA, row.names=T, quote=F,sep='\t')
}

STEP2_Gene <- function(){
        library('biomaRt')
        df <- read.table("txi.rsem.counts.txt", sep="\t", header=T, stringsAsFactor=F)
        mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")) #mmusculus_gene_ensembl

        genes <- df$X
        G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","entrezgene_id"),values=genes,mart= mart)

        df$hgnc_symbol = ""
        df["hgnc_symbol"] = lapply("hgnc_symbol", function(x) G_list[[x]][match(df$X, G_list$ensembl_gene_id)])
        df$entrezgene_id=""
        df["entrezgene_id"] = lapply("entrezgene_id", function(x) G_list[[x]][match(df$X, G_list$ensembl_gene_id)])
        df2 <- df[,c(1,58,59,2:57)]
        write.table(df2, "txi.rsem.counts.gene.txt", col.names=NA, row.names=T, quote=F,sep='\t')
}

STEP3_TMMnorm <- function(){
        library(edgeR)
        Data <- read.table("txi.rsem.counts.gene.txt", header=T, sep='\t')
        head(Data,2)
# X.1               X hgnc_symbol entrezgene_id s1 s2
#1   1 ENSG00000000003      TSPAN6          7105              74              34

        CountData <- Data[,seq(from=5,to=ncol(Data))]
        rownames(CountData) <- Data[,2]
        a <- CountData[,1:36]
        b <- CountData[,37:56]
        CountData2 <- cbind(b,a)

        Group=rep(c("control","cancer"),c(20,36))
        Group <-factor(Group,levels=c("control","cancer")) 

        y <- DGEList(counts=CountData2, group=Group)
        y <- calcNormFactors(y, method="TMM")
        #plotMDS(y, col = as.numeric(Group))
        TMM <- cpm(y, normalized.lib.sizes=TRUE)        #log=T

        #2. log2+1
        dat.log2 <- TMM
        for (i in 1:ncol(TMM)){
                dat.log2[,i] <- log(TMM[,i]+1,2)
        }
       #=========================================================
        itr <- c()
        for(i in 1:nrow(dat.log2)){
                if(sd(dat.log2[i,])==0)
                        {itr <- c(itr,i)}
        }
        length(itr)
        in.dat <- dat.log2[-itr,]
        #=========================================================

        TMM_ColName <- paste("TMMlog2_",gsub(".bam", "",colnames(dat.log2)),sep="")
        colnames(dat.log2) <- TMM_ColName
        
        ### DEG Analysis - 1. DESeq2 ###
        library(DESeq2)
        W=c()
        for (i in colnames(dat.log2)[1:20]){
                W = c(W,i)
        }
        B=c()
        for (i in colnames(dat.log2)[21:56]){
          B = c(B,i)
        }
        Sample=factor(c(W,B))
        Group <- relevel(Group, ref = "control")
        sampleData <- data.frame(Sample=Sample,Group)

        dds <- DESeqDataSetFromMatrix(countData = round(dat.log2), colData = sampleData, design = ~ Group)
        dds <- DESeq(dds)
        dds_results <- results(dds, contrast = c("Group","cancer", "control")) #Group cancer vs control
        Result2 <- cbind(Data[,2:4],dds_results, dat.log2, CountData2)
        write.csv(Result2, file="DESeq2_black_white.csv", quote=F)
           
       ### DEG Analysis - 2. limma ###
        library(limma)
        design <- model.matrix(~0 +Group)
        fit <- lmFit(dat.log2, design)
        contr <- makeContrasts(Groupcancer - Groupcontrol, levels = colnames(coef(fit)))
        tmp <- contrasts.fit(fit, contr)
        tmp <- eBayes(tmp)
        result <- topTable(tmp, sort.by = "P", n = Inf)
        write.table(result,file="TMMnorm_count.log2.lm.txt",sep="\t", col.names=NA,quote=F)
}

TMM_Wilcox <- function(){
        library(dplyr)
        dat.log3 <- as.data.frame(dat.log2)
        dat.log3$Avg.control <- rowMeans(dat.log3[,1:20])
        dat.log3$Avg.cancer <- rowMeans(dat.log3[,21:56])
        dat.log3 <- mutate(dat.log3, log2FC = (Avg.cancer - Avg.control))

        dat.log3 <- mutate(dat.log3, wilcox.test = apply(dat.log3[, 1:56], 1, function(x) wilcox.test(x[1:20], x[21:56])$p.value))      #wilcox.test(unlist(dat.log3[1, 1:20]), unlist(dat.log3[1, 21:56]))$p.value
        dat.log3 <- mutate(dat.log3, bonferroni = p.adjust(dat.log3$wilcox.test, method="bonferroni"))
        dat.log3 <- mutate(dat.log3, fdr = p.adjust(dat.log3$wilcox.test, method="fdr"))
     
        #Clipper
        library(Clipper)
        control <- dat.log2[,1:20]
        cancer <- dat.log2[,21:56]
        Clipper <- Clipper(score.exp = cancer, score.back = control, analysis = "differential",FDR=0.2)
        Clipper$discoveries
        Result <- cbind(Data[,2:4], dat.log3[,c(57:62)], Clipper$q,Clipper_p$q, dat.log3[,c(1:56)], CountData2)
        write.csv(Result, file="TMMlog2.Wilcox.csv", quote=F)
}


