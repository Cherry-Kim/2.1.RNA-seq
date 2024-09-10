volcano <- function(){
        data <- read.table('rsem.counts.DESeq2.NA.txt',header=T,sep='\t')
        topT <- data[,c(2,61:65)]
        # hgnc_symbol log2FoldChange     lfcSE         stat      pvalue       padj
#1       NLRP2    0.006522454 0.6942832  0.009394515 0.992504372 0.99957973
        colnames(topT) <- c("symbol","log2FC","lfcSE","stat","pvalue","padj")
      
        cut_lfc <- 0.58
        cut_fdr <- 0.05

        par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
        with(topT, plot(log2FC, -log10(pvalue), pch=20, main="PR vs. NR (ref=NR)", col='grey', cex=1.0, xlab=bquote(~log2(Fold~Change)), ylab=bquote(~-log10(pvalue))) )
        with(subset(topT, pvalue<0.05 & log2FC>cut_lfc), points(log2FC, -log10(pvalue), pch=20, col='red', cex=1.5))
        with(subset(topT, pvalue<0.05 & log2FC <(-cut_lfc)), points(log2FC, -log10(pvalue), pch=20, col='blue', cex=1.5))
        with(subset(topT, topT$Sig == TRUE), points(log2FC, -log10(pvalue), pch=20, col='green', cex=1.5))


        abline(v=0, col='black', lty=3, lwd=1.0)
        abline(v=-cut_lfc, col='black', lty=4, lwd=2.0)
        abline(v=cut_lfc, col='black', lty=4, lwd=2.0)
        abline(h=-log10(max(topT$pvalue[topT$pvalue<cut_pvalue], na.rm=TRUE)), col='black', lty=5, lwd=2.0)
        abline(h=-log10(max(topT$pvalue[topT$pvalue<0.05], na.rm=TRUE)), col='black', lty=2, lwd=2.0)

        ## label        
        #topT$gene <- rownames(topT)
        topT$gene <- topT$symbol
        Sig <- topT$Sig
        text(topT$log2FC[Sig], -log10(topT$pvalue)[Sig], vfont = NULL, lab=(topT$gene)[Sig], cex=0.6, col='black', pos=4)
        topT$gene[Sig]

  ## label2
  topT$Sig2 <- ifelse( (topT$pvalue <0.05 & topT$log2FC<(-cut_lfc)) | (topT$pvalue <0.05 & topT$log2FC>cut_lfc), TRUE, FALSE)
  table(topT$Sig2)
  topT$gene <- rownames(topT)
  Sig2 <- topT$Sig2
  text(topT$log2FC[Sig2], -log10(topT$pvalue)[Sig2], vfont = NULL, lab=(topT$gene)[Sig2], cex=0.6, col='black', pos=4)
  topT$gene[Sig2]
}

