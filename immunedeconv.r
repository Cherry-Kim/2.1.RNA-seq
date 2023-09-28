        #install.packages("remotes")
        #remotes::install_github("grst/immunedeconv")
        library(devtools)
        library(readr)
        library(dplyr)
        library(ggplot2)
        library(tidyr)
        library(immunedeconv)
        library(tibble)

        data <- read.table("TMMnorm.t.txt",header=T, stringsAsFactors=F)
        #Symbol sample1 sample2 sample3
        #A1BG  3.971484  3.560819  3.8813731
        
        counts <- data[,-1]
        rownames(counts) <- as.character(data[,1])
        res_mcp_counter = deconvolute(counts, "mcp_counter")
        res_mcp_counter = deconvolute(counts, "mcp_counter")
        #png(OUTPUT)
        res_mcp_counter %>%
                gather(sample, score, -cell_type) %>%
                ggplot(aes(x=sample, y=score, color=cell_type)) +
                geom_point(size=0.5) +
                facet_wrap(~cell_type, scales="free_x", ncol=3) +
                scale_color_brewer(palette="Paired", guide=FALSE) +
                coord_flip() +
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        #dev.off()
        write.table(res_mcp_counter, "res_mcp_counter.out", col.names=NA, row.names=T,quote=F, sep='\t')
        write.table(t(res_mcp_counter), "res_mcp_counter.t.out", col.names=NA, row.names=T,quote=F, sep='\t')
