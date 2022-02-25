TMM_norm <- function(){
        library(edgeR)
        data <- read.table("rsem.txt", header=T)
        head(data[1:3,1:3])
        ##----------------------------------------------------
                #      sample1_1 sample1_2 sample1_3 sample2_1 sample2_2 sample2_3
        #3715        0       1       0       0       0       0
        #6845        0       0       0       0       0      13
        #----------------------------------------------------

        group <- colnames(data[3:ncol(data)])
        d <- DGEList(counts = data[,3:ncol(data)], group=group)
        d <- calcNormFactors(d, method="TMM")   #Normalize your raw data using the TMM method
        cps <- cpm(d, normalized.lib.sizes=TRUE)

        re.dat <- cps
        dat.log2 <- re.dat
        for (i in 1:ncol(re.dat)){
                dat.log2[,i] <- log(re.dat[,i]+1,2)
        }

        rownames(dat.log2) <- data[,1]
        write.table(dat.log2, file="output_TMMnorm.txt",sep="\t", col.names=NA, quote=FALSE)

#       rownames(cps) <- as.character(data[,1])

#       itr <- c()
#       for(i in 1:nrow(dat.log2)){
#               if(sd(dat.log2[i,])==0)
#                       {itr <- c(itr,i)}
#       }
#       length(itr)
#       in.dat <- dat.log2[-itr,]

        return(dat.log2)
}

TMM_norm()
