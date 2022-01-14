import string,os,glob,sys

co=0
fp=glob.glon("*_1.fq.gz")
for fname in fp:
	co+=1
	sample=fname.split('_1.fq.gz')[0]
	
	print co, sample, "STEP1. Mapping using HISAT2"
	os.system('hisat2 -p 20 --dta --rna-strandness RF --novel-splicesite-outfile '+sample+'.novel_splicesites.txt -x hisat2-2.2.1/example/index/22_20-21M_snp  -1 '+sample+'_1.fq.gz -2 '+sample+'_2.fq.gz -S '+sample+'.sam')
	# Hisat index : https://daehwankimlab.github.io/hisat2/download/

	os.system('samtools view -@ 48 -bhS '+sample+'.sam > '+sample+'.bam')
	os.system('samtools sort -@ 48 -o '+sample+'.sroted.bam '+sample+'.bam')

	print co, sample, "STEP2. Expression profiling using STRINGTIE"
	os.system('stringtie '+sample+'.sroted.bam --rf -p 20 -l '+sample+' -o '+sample+'.gtf -G hg19.refGene.gtf -A '+sample+'.gene_abund.anno.tab -C '+sample+'.cov_refs.anno.gtf -e')
	#-A : gene abundance estimation output file 
	#-C : output a file with reference transcripts that are covered by reads


###### STEP3. DEG in R ########
library(tximport)
#devtools::install_github("tidyverse/readr")
library(tidyverse)
library(readr)
library(DESeq2)

dir="/home/Sample_hy/"
file_list <- dir(pattern="sample*")
files <- c()
for (i in 1:length(file_list)) {
+ temp <- paste(dir,file_list[i],"/t_data.ctab",sep="")
+ files <- append(files, temp)
+ }

tmp <- read_tsv(files[1])
tx2gene <- tmp[, c("t_name", "gene_name")]

txi <- tximport(files, type = "stringtie", tx2gene = tx2gene)
write.table(txi, "test.txt")

sampleNames <- c("sample1","sample2","sample3","sample4")
sampleGroup <- c("T","T","N","N")
sampleTable <- data.frame(sampleName=sampleNames, type=sampleGroup)
rownames(sampleTable) <- colnames(txi$counts)

dds  <- DESeqDataSetFromTximport(txi, colData = sampleTable, design = ~ type)
dds <- DESeq(dds)

