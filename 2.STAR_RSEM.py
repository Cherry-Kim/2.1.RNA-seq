#! /usr/bin/env Rscript

import string,sys,glob,os

STAR --runMode genomeGenerate --runThreadN 24 --genomeDir /home/hykim/REF/Human/hg38/STAR --genomeFastaFiles GRCh38.fa

STAR --genomeDir /home/hykim/REF/Human/hg38/STAR/ --readFilesIn PM-AA-0052-N-A1_1.fastq.gz PM-AA-0052-N-A1_2.fastq.gz --readFilesCommand zcat --outTmpDir /home/hykim/Project/Cancer/Colon/RNA-seq/RAW_fastq_files/temp_rsem7 --outFileNamePrefix test2 --quantMode TranscriptomeSAM

rsem-prepare-reference --gtf /home/hykim/REF/Human/hg38/STAR/GRCh38.gtf --star --star-path /home/program/STAR-2.7.2d/bin/Linux_x86_64 -p 8 /home/hykim/REF/Human/hg38/STAR/GRCh38.fa GRCh38

rsem-calculate-expression --bam  --paired-end  testAligned.toTranscriptome.out.bam  /home/hykim/REF/Human/hg38/STAR/GRCh38 b

### STEP3. PCA ###
from rpy2 import robjects as ro
r = ro.r
r.source("edgeR.R")
