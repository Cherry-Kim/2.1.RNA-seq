#! /usr/bin/env Rscript

import string,sys,glob,os

STAR --runMode genomeGenerate --runThreadN 24 --genomeDir /REF/STAR --genomeFastaFiles GRCh38.fa

STAR --genomeDir /REF/STAR/ --readFilesIn sample1_1.fastq.gz sample1_2.fastq.gz --readFilesCommand zcat --outTmpDir /home/temp_rsem --outFileNamePrefix test --quantMode TranscriptomeSAM

rsem-prepare-reference --gtf /REF/STAR/GRCh38.gtf --star --star-path /STAR-2.7.2d/bin/Linux_x86_64 -p 8 /REF/STAR/GRCh38.fa GRCh38

rsem-calculate-expression --bam  --paired-end  testAligned.toTranscriptome.out.bam  /REF/STAR/GRCh38 b

### STEP3. PCA ###
from rpy2 import robjects as ro
r = ro.r
r.source("edgeR.R")
