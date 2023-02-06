#! /usr/bin/env Rscript
#####################################################
# File Name: STAR-RSEM-Immunedecov.py
#####################################################
import string,sys,glob,os
from rpy2 import robjects as ro
r = ro.r

#-----------------------------------------------------------------------------------------------
wget http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.gtf.gz
mv Homo_sapiens.GRCh38.101.gtf GRCh38.gtf
mv Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38.fa
#-----------------------------------------------------------------------------------------------
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz
#wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz
#gunzip *.gz
#STAR --runMode genomeGenerate --runThreadN 48 --genomeDir . --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v34.primary_assembly.annotation.gtf  --sjdbOverhang 100
#rsem-prepare-reference --gtf gencode.v34.primary_assembly.annotation.gtf --star --star-path /home/STAR/ -p 48 GRCh38.primary_assembly.genome.fa GRCh38
#-----------------------------------------------------------------------------------------------

# apt install rna-star
# Generating genome indexes
STAR --runMode genomeGenerate  --runThreadN 48 --genomeDir .  --genomeFastaFiles GRCh38.fa --outFileNamePrefix GRCh38 --sjdbGTFfile GRCh38.gtf --sjdbOverhang 100
rsem-prepare-reference --gtf GRCh38.gtf --star --star-path /STAR-2.7.2d/bin/Linux_x86_64 -p 48 /REF/STAR/GRCh38.fa GRCh38

STAR --genomeDir /REF/STAR/ --readFilesIn sample1_1.fastq.gz sample1_2.fastq.gz --readFilesCommand zcat --outTmpDir /home/temp_rsem --outFileNamePrefix test --quantMode TranscriptomeSAM
rsem-calculate-expression --bam  --paired-end  testAligned.toTranscriptome.out.bam  /REF/STAR/GRCh38 sample

### STEP3. PCA ###
from rpy2 import robjects as ro
r = ro.r
r.source("edgeR.R")

import string,sys,glob,os
  
def STEP1_Preprocessing(sample):
    #os.system('fastqc -t 48 --nogroup '+sample+'.R1.fastq.gz')
    #os.system('fastqc -t 48 --nogroup '+sample+'.R2.fastq.gz')

    os.system('trim_galore --paired  --gzip --fastqc  -o trim_galore/ '+sample+'.R1.fastq.gz '+sample+'.R2.fastq.gz  -q 20 --length 20')
    ##os.system('multiqc ./')

def STEP2_STAR(sample,REF):
    os.system('STAR --runThreadN 48  --genomeDir '+REF+' --readFilesIn trim_galore/'+sample+'.R1_val_1.fq.gz trim_galore/'+sample+'.R2_val_2.fq.gz --readFilesCommand zcat --outTmpDir TEMP_'+sample+'  --outFileNamePrefix '+sample+' --quantMode TranscriptomeSAM')

def STEP3_RSEM(sample,REF,SP):
    os.system('/BIO1/RSEM-1.3.3/rsem-calculate-expression -p 48  --bam --paired-end '+sample+'.Aligned.toTranscriptome.out.bam '+REF+SP+' '+sample)

def Main():
    REF,SP = 'REF/mm10/','mm10'   #REF,SP = 'REF/GRCh38/','GRCh38' 
    co = 0
    file_list=os.listdir('./')
    f_list=[file for file in file_list if file.endswith(".R1.fastq.gz")]
    for fname in f_list:
        co += 1
        sample = fname.split('.R1.fastq.gz')[0]
        print ("###",co,sample)
        #STEP1_Preprocessing(sample)
        STEP2_STAR(sample,REF)
        #STEP3_RSEM(sample,REF,SP)
Main()

