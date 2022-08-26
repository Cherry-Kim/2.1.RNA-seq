import string,sys,os,glob
def STEP0_GRCh38(STAR_HUMAN_INDEX,REFERENCE_GENOME,GTF,remapNCBI):
    os.system('wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz')
    os.system('wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz')
    os.system('gunzip GRCh38.primary_assembly.genome.fa.gz gencode.v34.primary_assembly.annotation.gtf.gz')
    os.system('wget https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz')
    os.system('wget https://ftp.ncbi.nih.gov/snp/archive/b155/VCF/GCF_000001405.39.gz.tbi')
   
    #2. Reference data setup 
    os.system('samtools faidx '+ REFERENCE_GENOME)
    os.system('java -jar picard.jar CreateSequenceDictionary R='+REFERENCE_GENOME)
    os.system('STAR-2.5.4b/bin/Linux_x86_64/STAR --runMode genomeGenerate  --runThreadN 48 --genomeDir '+STAR_HUMAN_INDEX+' --genomeFastaFiles '+ REFERENCE_GENOME +'  --sjdbGTFfile '+GTF+'  --sjdbOverhang 100')

    os.system("""awk '{if($3=="exon") {print $1"\t"$4-100"\t"$5+100"\t"substr($16,2,length($16)-3)}}' """+GTF+""" | sort -k 1,1 -k2,2n | bgzip > GRCh38_exome.bed.gz""")
    os.system('tabix GRCh38_exome.bed.gz')
 
    os.system('bcftools annotate --threads 48 --output-type z --rename-chrs '+remapNCBI+' --output dbSNPbuild154Renamed.vcf.gz GCF_000001405.39.gz')
    os.system('tabix dbSNPbuild154Renamed.vcf.gz')

def STEP1_STAR(STAR_HUMAN_INDEX,SAMPLE):
    os.system('trim_galore --paired -j 48 --gzip  -o trim_galore/ --fastqc  '+ SAMPLE +'_R1.fastq.gz '+ SAMPLE +'_R2.fastq.gz  -q 20 --length 20')
    os.system('STAR --runThreadN 48 --genomeDir '+ STAR_HUMAN_INDEX+' --readFilesIn trim_galore/'+ SAMPLE +'_R1_val_1.fq.gz trim_galore/'+ SAMPLE +'_R2_val_2.fq.gz  --outFileNamePrefix '+SAMPLE +' --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outSAMattributes NM --twopassMode Basic --outFilterMultimapNmax 1 --outFilterMismatchNoverLmax 0.1')

def STEP2_DUPLICATE(SAMPLE):
    os.system('gatk MarkDuplicates --CREATE_INDEX true -I '+SAMPLE+'Aligned.sortedByCoord.out.bam  -O '+SAMPLE+'marked_duplicates.bam --VALIDATION_STRINGENCY SILENT  -M '+SAMPLE+'marked_dup_metrics.txt')

def STEP3_SplitNCigarReads(SAMPLE,REFERENCE_GENOME):
    os.system('gatk SplitNCigarReads  -R '+ REFERENCE_GENOME +' -I '+SAMPLE+'marked_duplicates.bam -O '+SAMPLE+'splitN.bam')

def STEP4_BQSR(SAMPLE,REFERENCE_GENOME,DBSNP):
    os.system('gatk AddOrReplaceReadGroups --CREATE_INDEX true -I '+SAMPLE+'splitN.bam -O '+SAMPLE+'.grouped.bam --RGID rnasq --RGLB lb --RGPL illumina --RGPU pu --RGSM '+SAMPLE)
    os.system('gatk BaseRecalibrator -I '+SAMPLE+'.grouped.bam  -R '+REFERENCE_GENOME+' --known-sites '+DBSNP+' -O '+SAMPLE+'.recal_data.table')
    os.system('gatk ApplyBQSR -R '+REFERENCE_GENOME+' -I '+SAMPLE+'.grouped.bam --use-original-qualities --add-output-sam-program-record --bqsr-recal-file '+SAMPLE+'.recal_data.table -O '+SAMPLE+'.recal_output.bam')

