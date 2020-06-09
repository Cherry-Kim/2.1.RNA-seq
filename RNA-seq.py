#!/usr/bin/python

import os,string,glob

sample1='Homo-M2_L'
sample2='Homo-M2_R'
sample3='HT-M3_L'
sample4='HT-M3_R'

REF='/BiO/REF/mm10/mm10'
gtf='/BiO/REF/mm10/mm.GRCm38.mod.gtf'
gtf_ind='/BiO/REF/mm10/mm10/mm.GRCm38.mod'

co=0
fp=glob.glob('*_1.fastq.gz')
for fname in fp:
    co+=1
    a=string.split(fname,'_')
    Sample=a[0]+'_'+a[1]
    print Sample

    print '### STEP0. Trimming with Trimmomatic  ###'
    os.system('java -jar /BiO/BioTools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 24 -phred33 '+Sample+'_1.fastq.gz '+Sample+'_2.fastq.gz '+Sample+'.r1.trim.fq '+Sample+'.r1.unpair.fq '+Sample+'.r2.trim.fq '+Sample+'.r2.unpair.fq ILLUMINACLIP:/BiO/BioTools/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:0:151:0 LEADING:0 TRAILING:0 SLIDINGWINDOW:4:0 MINLEN:0')

    print '### STEP1. Mapping with TopHat ###'
    os.system('mkdir -p results_tophat/'+Sample)
    os.system('tophat2 -p 24 -o results_tophat/'+Sample+' -G '+gtf+' --transcriptome-index '+gtf_ind+' '+REF+' '+Sample+'.r1.trim.fq '+Sample+'.r2.trim.fq')

    print '### STEP2. Assembly with Cufflinks ###'
    os.system('mkdir -p results_cufflinks/'+Sample)
    os.system('cufflinks -p 24 -o results_cufflinks/'+Sample+' results_tophat/'+Sample+'/accepted_hits.bam')

'''
print '### STEP3. Merge with Cuffmerge ###'
if not os.path.isdir('results_cuffmerge'):
    os.mkdir('results_cuffmerge')
os.system('echo "results_cufflinks/'+sample1+'/transcripts.gtf" > assemblies.txt')
os.system('echo "results_cufflinks/'+sample2+'/transcripts.gtf" >> assemblies.txt')
os.system('echo "results_cufflinks/'+sample3+'/transcripts.gtf" >> assemblies.txt')
os.system('echo "results_cufflinks/'+sample4+'/transcripts.gtf" >> assemblies.txt')

os.system('cuffmerge -g '+gtf+' -s '+REF+'.fa -o results_cuffmerge/ -p 24 assemblies.txt')

##os.system('cuffquant -p 24 -o results_cufflinks/'+Sample+' results_cufflinks/merged.gtf results_tophat/'+Sample+'/accepted_hits.bam')
##os.system('cuffnorm /BiO/REF/mm10/mm.GRCm38.mod.gtf -p 24 -o results_cuffnorm/ results_tophat/'+sample1+'/accepted_hits.bam results_tophat/'+sample2+'/accepted_hits.bam')

print '### STEP4. Cuffdiff ###'
os.system('cuffdiff -b '+REF+'.fa -L '+sample1+','+sample2+' -u results_cuffmerge/merged.gtg -o results_cuffdiff/ -p 24 results_tophat/'+sample1+'/accepted_hits.bam results_tophat/'+sample2+'/accepted_hits.bam')
os.system('awk \'{if ($14 == "yes") print}\' results_cuffdiff/gene_exp.diff > results_cuffdiff/gene_exp_filt2.diff')
##os.system('cat results_cuffdiff/gene_exp.diff | awk \'($14 == "yes")\' > results_cuffdiff/gene_exp_filt.diff')
'''
