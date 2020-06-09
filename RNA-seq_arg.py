#!/usr/bin/python
#REF = /BiO/Reference/galGal5/galGal5_mod2.fa
#REF2 = /BiO/Reference/galGal5/Gallus_gallus_mod2.gtf

import os,sys

def help():
	print ''
	print '[USAGE] *.py [sample]'
	print ''
#python 1.RNA-seq_CDL_CLL.py COT2-1D COT2-2D COT2-3D COT1-1L COT1-2L COT1-3L

def main(s1,s2,s3,s4,s5,s6):
#	print '### STEP0. Trimming with Trimmomatic S1 ###'
#	os.system('java -jar /BiO/hykim/NGS/1.RNA-seq/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 24 -phred33 /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.fastq.gz /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.fastq.gz /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.unpair.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.unpair.fastq ILLUMINACLIP:/BiO/hykim/NGS/1.RNA-seq/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'%(s1,s1,s1,s1,s1,s1))

#	print '### STEP0. Trimming with Trimmomatic S2 ###'
#	os.system('java -jar /BiO/hykim/NGS/1.RNA-seq/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 24 -phred33 /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.fastq.gz /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.fastq.gz /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.unpair.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.unpair.fastq ILLUMINACLIP:/BiO/hykim/NGS/1.RNA-seq/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'%(s2,s2,s2,s2,s2,s2))

#	print '### STEP0. Trimming with Trimmomatic S3 ###'
#	os.system('java -jar /BiO/hykim/NGS/1.RNA-seq/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 24 -phred33 /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.fastq.gz /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.fastq.gz /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.unpair.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.unpair.fastq ILLUMINACLIP:/BiO/hykim/NGS/1.RNA-seq/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'%(s3,s3,s3,s3,s3,s3))

#	print '### STEP0. Trimming with Trimmomatic S4 ###'
#	os.system('java -jar /BiO/hykim/NGS/1.RNA-seq/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 24 -phred33 /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.fastq.gz /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.fastq.gz /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.unpair.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.unpair.fastq ILLUMINACLIP:/BiO/hykim/NGS/1.RNA-seq/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'%(s4,s4,s4,s4,s4,s4))

#	print '### STEP0. Trimming with Trimmomatic S5 ###'
#	os.system('java -jar /BiO/hykim/NGS/1.RNA-seq/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 24 -phred33 /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.fastq.gz /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.fastq.gz /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.unpair.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.unpair.fastq ILLUMINACLIP:/BiO/hykim/NGS/1.RNA-seq/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'%(s5,s5,s5,s5,s5,s5))

#	print '### STEP0. Trimming with Trimmomatic S6 ###'
#	os.system('java -jar /BiO/hykim/NGS/1.RNA-seq/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 24 -phred33 /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.fastq.gz /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.fastq.gz /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.unpair.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.unpair.fastq ILLUMINACLIP:/BiO/hykim/NGS/1.RNA-seq/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36'%(s6,s6,s6,s6,s6,s6))
########################################################################

#	print '### STEP1. Mapping with TopHat Sample1 ###'
#	os.system('mkdir -p results_tophat/%s'%(s1))
#	os.system('tophat -p 24 -o results_tophat/%s -G /BiO/Reference/galGal5/Gallus_gallus_mod2.gtf /BiO/Reference/galGal5/galGal /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.trim.fastq'%(s1,s1,s1))

#	print '### STEP1. Mapping with TopHat Sample2 ###'
#	os.system('mkdir -p results_tophat/%s'%(s2))
#	os.system('tophat -p 24 -o results_tophat/%s -G /BiO/Reference/galGal5/Gallus_gallus_mod2.gtf /BiO/Reference/galGal5/galGal /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.trim.fastq'%(s2,s2,s2))

#	print '### STEP1. Mapping with TopHat Sample3 ###'
#	os.system('mkdir -p results_tophat/%s'%(s3))
#	os.system('tophat -p 24 -o results_tophat/%s -G /BiO/Reference/galGal5/Gallus_gallus_mod2.gtf /BiO/Reference/galGal5/galGal /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.trim.fastq'%(s3,s3,s3))

#	print '### STEP1. Mapping with TopHat Sample4 ###'
#	os.system('mkdir -p results_tophat/%s'%(s4))
#	os.system('tophat -p 24 -o results_tophat/%s -G /BiO/Reference/galGal5/Gallus_gallus_mod2.gtf /BiO/Reference/galGal5/galGal /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.trim.fastq'%(s4,s4,s4))

#	print '### STEP1. Mapping with TopHat Sample5 ###'
#	os.system('mkdir -p results_tophat/%s'%(s5))
#	os.system('tophat -p 24 -o results_tophat/%s -G /BiO/Reference/galGal5/Gallus_gallus_mod2.gtf /BiO/Reference/galGal5/galGal /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.trim.fastq'%(s5,s5,s5))
	
#	print '### STEP1. Mapping with TopHat Sample6 ###'
#	os.system('mkdir -p results_tophat/%s'%(s6))
#	os.system('tophat -p 24 -o results_tophat/%s -G /BiO/Reference/galGal5/Gallus_gallus_mod2.gtf /BiO/Reference/galGal5/galGal /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_1.trim.fastq /home/hykim/1.RNA-seq/1.L-Met_broiler_Heat/Liver/0.DATA/%s_2.trim.fastq'%(s6,s6,s6))
###################################################################################
	
#	print '### STEP2. Assembly with Cufflinks S1 ###'
#	os.system('mkdir -p results_cufflinks/%s'%(s1))
#	os.system('cufflinks -p 24 -o results_cufflinks/%s results_tophat/%s/accepted_hits.bam'%(s1,s1))
	
#	print '### STEP2. Assembly with Cufflinks S2 ###'
#	os.system('mkdir -p results_cufflinks/%s'%(s2))
#	os.system('cufflinks -p 24 -o results_cufflinks/%s results_tophat/%s/accepted_hits.bam'%(s2,s2))

#	print '### STEP2. Assembly with Cufflinks S3 ###'
#	os.system('mkdir -p results_cufflinks/%s'%(s3))
#	os.system('cufflinks -p 24 -o results_cufflinks/%s results_tophat/%s/accepted_hits.bam'%(s3,s3))

#	print '### STEP2. Assembly with Cufflinks S4 ###'
#	os.system('mkdir -p results_cufflinks/%s'%(s4))
#	os.system('cufflinks -p 24 -o results_cufflinks/%s results_tophat/%s/accepted_hits.bam'%(s4,s4))

#	print '### STEP2. Assembly with Cufflinks S5 ###'
#	os.system('mkdir -p results_cufflinks/%s'%(s5))
#	os.system('cufflinks -p 24 -o results_cufflinks/%s results_tophat/%s/accepted_hits.bam'%(s5,s5))

#	print '### STEP2. Assembly with Cufflinks S6 ###'
#	os.system('mkdir -p results_cufflinks/%s'%(s6))
#	os.system('cufflinks -p 24 -o results_cufflinks/%s results_tophat/%s/accepted_hits.bam'%(s6,s6))
##################################################################################

#	print '### STEP3. Merge with Cuffmerge ###'
#	os.system('echo "results_cufflinks/%s/transcripts.gtf" > assemblies.txt'%(s1))
#	os.system('echo "results_cufflinks/%s/transcripts.gtf" >> assemblies.txt'%(s2))
#	os.system('echo "results_cufflinks/%s/transcripts.gtf" >> assemblies.txt'%(s3))
#	os.system('echo "results_cufflinks/%s/transcripts.gtf" >> assemblies.txt'%(s4))
#	os.system('echo "results_cufflinks/%s/transcripts.gtf" >> assemblies.txt'%(s5))
#	os.system('echo "results_cufflinks/%s/transcripts.gtf" >> assemblies.txt'%(s6))

#	os.system('cuffmerge -g /BiO/Reference/galGal5/Gallus_gallus_mod2.gtf -s /BiO/Reference/galGal5/galGal5_mod2.fa -o results_cufflinks/ -p 24 assemblies.txt')
##################################################################################

#	print '### STEP4. Cuffdiff ###'
#	os.system('cuffdiff -b /BiO/Reference/galGal5/galGal5_mod2.fa -L C_D-Met,C_L-Met -u results_cufflinks/merged.gtf -o results_cuffdiff/ -p 24 results_tophat/%s/accepted_hits.bam,results_tophat/%s/accepted_hits.bam,results_tophat/%s/accepted_hits.bam results_tophat/%s/accepted_hits.bam,results_tophat/%s/accepted_hits.bam,results_tophat/%s/accepted_hits.bam'%(s1,s2,s3,s4,s5,s6))


if __name__ == '__main__':

	if len(sys.argv) < 7:
		sys.exit(help())

	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6])

