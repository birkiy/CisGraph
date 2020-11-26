

fastqc



STAR --runThreadN 6 \
--runMode genomeGenerate \
--genomeDir igenome \
--genomeFastaFiles ~/genomeAnnotations/hg19.fa










#!/bin/bash
#SBATCH --job-name=star
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=8
#SBATCH --export=all
#SBATCH -p long

home=/home/ualtintas/STARexample

STAR --genomeDir $home/igenome \
--runThreadN 8 \
--readFilesCommand zcat \
--readFilesIn $home/fastq/SRR5796680_1.fastq.gz fastq/SRR5796680_2.fastq.gz \
--outFileNamePrefix $home/results/SRR5796680  \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--runMode alignReads

STAR --genomeDir $home/igenome \
--runThreadN 8 \
--readFilesCommand zcat \
--readFilesIn $home/fastq/SRR5796681_1.fastq.gz fastq/SRR5796681_2.fastq.gz \
--outFileNamePrefix $home/results/SRR5796681  \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--runMode alignReads

STAR --genomeDir $home/igenome \
--runThreadN 8 \
--readFilesCommand zcat \
--readFilesIn $home/fastq/SRR5796682_1.fastq.gz fastq/SRR5796682_2.fastq.gz \
--outFileNamePrefix $home/results/SRR5796682  \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--runMode alignReads

STAR --genomeDir $home/igenome \
--runThreadN 8 \
--readFilesCommand zcat \
--readFilesIn $home/fastq/SRR5796683_1.fastq.gz fastq/SRR5796683_2.fastq.gz \
--outFileNamePrefix $home/results/SRR5796683  \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--runMode alignReads




#!/bin/bash
#SBATCH --job-name=cufflinks
#SBATCH --mem-per-cpu=30G
#SBATCH --cpus-per-task=2
#SBATCH --export=all
#SBATCH -p long

home=/home/ualtintas/STARexample


cd $home/results/SRR5796680
cufflinks $home/results/SRR5796680/SRR5796680Aligned.sortedByCoord.out.bam --library-type fr-firststrand

cd $home/results
mkdir SRR5796681
mv $home/results/SRR5796681Aligned.sortedByCoord.out.bam $home/results/SRR5796681/.
cd $home/results/SRR5796681/
cufflinks $home/results/SRR5796681/SRR5796681Aligned.sortedByCoord.out.bam --library-type fr-firststrand

cd $home/results
mkdir SRR5796682
mv $home/results/SRR5796682Aligned.sortedByCoord.out.bam $home/results/SRR5796682/.
cd $home/results/SRR5796682/
cufflinks $home/results/SRR5796682/SRR5796682Aligned.sortedByCoord.out.bam --library-type fr-firststrand

cd $home/results
mkdir SRR5796683
mv $home/results/SRR5796683Aligned.sortedByCoord.out.bam $home/results/SRR5796683/.
cd $home/results/SRR5796683/
cufflinks $home/results/SRR5796683/SRR5796683Aligned.sortedByCoord.out.bam --library-type fr-firststrand


assemblies.txt


/home/ualtintas/STARexample/results/SRR5796680/transcripts.gtf
/home/ualtintas/STARexample/results/SRR5796681/transcripts.gtf
/home/ualtintas/STARexample/results/SRR5796682/transcripts.gtf
/home/ualtintas/STARexample/results/SRR5796683/transcripts.gtf

cuffmerge -s /home/ualtintas/genomeAnnotations/hg19.fa -p 8 assemblies.txt -o /home/ualtintas/STARexample/results/out.stdout


home=/home/ualtintas/STARexample

cuffquant tmp_meta_asm.combined.gtf \
  $home/results/SRR5796680/SRR5796680Aligned.sortedByCoord.out.bam \
  $home/results/SRR5796681/SRR5796681Aligned.sortedByCoord.out.bam \
  $home/results/SRR5796682/SRR5796682Aligned.sortedByCoord.out.bam \
  $home/results/SRR5796683/SRR5796683Aligned.sortedByCoord.out.bam \
  -p 8 --library-type fr-firststrand


cuffdiff -o diff.out --labels AR,AR,IgG,IgG tmp_meta_asm.combined.gtf \
  $home/results/SRR5796680/SRR5796680Aligned.sortedByCoord.out.bam \
  $home/results/SRR5796681/SRR5796681Aligned.sortedByCoord.out.bam \
  $home/results/SRR5796682/SRR5796682Aligned.sortedByCoord.out.bam \
  $home/results/SRR5796683/SRR5796683Aligned.sortedByCoord.out.bam \
  -p 8 --library-type fr-firststrand
