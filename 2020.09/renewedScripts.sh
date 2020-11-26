###########################################
#                                         #
#          Umut Berkay Altıntaş           #
#                                         #
###########################################




#### DOWNLOADING #### => An optimized download script for many samples.

* Be aware that script uses a sampleList file that contains, at 1st col SRR IDs, at 2nd col SRX IDs with tab delimited.

* Notice that this script can download both paired and single end libraries. However, when single, will give error that could not find "_2.fastq" file. But still should work.

* Note that to use of high cpus as some files are hard to download, because of time-out error. You can try higher but server may halt, so do not push too much. Use wget FTP.

sbatch ~/scripts/cpDownload.sh CpDownList.txt



#!/bin/bash
#SBATCH --job-name=cpDownload
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=15
#SBATCH -p cosmos
#SBATCH --time=24:00:00

sampleList=$1

inpPath=/kuacc/users/ualtintas20/0-Downloading
outPath=/kuacc/users/ualtintas20/0-Downloading

while IFS="$(printf '\t')" read -r SRR SRX ; do

	parallel-fastq-dump --sra-id $SRR \
	--outdir $outPath/Chip-Seq/DENEME \
	--threads 15 --split-files ;

	cat $outPath/Chip-Seq/DENEME/$SRR"_1.fastq" >> $outPath/Chip-Seq/DENEME/$SRX"_1.fastq";
	cat $outPath/Chip-Seq/DENEME/$SRR"_2.fastq" >> $outPath/Chip-Seq/DENEME/$SRX"_2.fastq";

	rm $outPath/Chip-Seq/DENEME/$SRR"_1.fastq" ;
	rm $outPath/Chip-Seq/DENEME/$SRR"_2.fastq" ;

done < <(cat $inpPath/$sampleList)




#!/bin/bash
#SBATCH --job-name=cpDownload
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=15
#SBATCH -p cosmos
#SBATCH --time=24:00:00

sampleList=$1

inpPath=/kuacc/users/ualtintas20/0-Downloading
outPath=/kuacc/users/ualtintas20/0-Downloading

while IFS="$(printf '\t')" read -r SRX SRR LibraryType ExperimentType ; do

	if [ $LibraryType == "SINGLE" ] ; then


		parallel-fastq-dump --split-files --threads 15 --sra-id $SRR --outdir $inpPath/Downloads/$ExperimentType/$LibraryType;

		cat $inpPath/Downloads/$ExperimentType/$LibraryType/$SRR"_1.fastq" >> $outPath/Downloads/$ExperimentType/$LibraryType/$SRX"_1.fastq" ;

		rm $inpPath/Downloads/$ExperimentType/$LibraryType/$SRR"_1.fastq";


	elif [ $LibraryType == "PAIRED" ] ; then


		parallel-fastq-dump --split-files --threads 15 --sra-id $SRR --outdir $inpPath/Downloads/$ExperimentType/$LibraryType;

		cat $inpPath/Downloads/$ExperimentType/$LibraryType/$SRR"_1.fastq" >> $outPath/Downloads/$ExperimentType/$LibraryType/$SRX"_1.fastq";
		cat $inpPath/Downloads/$ExperimentType/$LibraryType/$SRR"_2.fastq" >> $outPath/Downloads/$ExperimentType/$LibraryType/$SRX"_2.fastq";

		rm $inpPath/Downloads/$ExperimentType/$LibraryType/$SRR"_1.fastq";
		rm $inpPath/Downloads/$ExperimentType/$LibraryType/$SRR"_2.fastq";


	fi


done < <(cat $sampleList)













#### MAPPING ####

* Can be used for both single and paired libraries, require LibraryType indication in MapList file.

* Be aware that this script uses MapList file which contains at 1st col SRX IDs, at 2nd col LibraryType (Paired, Single).

* Awarage server usage. Tried with 12-13 samples per run.

sbatch ~/scripts/Mapping.sh

#!/bin/bash
#SBATCH --job-name=Mapping
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=10
#SBATCH -p cosmos
#SBATCH --time=24:00:00

sampleList=$1

inpPath=/userfiles/ualtintas20/Chip-Seq
indPath=/kuacc/users/ualtintas20/genomeAnnotations/hg19/bowtie2HumanIndexes/bowtie2HumanIndexes
outPath=/kuacc/users/ualtintas20/2-Bowtie2/Chip-Seq/samFiles

while IFS="$(printf '\t')" read -r SRX LibraryType; do

	echo $SRX ;

	if [ $LibraryType == "Single" ]
	then
		echo $LibraryType "A";
		bowtie2 --threads 10 \
		-x $indPath \
		-U $inpPath/$SRX"_1.fastq.gz" \
		-S $outPath/$SRX".sam" ;

	elif [ $LibraryType == "Paired" ]
	then
		echo $LibraryType "B";
		bowtie2 --threads 10 \
		-x $indPath \
		-1 $inpPath/$SRX"_1.fastq.gz" -2 $inpPath/$SRX"_2.fastq.gz" \
		-S $outPath/$SRX".sam";
	else
		echo "NO"

	fi


done < <(cat $sampleList)




















##### PROCESSING #### => An optimized and renewed way to process Binary Alignment Mapping files (BAM)



#!/bin/bash
#SBATCH --job-name=Processing
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=15
#SBATCH -p cosmos
#SBATCH --time=24:00:00

sampleList=$1

samPath=/kuacc/users/ualtintas20/2-Bowtie2/Chip-Seq/samFiles
bamPath=/kuacc/users/ualtintas20/2-Bowtie2/Chip-Seq/bamFiles
tmpPath=/kuacc/users/ualtintas20/2-Bowtie2/Chip-Seq/tmpFiles
outPath=/kuacc/users/ualtintas20/2-Bowtie2/Chip-Seq/outFiles

while IFS="$(printf - )" read -r SRX ; do

	echo $SRX

	samtools view -@ 15 -b -o $bamPath/$SRX".bam" $samPath/$SRX".sam";

	mkdir $tmpPath/$SRX


	samtools view -@ 15 -F 4 -b -o $tmpPath/$SRX/$SRX".mapped.bam" $bamPath/$SRX".bam";


	samtools sort -@ 15 -n -o $tmpPath/$SRX/$SRX".mapped.collated.bam" $tmpPath/$SRX/$SRX".mapped.bam" ;
	samtools fixmate -m -O bam $tmpPath/$SRX/$SRX".mapped.collated.bam" $tmpPath/$SRX/$SRX".mapped.collated.fixed.bam" ;
	samtools sort -@ 15 -o $tmpPath/$SRX/$SRX".mapped.collated.fixed.possorted.bam" $tmpPath/$SRX/$SRX".mapped.collated.fixed.bam"
	samtools markdup -r $tmpPath/$SRX/$SRX".mapped.collated.fixed.possorted.bam" $tmpPath/$SRX/$SRX".mapped.collated.fixed.possorted.markdup.bam";

	samtools view -@ 15 -b -q 20 -o $outPath/$SRX".mapped.collated.fixed.possorted.markdup.mapq.bam" $tmpPath/$SRX/$SRX".mapped.collated.fixed.possorted.markdup.bam"

	samtools index -b $outPath/$SRX".mapped.collated.fixed.possorted.markdup.mapq.bam"  $outPath/$SRX".mapped.collated.fixed.possorted.markdup.mapq.bam.bai" ;

	rm -r $tmpPath/$SRX


done < <(cat $sampleList)




#### PEAK CALLING ####


#!/bin/bash
#SBATCH --job-name=PeakCalling
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=10
#SBATCH -p cosmos
#SBATCH --time=24:00:00

sampleList=$1

inpPath=/kuacc/users/ualtintas20/2-Bowtie2/Chip-Seq/outFiles
outPath=/kuacc/users/ualtintas20/5-PeakCall

extStr=".mapped.collated.fixed.possorted.markdup.mapq.bam"


while IFS="$(printf '\t')" read -r tSRX cSRX target cellType condition LibraryType; do

	mkdir $outPath/$tSRX
	echo $LibraryType

	if (( $LibraryType == "Single" ))
	then

		macs2 callpeak -t $inpPath/$tSRX$extStr -c $inpPath/$cSRX$extStr -B -f BAM -g hs --outdir $outPath/$tSRX -n $tSRX"."$target"."$condition"."$cellType -q 0.00001;

	elif (( $LibraryType == "Paired" ))
	then
		echo $LibraryType
		macs2 callpeak -t $inpPath/$tSRX$extStr -c $inpPath/$cSRX$extStr -B -f BAMPE -g hs --outdir $outPath/$tSRX -n $tSRX"."$target"."$condition"."$cellType -q 0.00001;
	fi


done < <(cat $sampleList)







#### KALLISTO ####





#!/bin/bash
#SBATCH --job-name=Kallisto
#SBATCH --mem-per-cpu=10G
#SBATCH --cpus-per-task=15
#SBATCH -p cosmos
#SBATCH --time=24:00:00

sampleList=$1

inpPath=/kuacc/users/ualtintas20/0-Downloading/RNA-Seq
outPath=/kuacc/users/ualtintas20/4-Kallisto/all
indPath=/kuacc/users/ualtintas20/genomeAnnotations/hg19/ensemblTranscriptome/hg19.kallisto.index



while IFS="$(printf '\t')" read -r SRX LibraryType ; do

	mkdir $outPath/$SRX ;



	if [ $LibraryType == "SINGLE" ] ; then

		echo $SRX $LibraryType "A";

		kallisto quant --index=$indPath \
		-b 100 --single -l 180 -s 20 --threads=15 --plaintext \
		--output-dir=$outPath/$SRX \
		$inpPath/singleEnd/$SRX"_1.fastq.gz" ;

	elif [ $LibraryType == "PAIRED" ] ; then

		echo $SRX $LibraryType "B";

		kallisto quant --index=$indPath \
		-b 100 --threads=15 --plaintext \
		--output-dir=$outPath/$SRX \
		$inpPath/pairedEnd/$SRX"_1.fastq.gz" $inpPath/pairedEnd/$SRX"_2.fastq.gz" ;


	fi


done < <(cat $sampleList)
