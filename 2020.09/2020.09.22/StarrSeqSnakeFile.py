

reference_genome="/home/tmorova/tools/hg19.fa"
bwa_dict="/home/tmorova/tools/bwa_hg19_index/hg19.fa"
blacklist_file="/home/tmorova/blacklisted_regions/hg19-blacklist.bed"
validation_region_path="/home/tmorova/STARR/coverage/ValidationRegions_ARBS-NControl-PControl.bed"
arbs_region_path="/home/tmorova/STARR/coverage/ARBS.bed"
pos_control_region_path="/home/tmorova/STARR/coverage/PosCont.bed"
neg_control_region_path="/home/tmorova/STARR/coverage/NegControl.bed"

analysis_section={"bigwig":["unique-reads.bigWig","raw-extended-reads.bigWig","norm-extended-reads.bigWig"],
"coverage":['extended-reads.count-table-deeptools.txt','extended-reads.count-table-deeptools.npz']}

#analysis_section_output_extention=["bigWig","count-table-deeptools.txt"]
#output_type=["unique-reads","raw-extended-reads","norm-extended-reads"]



def get_fastq():
	for sample_name in config["samples"]:
		for name,outputs in analysis_section.items():
			for i in outputs:
				#print("analysis/{}/{}-{}.{}".format(list(j)[0],i,list(j)[1],k))
				yield "analysis/{}/{}-{}".format(name,sample_name,i)

get_fastq()
rule all:
	input: get_fastq()

rule mapping:
	input:
		forward="raw-data/starrseq/{sample_name}_R1_001.fastq.gz",
		reverse="raw-data/starrseq/{sample_name}_R2_001.fastq.gz"
	output:
		bam=temp("analysis/mapping/{sample_name}.bam"),
		bam_index=temp("analysis/mapping/{sample_name}.bam.bai")
	log: "logs/log-mapping-{sample_name}.txt"
	conda:
		"code/starrseq-sm.env"
	threads: 32
	shell:
		"""
		bwa mem -t {threads} {bwa_dict} {input.forward} {input.reverse} 2>{log} | samtools sort -m 5G -@ {threads} - > {output.bam}
		samtools index -@ {threads} {output.bam} {output.bam_index}
		"""
rule bam_filtration:
	input:
		rules.mapping.output.bam
	output:
		bam="analysis/mapping/{sample_name}-cleaned.bam",
		bam_index="analysis/mapping/{sample_name}-cleaned.bam.bai"
	conda:
		"code/starrseq-sm.env"
	threads:16
	shell:
		"""
		cat <(samtools view -H {input}) <(samtools view -q 60 {input} | awk '$6 !~ "I|D"' - ) |
		samtools view -Sb - | samtools sort -@ {threads} -n - | samtools fixmate -@ {threads} - - |
		samtools sort - | samtools view -f2 -b - > {output.bam}
		samtools index -@ {threads} {output.bam} {output.bam_index}
		"""
rule generate_bigwig:
	input:
		bam=rules.bam_filtration.output.bam
	output:
		unique_bigwig="analysis/bigwig/{sample_name}-unique-reads.bigWig",
		extended_bigwig_raw="analysis/bigwig/{sample_name}-raw-extended-reads.bigWig",
		extended_bigwig_norm="analysis/bigwig/{sample_name}-norm-extended-reads.bigWig"
	log:
		log1="logs/log-generate_bigwig-{sample_name}-unique-reads.txt",
		log2="logs/log-generate_bigwig-{sample_name}-extended-reads-raw.txt",
		log3="logs/log-generate_bigwig-{sample_name}-extended-reads-norm.txt"
	conda:
		"code/starrseq-sm.env"
	threads:16
	shell:
		"""
		bamCoverage --bam {input.bam} -o {output.unique_bigwig} --samFlagInclude 64 --extendReads --ignoreDuplicates --blackListFileName {blacklist_file} -p {threads} 2>{log.log1}
		bamCoverage --bam {input.bam} -o {output.extended_bigwig_raw} --samFlagInclude 64 --extendReads --normalizeUsing None --blackListFileName {blacklist_file} -p {threads} 2>{log.log2}
		bamCoverage --bam {input.bam} -o {output.extended_bigwig_norm} --samFlagInclude 64 --extendReads --normalizeUsing RPKM --blackListFileName {blacklist_file} -p {threads} 2>{log.log3}
		"""
rule deeptools_count_overlapping_reads:
	input:
		bam=rules.bam_filtration.output.bam
	output:
		extended_reads_count_table="analysis/coverage/{sample_name}-extended-reads.count-table-deeptools.txt",
		temperory_numpy=temp("analysis/coverage/{sample_name}-extended-reads.count-table-deeptools.npz")
	log:
		log4="logs/log-deeptools_count_overlapping_reads-{sample_name}-extended-reads.txt"
	conda:
		"code/starrseq-sm.env"
	threads: 16
	shell:
		"""
		multiBamSummary BED-file --BED {validation_region_path} --bamfiles {input} --extendReads --samFlagInclude 64 --blackListFileName {blacklist_file} -p {threads} --outRawCounts {output.extended_reads_count_table} --labels {wildcards.sample_name} -out {output.temperory_numpy} 2>{log.log4}
		"""
