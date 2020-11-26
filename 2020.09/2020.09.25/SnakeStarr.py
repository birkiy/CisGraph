

reference_genome="/home/ualtintas/genomeAnnotations/hg19.fa"

bwa_dict="/home/ualtintas/genomeAnnotations/hg19.bwa.idx/hg19.bwa.idx"
blacklist_file="/home/ualtintas/genomeAnnotations/ENCFF001TDO.bed"
gr_region_path="/home/ualtintas/GR/GR.bed"

analysis_section={"bigwig":["unique-reads.bigWig","raw-extended-reads.bigWig","norm-extended-reads.bigWig"],
"coverage":['extended-reads.count-table-deeptools.txt','extended-reads.count-table-deeptools.npz']}

samples = ["SRR7122167", "SRR7122172", "SRR7122177", "SRR7122182", "SRR7122187", "SRR7122192", "SRR7122193", "SRR7122194", "SRR7122195", "SRR7122197", "SRR7122198", "SRR7122196", "SRR7122199", "SRR7122200", "SRR7122201", "SRR7122202", "SRR7122203"]



def get_fastq():
        for sample_name in samples:
                for name,outputs in analysis_section.items():
                        for i in outputs:
                                #print("analysis/{}/{}-{}.{}".format(list(j)[0],i,list(j)[1],k))
                                yield f"analysis/{name}/{sample_name}-{i}"


rule all:
        input: get_fastq()

rule mapping:
        input:
                forward="rawFasta/{sample_name}_1.fastq.gz",
                reverse="rawFasta/{sample_name}_2.fastq.gz"
        output:
                bam="analysis/mapping/{sample_name}.bam",
                bam_index="analysis/mapping/{sample_name}.bam.bai"
        log: "logs/log-mapping-{sample_name}.txt"
        threads: 40
        shell:
                """
                bwa mem -t {threads} {bwa_dict} {input.forward} {input.reverse} 2>{log} | samtools sort -m 15G -@ {threads} - > {output.bam}
                samtools index -@ {threads} {output.bam} {output.bam_index}
                """
rule bam_filtration:
        input:
                rules.mapping.output.bam
        output:
                bam="analysis/mapping/{sample_name}-cleaned.bam",
                bam_index="analysis/mapping/{sample_name}-cleaned.bam.bai"
        threads:16
        shell:
                """
                cat <(samtools view -H {input}) <(samtools view -q 60 {input} | awk '$6 !~ "I|D"' - ) |
                samtools view -Sb - | samtools sort -m 10G -@ {threads} -n - | samtools fixmate -m -@ {threads} - - |
                samtools sort -m 10G -@ {threads} - | samtools markdup -r - - | samtools view -f2 -b - > {output.bam}
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
        threads: 16
        shell:
                """
                multiBamSummary BED-file --BED {gr_region_path} --bamfiles {input} --extendReads --samFlagInclude 64 --blackListFileName {blacklist_file} -p {threads} --outRawCounts {output.extended_reads_count_table} --labels {wildcards.sample_name} -out {output.temperory_numpy} 2>{log.log4}
                """



#
# bwa mem -t 15 /home/ualtintas/genomeAnnotations/hg19.bwa.idx/hg19.bwa.idx /home/ualtintas/GR/STARRbegin/rawFasta/SRR7122167_1.fastq.gz /home/ualtintas/GR/STARRbegin/rawFasta/SRR7122167_2.fastq.gz 2>logs/log-mapping-SRR7122167.txt | samtools sort -m 5G -@ 15 - > analysis/mapping/SRR7122167.bam
