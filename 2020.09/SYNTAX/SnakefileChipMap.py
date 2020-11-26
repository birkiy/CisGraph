

reference_genome="/home/ualtintas/genomeAnnotations/hg19.fa"

idx="/home/ualtintas/genomeAnnotations/bowtieIdx/hg19"
blacklistFile="/home/ualtintas/genomeAnnotations/ENCFF001TDO.bed"

chrSizes="/home/ualtintas/genomeAnnotations/hg19.chrom.sizes"

analysis_section={
    "bigwig":[
        "forward.5end.bigWig",
        "reverse.5end.bigWig",
        "extended.reverse.raw.bigWig",
        "extended.reverse.raw.bigWig"
        ]
    }

samples = [
"AR.dht.rep1",
"AR.dht.rep2",
"AR.dmso.rep1",
"AR.dmso.rep2",
"ARID.r1881.rep1",
"ARID.r1881.rep2",
"ARID.veh.rep1",
"ARID.veh.rep2",
"BRG1.r1881.rep1",
"BRG1.r1881.rep2",
"BRG1.veh.rep1",
"BRG1.veh.rep2",
"FOXA1.dht.rep1",
"FOXA1.dht.rep2",
"FOXA1.dmso.rep1",
"FOXA1.dmso.rep2",
"HOB13.r1881.rep1",
"HOB13.r1881.rep2",
"HOB13.r1881.rep3",
"HOB13.veh.rep1",
"HOB13.veh.rep2",
"HOB13.veh.rep3",
"input.rep1",
"input.rep2",
"TLE3.r1881.rep1",
"TLE3.r1881.rep2",
"TLE3.veh.rep1",
"TLE3.veh.rep2",
"TRIM28.r1881.rep1",
"TRIM28.r1881.rep2",
"TRIM28.veh.rep1",
"TRIM28.veh.rep2",
"WDHD1.r1881.rep1",
"WDHD1.r1881.rep2",
"WDHD1.veh.rep1",
"WDHD1.veh.rep2"
]

# sample_name="AR.dht.rep1"
def get_fastq():
        for sample_name in samples:
                for name,outputs in analysis_section.items():
                        for i in outputs:
                                #print("analysis/{}/{}-{}.{}".format(list(j)[0],i,list(j)[1],k))
                                yield f"analysis/{name}/{sample_name}.{i}"


rule all:
        input: get_fastq()
        # input: "analysis/mapping/AR.dht.rep1.final.bam"

rule mapping:
        input: "rawData/{sample_name}.fastq.gz",
        output:
                sam="analysis/mapping/{sample_name}.sam"
        log: "logs/log.mapping.{sample_name}.txt"
        threads: 64
        shell:
                """
                echo "\n mapping {input}"

                gzip -dc {input} | bowtie --chunkmbs 512 -k 1 -m 1 -v 2 --best --strata {idx} --threads {threads} -q - -S {output.sam}
                """
rule bamProcess:
        input:
                rules.mapping.output.sam
                # "analysis/mapping/{sample_name}.sam"
        output:
                bam="analysis/mapping/{sample_name}.final.bam",
                bamIndex="analysis/mapping/{sample_name}.final.bam.bai"
        threads:16
        shell:
                """
                echo "\n bamProcess {input} \n"

                cat <(samtools view -H {input}) <(samtools view -q 30 -F 1804 {input} | awk '$6 !~ "I|D"' - ) | \
                samtools sort -@ {threads} -m 10G -n - | \
                samtools fixmate -m -@ {threads} - - | \
                samtools sort -@ {threads} -m 10G - | \
                samtools markdup -r - - | \
                samtools view -b - > {output.bam}
            	samtools index -@ {threads} {output.bam} {output.bamIndex}

                samtools view -b  {output.bam} | wc -l
                """
rule generateBigwig:
        input:
                rules.bamProcess.output.bam
        output:
                forwardBG="analysis/bigwig/{sample_name}.forward.5end.bedGraph",
                reverseBG="analysis/bigwig/{sample_name}.reverse.5end.bedGraph",
                forward5="analysis/bigwig/{sample_name}.forward.5end.bigWig",
                reverse5="analysis/bigwig/{sample_name}.reverse.5end.bigWig",
                extendedForwardRaw="analysis/bigwig/{sample_name}.extended.forward.raw.bigWig",
                extendedReverseRaw="analysis/bigwig/{sample_name}.extended.reverse.raw.bigWig"
        log:
                log1="logs/log.generateBigwig.{sample_name}.forward.txt",
                log2="logs/log.generateBigwig.{sample_name}.forward.txt"
        threads:16
        shell:
                """
                echo "\n generateBigwig {input}"

                echo "BedTools Start!"
                bedtools genomecov -5 -bg -strand + -ibam {input} | bedtools sort > {output.forwardBG}
                wc -l {output.forwardBG}
                bedGraphToBigWig {output.forwardBG} {chrSizes} {output.forward5}

                bedtools genomecov -5 -bg -strand - -ibam {input} | bedtools sort > {output.reverseBG}
                bedGraphToBigWig {output.reverseBG} {chrSizes} {output.reverse5}

                echo "DeepTools Start! \n"

                bamCoverage --bam {input} -o {output.extendedForwardRaw} \
                    --filterRNAstrand forward \
                    --samFlagInclude 64 \
                    --extendReads 150 \
                    --normalizeUsing None \
                    --blackListFileName {blacklistFile} \
                    -p {threads} 2>{log.log1}

                bamCoverage --bam {input} -o {output.extendedReverseRaw} \
                    --filterRNAstrand reverse \
                    --samFlagInclude 64 \
                    --extendReads 150 \
                    --normalizeUsing None \
                    --blackListFileName {blacklistFile} \
                    -p {threads} 2>{log.log2}
                """
#


# samtools view {output.sam} -@ {threads} -Sb | samtools sort -m 10G -@ {threads} - > {output.bam}
#
# samtools index -@ {threads} {output.bam} {output.bamIndex}


# rule callPeaks:
#         input:
#                 bam = rules.generateBigwig.output.bam
#                 inputBias = analysis/mapping/{sample_name}.final.bam
#
#
#
# rule deeptools_count_overlapping_reads:
#         input:
#                 bam=rules.bam_filtration.output.bam
#         output:
#                 extended_reads_count_table="analysis/coverage/{sample_name}-extended-reads.count-table-deeptools.txt",
#                 temperory_numpy=temp("analysis/coverage/{sample_name}-extended-reads.count-table-deeptools.npz")
#         log:
#                 log4="logs/log-deeptools_count_overlapping_reads-{sample_name}-extended-reads.txt"
#         threads: 16
#         shell:
#                 """
#                 multiBamSummary BED-file --BED {gr_region_path} --bamfiles {input} --extendReads --samFlagInclude 64 --blackListFileName {blacklist_file} -p {threads} --outRawCounts {output.extended_reads_count_table} --labels {wildcards.sample_name} -out {output.temperory_numpy} 2>{log.log4}
#                 """
#


#
# bwa mem -t 15 /home/ualtintas/genomeAnnotations/hg19.bwa.idx/hg19.bwa.idx /home/ualtintas/GR/STARRbegin/rawFasta/SRR7122167_1.fastq.gz /home/ualtintas/GR/STARRbegin/rawFasta/SRR7122167_2.fastq.gz 2>logs/log-mapping-SRR7122167.txt | samtools sort -m 5G -@ 15 - > analysis/mapping/SRR7122167.bam
