

TFs = ["AR", "FOXA1"]

conditions = ["dht", "dmso"]

directions = ["forward", "reverse"]

rep

rule





[ARpeaks, FOXA1peaks]


#
#
# def get_fastq():
#     for condition in conditions:
#         for tf in TFs:
#             for name,outputs in analysis_section.items():
#                 for i in outputs:
#                     yield f"analysis/{name}/{sample_name}.{i}"
#
#
#
#
# def getBam():
#     for condition in conditions:
#         for tf in TFs:
#             yield {
#                 "rep1":
#                     f"analysis/mapping/{tf}.{condition}.rep1.final.bam",
#                 "rep2":
#                     f"analysis/mapping/{tf}.{condition}.rep2.final.bam"
#                 }

rule all:
        input:
                expand("analysis/peaks/{tf}.{condition}.{ext}",
                        tf=["AR", "FOXA1"],
                        condition=["dht", "dmso"],
                        rep=["rep1", "rep2"],
                        ext=["_peaks.xls", "_peaks.narrowPeak"])



# rule mergeRep:
#         input:
#                 rep1="analysis/bigwig/{sample_name}.{condition}.rep1.forward.5end.bigWig",
#                 rep2="analysis/bigwig/{sample_name}.{condition}.rep2.forward.5end.bigWig",
#         output:
#                 bigwig="analysis/bigwig/{sample_name}.{condition}.merged.forward.5end.bigWig"
#
#         shell:
#                 """
#                 bigwigCompare -b1 {input.rep1} -b2 {input.rep2} --ratio=add  -o {}
#                 """

rule generatePseudoReplicates:
        input:
                getBam()
                # rep1="analysis/mapping/{tf}.{condition}.rep1.bam"
                # rep2="analysis/mapping/{tf}.{condition}.rep2.bam"
        output:
                bam=temp("analysis/mapping/{tf}.{condition}.merged.bam"),
                header=temp("analysis/mapping/{tf}.{condition}.merged.header.sam"),
                prBam0="analysis/mapping/{tf}.{condition}.pr0.bam",
                prBam1="analysis/mapping/{tf}.{condition}.pr1.bam"
        message: "Executing generatePseudoReplicates rule for {wildcards.tf} {wildcards.condition}"
        threads: 16
        shell:
                """
                samtools merge -u {output.bam} {input.rep1} {input.rep2} -@ {threads}

                samtools view -H {output.bam} > {output.header}

                nlines=$(samtools view {output.bam} | wc -l )
                nlines=$(( (nlines + 1) / 2 )) # half that number

                samtools view {output.bam} | shuf - | split -d -l ${nlines} - "{wildcards.tf}.{wildcards.condition}"

                cat header {wildcards.tf}.{wildcards.condition}00 | samtools view -bS - > {output.prBam0}
                cat header {wildcards.tf}.{wildcards.condition}01 | samtools view -bS - > {output.prBam1}

                """







rule generateBigwig:
        input:
                rules.generatePseudoReplicates.output.prBam0

        output:
                # TODO::Try to do expand like function. Is it possible?
                forwardBG=temp("analysis/bigwig/{tf}.{condition}.merged.forward.5end.bedGraph"),
                reverseBG=temp("analysis/bigwig/{tf}.{condition}.merged.reverse.5end.bedGraph"),
                forward5="analysis/bigwig/{tf}.{condition}.merged.forward.5end.bigWig",
                reverse5="analysis/bigwig/{tf}.{condition}.merged.reverse.5end.bigWig"
        message: "Executing generateBigwig rule for {tf}.{condition}"
        threads:16
        shell:
                """
                bedtools genomecov -5 -bg -strand + -ibam {input} | bedtools sort > {output.forwardBG}
                bedGraphToBigWig {output.forwardBG} {chrSizes} {output.forward5}

                bedtools genomecov -5 -bg -strand - -ibam {input} | bedtools sort > {output.reverseBG}
                bedGraphToBigWig {output.reverseBG} {chrSizes} {output.reverse5}
                """

rule PeakCalling:
        input:
                chip="analysis/mapping/{tf}.{condition}.{rep}.bam",
                control="analysis/mapping/input.{rep}.bam"
        output:
                getOut()
        message: "Executing peakCalling rule for {tf}.{condition}"
        threads:16
        shell:
                """
                macs2 callpeak \
                    -t {input.chip} \
                    -c {input.control} \
                    -n {tf}.{condition}.{rep} \
                    -B -f BAM -g hs -q 0.001 \
                    --outdir analysis/peaks \
                    --shift -75 --extsize 150

                """


samtools merge -r analysis/mapping/input.merged.bam analysis/mapping/input.rep1.final.bam analysis/mapping/input.rep2.final.bam -@ 10

rule macs2:
    input: trt = f'{testDir}'+'/{sample}.bam', inp = expand('{inpfile}',inpfile=inplist)
    output: xls = f'{pcDir}' + '/{sample}_peaks.xls' ,pk = f'{pcDir}' + '/{sample}_peaks.narrowPeak'
    params: nm = '{sample}' ,
        pcout = pcDir ,
        gs = "hs" ,
        th = 0.05
    threads: 40
    conda: "macs2.yaml"
    log: f'{logDir}' + '/macs2_{sample}.log'
    shell: "macs2 callpeak -t {input.trt} -c {input.inp} -f BAM -g {params.gs} -q {params.th} --outdir {params.pcout} -n {params.nm} 2> {log}"
