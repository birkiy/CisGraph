




bedtools genomecov -d -ibam results/mapping/SRX4108929.1.control.final.bam -scale `bc -l <<< 1000000/1613426` | bedtools sort > deneme.genomcov.bedgraph



blacklistFile = config["hg19"]["blacklist"]

rule BaseDepth:
        input:
                "results/mapping/{sampleName}.merged.final.bam"
        output:
                "results/bigWig/{type}/{sampleName}.bedGraph"
        shell:
                """
                bedtools genomecov -bg -strand + -ibam {input.bam} -scale `bc -l <<< 1000000/$N` | bedtools sort > {output.BG}
