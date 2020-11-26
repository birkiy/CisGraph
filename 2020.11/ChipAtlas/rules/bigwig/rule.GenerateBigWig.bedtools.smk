


chrSize = config["hg19"]["chrSize"]

blacklistFile = config["hg19"]["blacklist"]

rule GenerateBigWig:
        input:
                bam="results/mapping/{stage}/{raw}.bam"
        output:
                BG=temp("results/bigwig/{stage}/{raw}.bedGraph"),
                BGblk=temp("results/bigwig/{stage}/{raw}.blk.bedGraph"),
                BW="results/bigwig/genomcov/{stage}/{raw}.bigWig",
        threads:
                4
        message:
                "Executing GenerateBigWig rule with bedtools for {wildcards}."
        shell:
                """
                N=$((`samtools view {input.bam} | wc -l `))
                echo "Number of reads for coverage:"$N
                bedtools genomecov -bg -ibam {input.bam} -scale `bc -l <<< 1000000/$N` | bedtools sort > {output.BG}

                echo "Number of regions covered:"`wc -l {output.BG}`
                bedtools intersect -v -a {output.BG} -b {blacklistFile} > {output.BGblk}

                echo "Number of regions w/o blacklist:"`wc -l {output.BGblk}`
                bedGraphToBigWig {output.BGblk} {chrSize} {output.BW}
                """
