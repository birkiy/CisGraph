
chrSize = config["reference"]["chrSize"]

rule generateBigwig:
        input:
                # rules.mapping.bamProcess.output
                "results/mapping/processed/{raw}.final.bam"
        output:
                forwardBG=temp("results/bigwig/{raw}.+.5end.bedGraph"),
                forwardBW="results/bigwig/{raw}.+.5end.bigWig",
                reverseBG=temp("results/bigwig/{raw}.-.5end.bedGraph"),
                reverseBW="results/bigwig/{raw}.-.5end.bigWig"

        log:
                log1="logs/log.generateBigwig.{wildcards.raw}.txt"
        threads:
                16
        message:
                "Executing generateBigwig rule for {wildcards.raw}"
        shell:
                """
                echo "BedTools Start!"
                bedtools genomecov -5 -bg -strand + -ibam {input} | bedtools sort > {output.forwardBG}
                wc -l {output.forwardBG}
                bedGraphToBigWig {output.forwardBG} {chrSize} {output.forwardBW}

                bedtools genomecov -5 -bg -strand - -ibam {input} | bedtools sort > {output.reverseBG}
                wc -l {output.reverseBG}
                bedGraphToBigWig {output.reverseBG} {chrSize} {output.reverseBW}

                """
#
