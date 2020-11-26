
blacklistFile = config["hg19"]["blacklist"]




rule GenerateBigWigBamCom:
        input:
                bam="results/mapping/{stage}/{raw}.bam",
                control=getControl
        output:
                BW="results/bigwig/{type}/bamcom/{stage}/{raw}.bigWig",

        threads:
                16
        message:
                "Executing GenerateBigWigBamCom rule for {wildcards.raw} with type of {wildcards.type} at {wildcards.stage} stage."
        params:
                getParamsBamCom
        shell:
                """
                bamCompare -b1 {input.bam} -b2 {input.control} -o {output.BW} \
                --extendReads 150 \
                --centerReads \
                -p {threads} \
                {params}
                """
