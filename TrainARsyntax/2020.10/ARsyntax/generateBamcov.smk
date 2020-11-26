
def getParamBamCov(wildcards):
    if wildcards.type == "unique":
        return ""
    elif wildcards.type in ["RPKM", "None", "BPM"]:
        return "--normalizeUsing {wildcards.type}"


rule bamCov:
        input:
                bam="results/mapping/processed/{raw}.final.bam"
        output:
                out0="results/bigwig/{raw}.extended.{type}.bigWig",
                out1="results/bigwig/{raw}.extended.forward.{type}.bigWig",
                out2="results/bigwig/{raw}.extended.reverse.{type}.bigWig"
        threads:
                16
        message:
                "Executing bamCov rule for {wildcards.raw}"
        params:
                getParamBamCov
        shell:
                """
                bamCoverage --bam {input.bam} -o {output.out0} \
                    --samFlagInclude 64 \
                    --extendReads 150 \
                    --blackListFileName {blacklistFile} \
                    -p {threads}

                bamCoverage --bam {input.bam} -o {output.out1} \
                    --filterRNAstrand forward \
                    --samFlagInclude 64 \
                    --extendReads 150 \
                    --blackListFileName {blacklistFile} \
                    -p {threads}

                bamCoverage --bam {input.bam} -o {output.out2} \
                    --filterRNAstrand reverse \
                    --samFlagInclude 64 \
                    --extendReads 150 \
                    --blackListFileName {blacklistFile} \
                    -p {threads}
                """
