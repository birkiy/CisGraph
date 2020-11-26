


blacklistFile = config["blacklist"]

def getParamBamCov(wildcards):
        print(wildcards.type)
        if wildcards.type == "unique":
                return ""
        elif wildcards.type in ["RPKM", "None", "BPM"]:
                return f"--normalizeUsing {wildcards.type}"



rule bamCov:
        input:
                bam="results/mapping/processed/{raw}.final.bam"
        output:
                out0="results/bigwig/{raw}.{type}.noComp.bigWig"
        threads:
                16
        message:
                "Executing bamCov rule for {wildcards.raw}"
        params:
                getParamBamCov
        shell:
                """
                echo "hey"
                echo {params}

                bamCoverage --bam {input.bam} -o {output.out0} \
                --extendReads 150 \
                --blackListFileName {blacklistFile} \
                -p {threads} \
                {params}
                """



def pools(wildcards):
        return checkpoints.pool.get(**wildcards).output.bam


rule bamCom:
        input:
                #bam=pools,
                bam="results/mapping/processed/{raw}.final.bam",
                inp="results/mapping/processed/input.GSE114063.merged.final.bam"
        output:
                out0="results/bigwig/{raw}.{type}.SES.bigWig",

        threads:
                16
        message:
                "Executing bamCom rule for {wildcards.raw}"
        params:
                getParamBamCov
        shell:
                """
                echo "hey"
                echo {params}

                bamCompare -b1 {input.bam} -b2 {input.inp} -o {output.out0} \
                --blackListFileName {blacklistFile} \
                --extendReads 150 \
                --scaleFactorsMethod SES \
                --operation ratio \
                -p {threads}\
                --scaleFactors 1.2:0.5 \
                {params}
                """
