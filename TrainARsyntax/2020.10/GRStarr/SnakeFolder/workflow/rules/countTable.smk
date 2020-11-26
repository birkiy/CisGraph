




bedGR = config["bedGR"]

blacklistFile = config["blacklist"]

def getAll():
        samples = sampleDF[["sampleName", "rep"]].agg(".".join, axis=1)
        return expand("results/mapping/processed/{sample}.final.bam", sample=samples)


rule countTable:
        input:
                getAll()
        output:
                rawCounts="results/coverage/countTableRaw.txt",
                matrix=temp("results/coverage/count-table-deeptools.npz")
        threads:
                32
        message:
                "Executing countMatrix rule"
        shell:
                """
                multiBamSummary BED-file \
                --BED {bedGR} \
                --bamfiles {input} \
                --extendReads --samFlagInclude 64 \
                --blackListFileName {blacklistFile} \
                -p {threads} \
                --outRawCounts {output.rawCounts} -out {output.matrix}
                """
