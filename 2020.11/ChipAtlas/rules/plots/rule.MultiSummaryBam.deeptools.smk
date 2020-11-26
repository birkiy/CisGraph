

blacklistFile = config["hg19"]["blacklist"]


def getAll(wildcards):
        if wildcards.samples == "merged":
                samples = sorted(set((sampleDF["SampleName"]).to_list() + (controlDF["ControlName"] + ".control").to_list()))
        elif wildcards.samples == "final":
                samples = sorted(set(list(sampleDF["Raw"]) + list(controlDF["Raw"] + ".control")))
        return expand("results/mapping/{{samples}}/{sample}.bam", sample=samples)


rule MultiBamSummary:
        input:
                samples=getAll,
        output:
                rawCounts="results/coverage/countTableRaw.{samples}.txt",
                matrix="results/coverage/count-table-deeptools.{samples}.npz"
        threads:
                32
        message:
                "Executing MultiBamSummary rule"
        shell:
                """
                multiBamSummary bins \
                --bamfiles {input.samples} \
                --extendReads 150 \
                --centerReads \
                -p {threads} \
                --outRawCounts {output.rawCounts} -out {output.matrix}
                """
