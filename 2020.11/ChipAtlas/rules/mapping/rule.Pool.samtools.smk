

def getRepsToPool(wildcards):
        """
        Used by:
                -> Pool
        """
        if wildcards.sampleName.find("control") > 0:
                sampleName = wildcards.sampleName.rsplit(".", 1)[0]
                #print("ibne", sampleName, wildcards)
                reps = list(controlDF.loc[controlDF["ControlName"] == sampleName, "Replicate"])
                return expand("results/mapping/final/{sampleName}.{rep}.control.bam", rep=reps, sampleName=[sampleName])
        else:
                reps = list(sampleDF.loc[sampleDF["SampleName"] == wildcards.sampleName, "Replicate"])
                return expand("results/mapping/final/{{sampleName}}.{rep}.bam", rep=reps)



rule Pool:
        input:
                getRepsToPool
        output:
                bam="results/mapping/merged/{sampleName}.bam",
                bai="results/mapping/merged/{sampleName}.bam.bai"
        message:
                "Executing pool rule for {wildcards}"
        threads:
                16
        shell:
                """
                samtools merge -O BAM -@ {threads} {output.bam} {input}
                samtools index  -@ {threads} {output.bam} {output.bai}
                """
