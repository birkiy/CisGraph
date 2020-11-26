

# def getRepsToPool(wildcards):
#         return [
#                 f"results/mapping/processed/{raw}.{rep}.final.bam"
#                 for rep in list(sampleDF.loc[sampleDF["sampleName"] == wildcards.raw, "rep"])
#                 for raw in [wildcards.raw]
#                 ]



def getRepsToPool(wildcards):
        print("hey")
        reps = list(sampleDF.loc[sampleDF["sampleName"] == wildcards.raw, "rep"])
        print(reps)
        return expand("results/mapping/processed/{{raw}}.{rep}.final.bam", rep=reps)


rule pool:
        input:
                getRepsToPool
        output:
                bam="results/mapping/processed/{raw}.merged.final.bam",
                bai="results/mapping/processed/{raw}.merged.final.bam.bai"
        message:
                "Executing pool rule for {wildcards.raw}"
        threads:
                16
        shell:
                """
                #Merge treatment BAMS
                samtools merge -@ {threads} -u {output} {input}

                samtools index  -@ {threads} {output.bam} {output.bai}
                """
