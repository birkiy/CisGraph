



def getRepsToConcat(wildcards):
        lib = srxDF.loc[srxDF["SRX"] == wildcards.SRX, "Library"].to_list()[0]
        SRR = srxDF.loc[srxDF["SRX"] == wildcards.SRX, "SRR"].to_list()

        if wildcards.run == "U" or wildcards.run == "R1":
                return expand("raw/{srr}_1.fastq.gz", srr=SRR)
        elif wildcards.run == "R2":
                return expand("raw/{srr}_2.fastq.gz", srr=SRR)



rule ConcatFastq:
        input:
                getRepsToConcat
        output:
                "raw/{SRX}.{run}.fastq.gz"
        message:
                "Executing ConcatFastq rule for {wildcards}"
        threads:
                4
        shell:
                """
                zcat {input} > {output}
                """
