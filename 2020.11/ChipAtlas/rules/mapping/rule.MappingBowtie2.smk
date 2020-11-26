
idx = config["hg19"]["idx"]["bowtie2"]


def getFastqs(wildcards):
        """
        Used by:
                -> MappingBowtie2
        """
        outputD = dict()
        lib = getLib(wildcards)
        if lib == "Single":
                outputD["U"] = f"links/{wildcards.raw}.U.fastq.gz"
        elif lib == "Paired":
                outputD["R1"] = f"links/{wildcards.raw}.R1.fastq.gz"
                outputD["R2"] = f"links/{wildcards.raw}.R2.fastq.gz"
        return outputD


rule MappingBowtie2:
        input:
                unpack(getFastqs)
        output:
                bam="results/mapping/raw/{raw}.bam"
        message:
                "Executing MappingBowtie2 rule with bowtie2 for {wildcards.raw}"
        threads:
                16
        run:
                lib = getLib(wildcards)

                if lib == "Single":
                        shell("""
                        bowtie2 -x {idx} --threads {threads} -q {input.U} {params} | \
                        samtools view -bS - > {output.bam}
                        """)
                elif lib == "Paired":
                        shell("""
                        bowtie2 -x {idx} --threads {threads} -q -1 {input.R1} -2 {input.R2} {params} | \
                        samtools view -bS - > {output.bam}
                        """)
