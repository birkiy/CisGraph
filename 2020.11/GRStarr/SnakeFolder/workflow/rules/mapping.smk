


idx = config["reference"]["idx"]



rule mapping:
        input:
                forwardR="rawData/{raw}.R1.fastq.gz",
                reverseR="rawData/{raw}.R2.fastq.gz"
        output:
                sam=temp("results/mapping/{raw}.sam"),
                bam="results/mapping/{raw}.bam"
        message:
                "Executing mappingSingle rule for {wildcards.raw}"
        threads:
                16
        shell:
                """
                bwa mem -t {threads} {idx} {input.forwardR} {input.reverseR} > {output.sam}
                samtools sort -m 10G -@ {threads}  {output.sam}  > {output.bam}
                """
