

def bowtie2_inputs(wildcards):
    if (seq_type == "pe"):
        return expand("{reads}_{strand}.fastq", strand=["R1", "R2"], reads=wildcards.reads)
    elif (seq_type == "se"):
        return expand("{reads}.fastq", reads=wildcards.reads)

def bowtie2_params(wildcards, input):
    if (seq_type == "pe"):
        return f'-1 {input.reads[0]} -2 {input.reads[1]}'
    else:
        return f'-U {input.reads}'

params:
file_args=bowtie2_params

#



idx = config["reference"]["idx"]


rule mapping:
        input:
                "rawData/{raw}.fastq.gz"
        output:
                "results/mapping/{raw}.bam"
        log:
                "logs/log.mapping.{raw}.txt"
        message:
                "Executing mappingSingle rule for {wildcards.raw}"
        threads:
                64
        shell:
                """
                echo "\n mapping {input}"

                gzip -dc {input} | \
                bowtie --chunkmbs 512 -k 1 -m 1 -v 2 --best --strata {idx} --threads {threads} -q - -S| \
                samtools view -bS - > {output}
                """
