# BamProcess rule

rule BamProcess:
        input:
                bam="results/mapping/raw/{raw}.{rep}.bam"

        output:
                bam="results/mapping/final/{raw}.{rep}.bam",
                bai="results/mapping/final/{raw}.{rep}.bam.bai",
                nameSorted="results/mapping/namesorted/{raw}.{rep}.bam"
        threads:
                16
        message:
                "Executing BamProcess rule for {wildcards}"
        params:
                "" # Note that you can add parameters as "-q 30 -F 1804"
        shell:
                """
                N=$((`samtools view {input.bam} | wc -l `))
                echo "Number of reads before bam filtration:"$N

                cat <(samtools view -H {input.bam}) <(samtools view {params} {input.bam}) | \
                samtools fixmate -r -m -@ {threads} - - | \
                samtools sort -@ {threads} -m 10G - | \
                samtools markdup - - | \
                samtools view -b -r - > {output.bam}

                samtools index  -@ {threads} {output.bam} {output.bai}

                samtools sort -@ {threads} -m 10G -n {output.bam} -o {output.nameSorted}

                N=$((`samtools view {output.bam} | wc -l `))
                echo "Number of reads after bam filtration:"$N
                """
