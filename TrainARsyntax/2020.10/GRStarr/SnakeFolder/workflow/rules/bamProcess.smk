


rule bamProcess:
        input:
                "results/mapping/{raw}.bam"
        output:
                bam="results/mapping/processed/{raw}.final.bam",
                bai="results/mapping/processed/{raw}.final.bam.bai"

        threads:
                16
        message:
                "Executing bamProcess rule for {wildcards.raw}"
        shell:
                """
                echo "\n bamProcess {input} \n"

                samtools index -@ {threads} {input} {output.bai}

                cat <(samtools view -H {input}) <(samtools view -q 60 -F 1804 {input} | awk '$6 !~ "I|D"' - ) | \
                samtools sort -@ {threads} -m 10G -n - | \
                samtools fixmate -m -@ {threads} - - | \
                samtools sort -@ {threads} -m 10G - | \
                samtools markdup -r - - | \
                samtools view -b - > {output.bam}

                samtools index  -@ {threads} {output.bam} {output.bai}

                samtools view -b  {output.bam} | wc -l
                """
