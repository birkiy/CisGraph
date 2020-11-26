




rule bamProcess:
        input:
                rules.mapping.output
        output:
                "results/mapping/processed/{raw}.final.bam"
        threads:
                16
        message:
                "Executing bamProcess rule for {wildcards.raw}"
        shell:
                """
                echo "\n bamProcess {input} \n"

                cat <(samtools view -H {input}) <(samtools view -q 30 -F 1804 {input} | awk '$6 !~ "I|D"' - ) | \
                samtools sort -@ {threads} -m 10G -n - | \
                samtools fixmate -m -@ {threads} - - | \
                samtools sort -@ {threads} -m 10G - | \
                samtools markdup -r - - | \
                samtools view -b - > {output}


                samtools view -b  {output} | wc -l
                """



# samtools index -@ {threads} {output} "{output}.bai"
