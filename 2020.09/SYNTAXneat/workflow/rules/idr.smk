
#
# sampleId, repId = glob_wildcards("results/peak/idr/{{sample}}.{repId}.sorted.peaks.narrowPeak")
#
# print(repId)

repAvail = False

def getIdr(wildcards):
    global repAvail
    reps = list(set(samples.loc[tuple(wildcards.sample.split(".")), "rep"]))
    if len(reps) < 2:
        repAvail = False
    else:
        repAvail = True
    return expand("results/peak/idr/{sample}.{rep}.sorted.peaks.narrowPeak", rep=reps, **wildcards)



rule idr:
        input:
                getIdr                # getIdrInput
                # expand("results/peak/idr/{sample}.{{repIdr}}.sorted.peaks.narrowPeak",
                # sample=list(set(samples.loc[(wildcards.tf, wildcards.condition), "rep"])))
        output:
                "results/peak/idr/{sample}.idr"
        log:
                "logs/idr.{sample}.log"
        message:
                "Executing idr rule for {wildcards.sample}"
        run:
               if repAvail:
                        shell(
                        """
                        idr --samples {input} \
                        --input-file-type narrowPeak \
                        --rank p.value \
                        --output-file {output} \
                        --plot
                        """
                        )
               else:
                        shell(
                        """
                        echo "Replicates are not available for {wildcards.sample}" > {output}
                        """)






# shell:
#         if False:
#             """
#             idr --samples {input} \
#             --input-file-type narrowPeak \
#             --rank p.value \
#             --output-file {output} \
#             --plot
#             """
#         else:
#             """
#             echo "Replicates are not available for {wildcards.sample} > {output}"
#

mkdir hey
N=$(( $RANDOM % 10))
for j in $(seq 1 $N); do echo -n $j > hey/$j.txt; done
