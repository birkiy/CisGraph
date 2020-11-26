

def getControl(wildcards):
    controlName = list(set(controlDF.loc[tuple(wildcards.sortedPeak.split("."))]["controlName"]))
    yieldControl =  f"results/mapping/processed/control.{controlName[0]}.final.bam"
    yield yieldControl


rule peakCalling:
        input:
                case="results/mapping/processed/{sortedPeak}.final.bam",
                control=getControl
                # ["results/mapping/{tf}.{condition}.{rep}.final.bam",
                # "results/mapping/input.{inputName}.final.bam"]
        output:
                peak="results/peak/{sortedPeak}_peaks.narrowPeak",
                sortedPeak="results/peak/idr/{sortedPeak}.sorted.peaks.narrowPeak"
        log:
                "logs/peakCalling.{wildcards.sortedPeak}.log"
        message:
                "Executing peakCalling rule for {wildcards.sortedPeak}"
        shell:
                """
                macs2 callpeak \
                    -t {input.case} \
                    -c {input.control} \
                    -n {wildcards.sortedPeak} \
                    -B -f BAM -g hs -q 0.05 \
                    --outdir results/peak \
                    --shift -75 --extsize 150 2> {log}

                sort -k8,8nr {output.peak} > {output.sortedPeak}
                """
