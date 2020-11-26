

tfs=["AR", "FOXA1"]
tfs=["ARID", "BRG1", "HOB13", "TLE3", "TRIM28", "WDHD1"]
conditions=["dht", "dmso"]
conditions=["r1881", "veh"]
reps=["rep1", "rep2"]




{k:v for k in K for K in zip(tfs, conditions) for v in reps }


def getOut():
    for condition in conditions:
        for tf in tfs:
            for rep in reps:
                yield [
                    f"analysis/peaks/{tf}.{condition}.{rep}_peaks.xls",
                    f"analysis/peaks/{tf}.{condition}.{rep}_peaks.narrowPeak"
                ]

def getOut():
    for outFile, reps in sampleDict.items():






getOut()
rule all:
        input:
                getOut()

rule PeakCalling:
        input:
                chip="analysis/mapping/{tf}.{condition}.{rep}.final.bam",
                control=expand("analysis/mapping/input.{rep}.final.bam",  proxy=[] if {{tf}} in ["AR", "FOXA1"] else [None])
        output:
                xls="analysis/peaks/{tf}.{condition}.{rep}_peaks.xls",
                peak="analysis/peaks/{tf}.{condition}.{rep}_peaks.narrowPeak",
                sortedPeak="analysis/peaks/idr/{tf}.{condition}.{rep}.sorted.peaks.narrowPeak"
        log: "logs/{tf}.{condition}.{rep}.log"
        message: "Executing peakCalling rule for {wildcards.tf}.{wildcards.condition}.{wildcards.rep}"
        shell:
                """
                macs2 callpeak \
                    -t {input.chip} \
                    -c {input.control} \
                    -n {wildcards.tf}.{wildcards.condition}.{wildcards.rep} \
                    -B -f BAM -g hs -q 0.05 \
                    --outdir analysis/peaks \
                    --shift -75 --extsize 150 2> {log}

                sort -k8,8nr {output.peak} > {output.sortedPeak}
                """
