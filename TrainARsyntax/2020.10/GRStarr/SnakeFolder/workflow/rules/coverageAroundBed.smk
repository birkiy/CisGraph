



bedGR = config["bedGR"]




def getPooledBw():
        samples = set(sampleDF["sampleName"])
        return expand("results/bigwig/{sample}.merged.{{type}}.bigWig", sample=samples)




rule coverageAroundBed:
        input:
                bed=expand("peaks/{bed}", bed=["common.GR.peaks.bed", "negativeControlGR.final.bed"]),
                BW=getPooledBw()
        output:
                npz="results/coverage/coverage.{type}.npz",
                pdf="results/coverage/plot.{type}.pdf"
        threads:
                16
        message:
                "Executing coverageAroundSummits rule"
        shell:
                """
                computeMatrix reference-point -S \
                {input.BW} \
                -R {input.bed} \
                --referencePoint=center\
                -a 1000 -b 1000 \
                --sortRegions descend -p {threads} \
                -o {output.npz}

                plotHeatmap -m {output.npz} \
                --whatToShow 'heatmap and colorbar' \
                --colorMap "Blues" --missingDataColor 1 \
                --sortRegions descend \
                --heatmapHeight 20 --heatmapWidth 7 \
                -out {output.pdf}
                """
