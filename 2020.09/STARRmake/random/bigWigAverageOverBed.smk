
#
# "bigWigAverageOverBed v2 - Compute average score of big wig over each bed, which may have introns.\n"
#   "usage:\n"
#   "   bigWigAverageOverBed in.bw in.bed out.tab\n"
#   "The output columns are:\n"
#   "   name - name field from bed, which should be unique\n"
#   "   size - size of bed (sum of exon sizes\n"
#   "   covered - # bases within exons covered by bigWig\n"
#   "   sum - sum of values over all bases covered\n"
#   "   mean0 - average over bases with non-covered bases counting as zeroes\n"
#   "   mean - average over just covered bases\n"
#   "Options:\n"
#   "   -stats=stats.ra - Output a collection of overall statistics to stat.ra file\n"
#   "   -bedOut=out.bed - Make output bed that is echo of input bed but with mean column appended\n"
#   "   -sampleAroundCenter=N - Take sample at region N bases wide centered around bed item, rather\n"
#   "                     than the usual sample in the bed item.\n"



bedGR = config["bedGR"]
rule bigWigOverBedRule:
        input:
                "results/bigwig/{raw}.extended.RPKM.bigWig",

        output:
                "results/coverage/{raw}.avg.tab"
        message:
                "Executing bigWigOverBedRule rule for {wildcards.raw}"
        shell:
                """
                bigWigAverageOverBed {input} {bedGR} {output}
                """
