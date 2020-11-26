

# Activate the conda environment
source activate bpnet


# to generate common bw for replicate
bpnet align2bigwig analysis/mapping/AR.dht.rep1.final.bam analysis/mapping/AR.dht.rep2.final.bam --fragment-point=5prime --strand-spec AR.dht



bamCompare
