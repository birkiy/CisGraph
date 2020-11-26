

bws = [
        "A549.H3K27ac.dex.0h.rep1",
        "A549.H3K27ac.dex.0h.rep2",
        "A549.H3K27ac.dex.0h.rep3",
        "A549.H3K27ac.dex.8h.rep1",
        "A549.H3K27ac.dex.8h.rep2",
        "A549.H3K4me1.dex.0h.rep1",
        "A549.H3K4me1.dex.0h.rep2",
        "A549.H3K4me1.dex.0h.rep3",
        "A549.H3K4me1.dex.8h.rep1",
        "A549.H3K4me1.dex.8h.rep2",
        "A549.H3K4me1.dex.8h.rep3",
        "A549.H3K4me3.dex.0h.rep1",
        "A549.H3K4me3.dex.0h.rep2",
        "A549.H3K4me3.dex.0h.rep3",
        "A549.H3K4me3.dex.8h.rep1",
        "A549.H3K4me3.dex.8h.rep2",
        "A549.POLR2A.etoh.1h.rep1",
        "A549.POLR2A.etoh.1h.rep2",
        "A549.POLR2A.dex.1h.rep1",
        "A549.POLR2A.dex.1h.rep2",
        "A549.GR.dex.0h.rep1",
        "A549.GR.dex.0h.rep2",
        "A549.GR.dex.8h.rep1",
        "A549.GR.dex.8h.rep2",
        "A549.GR.dex.8h.rep3"
        ]



desiredOutputList = [
        f"coverage/{raw}.{type}.avg.tab"
        for raw in bws
        for type in ["con", "ind", "non"]
]

def getInput():
        for out in desiredOutputList:
                yield out
