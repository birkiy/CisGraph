import pybedtools

# set up 3 different bedtools
a = pybedtools.BedTool('a.bed')
b = pybedtools.BedTool('b.bed')
c = pybedtools.BedTool('c.bed')

ARhi0 = pybedtools.BedTool("ARhi.0M_peaks.narrowPeak")
ARhi1 = pybedtools.BedTool("ARhi.1nM_peaks.narrowPeak")
ARmo0 = pybedtools.BedTool("ARmo.0M_peaks.narrowPeak")
ARmo1 = pybedtools.BedTool("ARmo.1nM_peaks.narrowPeak")
PC80 = pybedtools.BedTool("PC8.0M_peaks.narrowPeak")
PC81 = pybedtools.BedTool("PC8.1nM_peaks.narrowPeak")


pybedtools.contrib.venn_maker(ARhi0,ARhi1,ARmo0,ARmo1,PC80,PC81)


(a-b-c).count()  # unique to a
(a+b-c).count()  # in a and b, not c
(a+b+c).count()  # common to all

intervene pairwise -i *.narrowPeak --names=ARhi0,ARhi1,ARmo0,ARmo1,PC80,PC81
