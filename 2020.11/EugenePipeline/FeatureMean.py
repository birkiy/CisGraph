

import pyBigWig as BW
import numpy as np
import os
import pandas as pd



from matplotlib import pyplot, patches
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import gridspec, colors
# from pygraphviz import *
from statannot import add_stat_annotation


BEDS = {
	"GRE": "/groups/lackgrp/ll_members/berkay/GenomicRegions/GR/totalGRE.cut.bed",
	"ARE": "/groups/lackgrp/ll_members/berkay/GenomicRegions/AR/ARBS.bed",
	"TSS": "/groups/lackgrp/ll_members/berkay/GenomicRegions/TSS/TSS.hg19.Idx.cut.bed"
}

path = "results/bigwig/genomcov/merged/"
bigWig = {}
i = 0

for file in os.listdir(path):
	# if (file.find("genomcov") == -1) or (file.find("merged") == -1):
	# 	print(f"{file}: passed")
	# 	continue
	print(f"{file}: started")
	TF = file[:-7]
	fullPath = f"results/bigwig/genomcov/merged/{file}"
	bw  =  BW.open(fullPath)
	for BED in BEDS:
		# if BED in ["TSS", "GRE"]:
		# 	continue
		for line in open(BEDS[BED]):
			i += 1
			row = line.strip().split()
			chr, start, end, name = row
			try:
				center = int(start) + int(np.floor((int(end) - int(start))/2))
				start = center - 350
				end = center + 350
				# start += center - 350
				# end -= center - 350
				region = np.array(
					bw.values(
						chr, start, end
					)
				)
				nBin = (end - start) // 50
				signal = np.mean(np.nan_to_num(np.split(region, nBin)), axis=1)
				if name in bigWig.keys():
					bigWig[name] = np.append(bigWig[name],signal)
				else:
					bigWig[name] = signal
				if i % 10000 == 0:
					print(i, TF)
			except:
				print(chr, start, end, name)
				pass









index = bigWig.keys()
data = bigWig.values()
BigWig = pd.DataFrame(data, index=index).fillna(0)

# updated 26.11.2020
BigWig.to_csv("bigwigs.occupancyScores.tsv", sep="\t")

Features = pd.DataFrame(BigWig , index=BigWig.index)


nBin = int(nBin)
N,m = Features.shape
nBw = m // nBin
nFeature = nBw * nBin




FeatureIdx = Features.index
#Features = Features.iloc[:,:nFeature].values.reshape(N, nBw, nBin)

features = Features.values.reshape(N, nBw, nBin)




path = "results/bigwig/genomcov/merged/"
Big = []
Big2 = []
for filename in os.listdir(path):
	start = filename[:-7]
	co_factor = filename[:-7]
	Big += [str(co_factor)]
	Big2 += [filename[:-7]]



def scaling (array, factor):
	"""
	Here we perform background scaling
	"""
	mean = factor[0]
	array_2 = array+mean
	#return 1/(1+ np.exp(-power))
	return array/array_2

def SES(array):
	"""
	SES normalization of the data.
	Returns the medium value for scaling and cutoff value

	"""
	values = np.zeros(4)
	play = array.flatten()
	#play2 = array2.flatten()
	#plays = np.append(play,play2,axis=0)
	s_c = sum(play)
	plays2= play/s_c
	plays2 = np.sort(plays2)
	plays3 = np.cumsum(plays2)
	x= np.linspace(0,1,num=len(plays2))
	a=np.argmax(x-plays3)
	scale_min = plays2[a]*s_c
	mean= np.percentile(plays2[a:]*s_c,50)
	return mean, scale_min



def scaling3(array):
	"""
	Same as SES scaling, but only outputs the final scaled scores.

	"""
	values = np.zeros(4)
	play = array.flatten()
	#play2 = array2.flatten()
	#plays = np.append(play,play2,axis=0)
	s_c = sum(play)
	plays2= play/s_c
	plays2 = np.sort(plays2)
	plays3 = np.cumsum(plays2)
	#print(plays3[-1])
	x= np.linspace(0,1,num=len(plays2))
	a=np.argmax(x-plays3)
	scale_min = plays2[a]*s_c
	#print(plays[:a])
	mean= np.percentile(plays2[a:]*s_c,50)
	values[0] = mean
	scale = np.percentile(plays2[a:]*s_c,95)
	print('Background Score',scale_min,'Max Score',scale,'Median Score',mean)
	#scale = np.percentile(plays,95)
	#array[array > scale] = scale
	#array[array<scale_min] = 0
	#array = array/scale
	array = scaling(array,[mean,mean-scale_min])
	#array2[array2 > scale] = scale
	#array2[array2<scale_min] = 0
	#array2 = array2/scale
	#array2 = scaling(array2,[mean,mean-scale_min])
	return array






play = features.reshape(N,nBw,nBin)

BigWigs2 = np.zeros((N, nBw, nBin))
Promoter_scale = np.zeros((N, nBw))

SES_c =np.zeros((90,2))

Promoters = features

for i in range(nBw):
	#Background_DHT = np.mean(Background_check.iloc[:,DHT])
	#Background_ETOH = np.mean(Background_check.iloc[:,ETOH])
	SES_c[i,:] = SES(play[:,i,:])
	BigWigs2[:,i,:] = scaling3(play[:,i,:])
	Promoter_scale[:,i] = np.mean(scaling(Promoters[:,i,:],SES_c[i,:]),axis=1)

Promoters_index = FeatureIdx

Data = pd.DataFrame(Promoter_scale,index=Promoters_index,columns=Big)


BigWigs2 = np.mean(BigWigs2, axis=2)

# Data= pd.DataFrame(BigWigs2, index=BigWigs.index,columns=Big )
# Data.to_csv('ARBS_chip_seq.csv')
# Data.to_csv('Promoters_chip_seq.csv')

chip = Data



BEDS = {
	"GRE": "/groups/lackgrp/ll_members/berkay/GenomicRegions/GR/totalGRE.cut.bed",
	"nGR": "/groups/lackgrp/ll_members/berkay/GenomicRegions/GR/negativeControlGR.final.bed",
	"ARE": "/groups/lackgrp/ll_members/berkay/GenomicRegions/AR/ARBS.bed",
	"TSS": "/groups/lackgrp/ll_members/berkay/GenomicRegions/TSS/TSS.hg19.Idx.cut.bed",
	"con": "~/ARBSs/regions/cons-arbs.bed",
	"ind": "~/ARBSs/regions/ind-arbs.bed",
	"non": "~/ARBSs/regions/Non-Active-ARBS.bed",
	"nAR": "~/ARBSs/regions/negativeControl.ARBS.bed"
}

chip["nodeClass"] = "NA"

for BED in BEDS:
	bed = pd.read_csv(BEDS[BED], sep="\t", names=["chr", "start", "end", "name"])
	chip.loc[chip.index.isin(bed["name"]), "nodeClass"] = BED



fig = plt.figure(figsize=[20,20])
gs = gridspec.GridSpec(ncols=6, nrows=5)
plt.subplots_adjust(wspace=0.4)

i = 0
for bw in Big:
	params = dict(x='nodeClass',
              y=bw)
	boxPairs = [("nAR", "con"),
	            ("con", "ind"),
	            ("con", "non"),
	            ("ind", "non"),
	            ("con", "TSS")
	            ]
	ax1 = fig.add_subplot(gs[i])
	ax1 = sns.boxplot(**params, data=Data)
	plt.ylabel(f"{bw}", fontsize=14)
	plt.xlabel("Node Classes", fontsize=14)
	for tick in ax1.xaxis.get_major_ticks():
	    tick.label.set_fontsize(13)
	#add_stat_annotation(ax1, **params, data=Data, box_pairs=boxPairs, test='Mann-Whitney')
	i += 1


fig.savefig("chip.pdf")
