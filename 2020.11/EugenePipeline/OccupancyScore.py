
import pyBigWig as BW
import numpy as np
import os
import pandas as pd



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
	for line in open("/home/ualtintas/genomeAnnotations/Regions/TSS.hg19.Idx.bed"):
		i += 1
		row = line.strip().split()
		chr, start, end, name, strand = row
		try:
			region = np.array(
				bw.values(
					chr, int(start), int(end)
				)
			)
			nBin = (int(end) - int(start))/50
			signal = np.mean(np.split(region, nBin), axis=1)
			if name in bigWig.keys():
				bigWig[name] = np.append(bigWig[name],signal)
			else:
				bigWig[name] = signal
			if i % 10000 == 0:
				print(i, TF)
		except:
			print(chr, start, end, name, strand)
			pass






	# for i in range(tss.shape[0]):
	# 	chr = tss.loc[i, "chr"]
	# 	start = tss.loc[i, "start"] - 25
	# 	end = tss.loc[i, "end"] + 25
	# 	name = tss.loc[i, "name"]
	# 	region = np.array(
	# 		bw.values(
	# 			chr, start, end
	# 		)
	# 	)
	# 	nBin = (end - start)/50
	# 	signal = np.mean(np.split(region, nBin), axis=1)
	# 	if name in bigWig.keys():
	# 		bigWig[name] = np.append(bigWig[name],signal)
	# 	else:
	# 		bigWig[name] = signal
	# 	if i % 10000 == 0:
	# 		print(i)
	#




index = bigWig.keys()
data = bigWig.values()
BigWig = pd.DataFrame(data, index=index).fillna(0)

BigWig.to_csv("bigwigs.occupancyScores.tsv", sep="\t")

Features = pd.DataFrame(BigWig , index=BigWig.index)


# Features = np.random.rand(1000, 1350)


nBin = int(nBin)
N,m = Features.shape
nBw = m // nBin
nFeature = nBw * nBin




FeatureIdx = Features.index
#Features = Features.iloc[:,:nFeature].values.reshape(N, nBw, nBin)

features = Features.values.reshape(N, nBw, nBin)


def scaling(mat, factor):
	"""
	Here we perform background scaling
	"""
	return mat / (mat + factor)

def SES(T, C): #
	"""
	T: ChIP experiment
	C: Bias


	SES normalization of the data.
	Returns the medium value for scaling and cutoff value
	"""
	p = np.cumsum(np.sort(T.flatten()))
	q = np.cumsum(np.sort(C.flatten()))
	diff = np.abs(p / p[-1] - q / q[-1])
	maxIndex = int(np.argmax(diff) * 0.8)
	cumSum = np.array(
        [
            float(p[maxIndex]),
            float(q[maxIndex])
            ]
        )
	sizeFactorsSES = cumSum.min()
	meanSES = [
        np.mean(
            np.sort(
                T.flatten()
                )[:maxIndex]
            ),
        np.mean(
            np.sort(
                C.flatten()
                )[:maxIndex]
            )
        ]
	medianScore = np.percentile(T.flatten()[maxIndex:], 50)
	maxScore = np.percentile(T.flatten()[maxIndex:], 95)
	# print("hello")
	return sizeFactorsSES, medianScore, maxScore

def finalSES(mat, C):
	"""
	Same as SES scaling, but only outputs the final scaled scores.
	"""
	scaleToMin, medianScore, maxScore = SES(mat, C)
	print(f"Background Score: {scaleToMin} \tMax Score: {maxScore} \tMedian Score: {medianScore}")
	return scaling(mat, medianScore)


scoresSES = np.zeros((nBw,3))
tmp = np.zeros((N, nBw, nBin))
occupancyScores = np.zeros((N, nBw))
C = np.linspace(0,1,num=len(features[:,i,:].flatten()))
C = C.reshape((N, nBin))

for i in range(nBw):
	scoresSES[i,:] = SES(features[:,i,:], C)
	tmp[:,i,:] = finalSES(features[:,i,:], C)
	occupancyScores[:,i] = np.mean(
		scaling(tmp[:,i,:],scoresSES[i,1])
		,axis=1
	)


occupancyScores = pd.DataFrame(occupancyScores, index=FeatureIdx).fillna(0)

occupancyScores.to_csv(f"{path}/occupancyScores.TSS.hg19.tsv", sep="\t", index=False, header=False)
