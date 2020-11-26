import pyBigWig as BW
import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt

"""
Generates chipseq profiles given the bed file

Defining path to location of bigwigs
"""
path = 'bigwigs/'
Big = {}

"""
Promoters means new promoter file

"""
for filename in os.listdir('bigwigs'):

	"""
	co_factor is the name of factor from bigwig file. 
	From there, we are going to define the path to each file
	"""
	co_factor = filename[:-7]


	Pathline = path + str(co_factor) + '.bigWig'
	print(Pathline)
	bw  =  BW.open(Pathline)
	
	"""
	Opening the bed file, to read in region data. 
	File is too large to open on its own, so I read in line by line.
	Depending on the work, load in a different bed file.
	"""
	for line in open('TSS.hg19.Idx.bed'):
		
		
		"""
		Read in chromosome location, start and end locations for each file in the 
		"""
		cols = line.strip().split()
		#print(cols)
		vals = np.array(bw.values(cols[0], int(cols[1]), int(cols[2])+50))
		"""
		Split each region into 50 bp chunks, and fill in Nan values and find average value
		"""
		vals2  = np.split(vals, 15)
		vals2 = np.nan_to_num(vals2)
		vals2 = np.mean(vals2,axis=1)

		#
		#print(np.mean(vals2,axis=0))
		"""
		Save the bigwig data 
		"""
		key =  cols[3]
		if key in Big.keys():
			Big[key] = np.append(Big[key],vals2)
		else:
			Big[key] = vals2


#print(Big['NegativeControl-AR_motif_2572'])
index = Big.keys()
data = Big.values()
BigWig = pd.DataFrame(data, index= index)
print(np.shape(BigWig))
BigWig.to_csv('Promoters_mean.csv')


BigWigs  = pd.read_csv('Promoters_mean.csv', index_col=0)

