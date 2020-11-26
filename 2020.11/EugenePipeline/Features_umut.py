
import pyBigWig as BW
import numpy as np
import os
import pandas as pd
from Functions import scaled_counts
from sklearn.metrics import mean_squared_error,mutual_info_score
from sklearn import linear_model
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn import preprocessing
from sklearn.decomposition import PCA
from scipy.stats import spearmanr,binned_statistic
from matplotlib import cm

"""
Defining the path and recover the ordering of the factor
"""
path = 'bigwigs/'
Big = []
Big2 = []
for filename in os.listdir('bigwigs'):
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
	
"""
Loads in values from ARBS regions and the new promoter regions
-90_Bigwigs_mean_v2; ARBS regions
-Promoters_mean; Promoters
Fills in zeros if needed
"""
BigWigs  = pd.read_csv('90_Bigwigs_mean_v2.csv', index_col=0)
Promoters = pd.read_csv('Promoters_mean.csv',index_col=0)

BigWigs = BigWigs.fillna(0)
Promoters = Promoters.fillna(0)



Features = pd.DataFrame(BigWigs , index = BigWigs.index)

#print(new_counts['Log_Values'].iloc[:10])
#print(BigWigs.iloc[:10,:10])
BigWigs = Features[(Features.index.str.endswith('ARBS')) | (Features.index.str.startswith('NegativeControl'))]

feature = BigWigs.iloc[:,:1350]
print(feature.shape)

N,m = feature.shape
N_2,m_2 = Promoters.shape

"""
Reshape the dataframes
"""


Promoters_index= Promoters.index
Promoters = Promoters.iloc[:,:1350].values.reshape(N_2,90,15)


play  =  feature.values.reshape(N,90,15)
Promoter_scale = np.zeros((N_2,90))
#BigWigs2 = np.max(play,axis=2)
BigWigs2 = np.zeros((N,90,15))

SES_c =np.zeros((90,2))

"""
Normalize using SES method

First determine the mean using the ARBS regions data
Finally apply scaling to promoter data
"""

for i in range(90):

	#Background_DHT = np.mean(Background_check.iloc[:,DHT])
	#Background_ETOH = np.mean(Background_check.iloc[:,ETOH])
	SES_c[i,:] = SES(play[:,i,:])

	BigWigs2[:,i,:] = scaling3(play[:,i,:])
	Promoter_scale[:,i] = np.mean(scaling(Promoters[:,i,:],SES_c[i,:]),axis=1)
	
	
BigWigs2 = np.mean(BigWigs2, axis=2)

Data= pd.DataFrame(BigWigs2, index=BigWigs.index,columns=Big )
Data.to_csv('ARBS_chip_seq.csv')
Data = pd.DataFrame(Promoter_scale,index=Promoters_index,columns=Big)
Data.to_csv('Promoters_chip_seq.csv')
