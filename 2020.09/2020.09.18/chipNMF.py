

from time import time

from sklearn.feature_extraction.text import TfidfVectorizer, CountVectorizer
from sklearn.decomposition import NMF, LatentDirichletAllocation
from sklearn.datasets import fetch_20newsgroups

from Functions.Packages import *

chipARBS = pd.read_csv(f"{dataRoot}/ARBS_chip_seq_2.csv")


n_samples = 2000
n_features = 1000
n_components = 10
n_top_words = 20




data, _ = fetch_20newsgroups(shuffle=True, random_state=1,
                             remove=('headers', 'footers', 'quotes'),
                             return_X_y=True)

data_samples = data[:n_samples]

# convert the text to a tf-idf weighted term-document matrix

vectorizer = TfidfVectorizer(max_features=2000, min_df=10, stop_words='english')

X = vectorizer.fit_transform(chipARBS)

idx_to_word = np.array(vectorizer.get_feature_names())

# apply NMF


nmf = NMF(n_components=20, solver="mu")

W = nmf.fit_transform(X)

H = nmf.components_
 













tfidf_vectorizer = TfidfVectorizer(max_df=0.95, min_df=2,
                                   max_features=n_features,
                                   stop_words='english')


t0 = time()
tfidf = tfidf_vectorizer.fit_transform(data_samples)
print("done in %0.3fs." % (time() - t0))


# Use tf (raw term count) features for LDA.
print("Extracting tf features for LDA...")
tf_vectorizer = CountVectorizer(max_df=0.95, min_df=2,
                                max_features=n_features,
                                stop_words='english')
t0 = time()
tf = tf_vectorizer.fit_transform(data_samples)
print("done in %0.3fs." % (time() - t0))
print()



# Fit the NMF model
print("Fitting the NMF model (Frobenius norm) with tf-idf features, "
      "n_samples=%d and n_features=%d..."
      % (n_samples, n_features))
t0 = time()
nmf = NMF(n_components=n_components, random_state=1,
          alpha=.1, l1_ratio=.5).fit(tfidf)
print("done in %0.3fs." % (time() - t0))
