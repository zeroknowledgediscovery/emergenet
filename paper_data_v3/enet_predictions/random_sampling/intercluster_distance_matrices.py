import os 
import sys
import numpy as np
import pandas as pd
from quasinet.qnet import qdistance
from domseq import load_model
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.manifold import MDS
from collections import Counter

# SUBTYPE: 'h1n1', 'h3n2'
# YEAR: '15_16,' '16_17,' '17_18', '18_19', '19_20'
X, SUBTYPE, YEAR = sys.argv

NAME = SUBTYPE + '_' + YEAR 
OUT_DIR = 'results/' + NAME + '/'

# load data
seq_df = pd.read_csv('data/data_dm/seq_df_' + NAME + '.csv')
dm = pd.read_csv('data/data_dm/dm_' + NAME + '.csv')

# clustering
embedding = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
dm_embed = embedding.fit_transform(dm)
bandwidth = estimate_bandwidth(dm_embed, quantile=0.3, random_state=42)
clustering = MeanShift(bandwidth=bandwidth)
clustering_predictions = clustering.fit_predict(dm_embed)

# sort unique clusters by size, and get top 2 largest
label_counts = Counter(clustering_predictions)
unique_clusters = sorted(label_counts.items(), key=lambda x: x[1], reverse=True)
classes = unique_clusters[:2]

# get cluster sequences
cluster_1 = dm.iloc[clustering_predictions == classes[0][0], clustering_predictions == classes[0][0]].values.flatten()
cluster_2 = dm.iloc[clustering_predictions == classes[1][0], clustering_predictions == classes[1][0]].values.flatten()
intercluster = dm.iloc[clustering_predictions == classes[0][0], clustering_predictions == classes[1][0]].values.flatten()

# save data
np.save(OUT_DIR + NAME + '_cluster1', cluster_1)
np.save(OUT_DIR + NAME + '_cluster2', cluster_2)
np.save(OUT_DIR + NAME + '_intercluster', intercluster)