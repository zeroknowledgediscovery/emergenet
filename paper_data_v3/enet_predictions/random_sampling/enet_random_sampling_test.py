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

ENET_DIR = 'enet_models/'
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
wanted_names_1 = dm.columns[clustering_predictions == classes[0][0]]
seq_df_1 = seq_df.iloc[wanted_names_1]
wanted_names_2 = dm.columns[clustering_predictions == classes[1][0]]
seq_df_2 = seq_df.iloc[wanted_names_2]

# random sample of sequence pairs from each cluster
sampled_seqs_1 = []
sampled_seqs_2 = []
for i in range(10):
    df_1 = seq_df_1.sample(2, random_state=i).values
    df_2 = seq_df_2.sample(2, random_state=i).values
    sampled_seqs_1.append(df_1[0][1] + ', ' + df_1[1][1])
    sampled_seqs_2.append(df_2[0][1] + ', ' + df_2[1][1])

# compute qdistances under different Enets
qdists_1 = []
qdists_2 = []
for n in range(1, 101):
    enet = load_model(ENET_DIR + NAME + '/enet_' + str(n) + '.joblib')
    qdist_1 = []
    qdist_2 = []
    for i in range(10):
        seqs_1 = seq_df_1.sample(2, random_state=i)['sequence'].values
        seqs_2 = seq_df_2.sample(2, random_state=i)['sequence'].values
        qdist_1.append(qdistance(np.array(list(seqs_1[0])), np.array(list(seqs_1[1])), enet, enet))
        qdist_2.append(qdistance(np.array(list(seqs_2[0])), np.array(list(seqs_2[1])), enet, enet))
    qdists_1.append(qdist_1)
    qdists_2.append(qdist_2)

# save data
data_df_1 = pd.DataFrame(np.array(qdists_1).T)
data_df_1['seqs'] = sampled_seqs_1
data_df_1 = data_df_1.set_index('seqs')
data_df_1.to_csv(OUT_DIR + NAME + '_cluster1_data.csv')
data_df_1.T.describe().T.to_csv(OUT_DIR + NAME + '_cluster1_summary.csv')

data_df_2 = pd.DataFrame(np.array(qdists_2).T)
data_df_2['seqs'] = sampled_seqs_2
data_df_2 = data_df_2.set_index('seqs')
data_df_2.to_csv(OUT_DIR + NAME + '_cluster2_data.csv')
data_df_2.T.describe().T.to_csv(OUT_DIR + NAME + '_cluster2_summary.csv')