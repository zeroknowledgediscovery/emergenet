import os 
import sys
import numpy as np
import pandas as pd
from domseq import DomSeq, save_model, load_model

# HEMISPHERE: 'north', 'south'
# SUBTYPE: 'h1n1', 'h3n2'
# N: 0-20
X, HEMISPHERE, SUBTYPE, N = sys.argv
n = int(N)

# construct seasons
YEARS = []
if HEMISPHERE == 'north':
    for i in np.arange(1, 23):
        YEAR = ''
        if i < 10:
            YEAR += '0' + str(i)
        else:
            YEAR += (str(i))
        if i + 1 < 10:
            YEAR += '_0' + str(i + 1)
        else:
            YEAR += '_' + str(i + 1)
        YEARS.append(YEAR)
else:
    for i in np.arange(2, 24):
        if i < 10:
            YEARS.append('0' + str(i))
        else:
            YEARS.append(str(i))
            
# directory paths
OUT_DIR = HEMISPHERE + '_' + SUBTYPE
NAME = OUT_DIR + '_ha_' + YEARS[n] # for accessing data, Enet, and storing distance matrix
NAME1 = OUT_DIR + '_ha_' + YEARS[n + 1] # for PRED_DIR and pred_df
PRED_DIR = 'results/enet_predictions/seasonal_predictions/' + OUT_DIR + '_ha/' + NAME1 + '_predictions.csv'
DIST_MATRIX_DIR = 'results/enet_predictions/distance_matrices/' + OUT_DIR + '_ha/' + NAME + '.csv'
DATA_DIR = 'raw_data/merged/' + OUT_DIR + '/' + NAME + '.csv' # notice no '_ha' after OUT_DIR
ENET_DIR = 'enet_models/' + OUT_DIR + '_ha/' + NAME + '.joblib'

# length to truncate sequences
TRUNC = 565
    
# initialize the DomSeq
domseq = DomSeq(seq_trunc_length=TRUNC, random_state=42)
# load data
seq_df = pd.read_csv(DATA_DIR)
# load enet
enet = load_model(ENET_DIR)
# compute prediction sequence
pred_accs, pred_names, pred_seqs, cluster_sizes = domseq.predict_domseq(seq_df, enet, 3, 1000, DIST_MATRIX_DIR)
    
# dataframe to store predictions
pred_df = pd.DataFrame({'season':3*[YEARS[n + 1]], 'acc':pred_accs,
                        'name':pred_names, 'sequence':pred_seqs, 'cluster_size':cluster_sizes})
pred_df = pred_df.sort_values(by='cluster_size',ascending=False)
pred_df.to_csv(PRED_DIR, index=False)