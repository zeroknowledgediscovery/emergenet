import os 
import sys
import pandas as pd
from domseq import DomSeq, save_model, load_model

# N: random state
# SUBTYPE: 'h1n1', 'h3n2'
# YEAR: '15_16,' '16_17,' '17_18', '18_19', '19_20'
X, N, SUBTYPE, YEAR = sys.argv
# directory path, length to truncate sequences, and load data
NAME = 'north_' + SUBTYPE + '_' + YEAR
ENET_DIR = 'enet_models/' + SUBTYPE + '_' + YEAR
os.makedirs(ENET_DIR, exist_ok=True)
seq_df = pd.read_csv('data/north_' + SUBTYPE + '/' + NAME + '.csv')
    
TRUNC = 565

# initialize the DomSeq
domseq = DomSeq(seq_trunc_length=TRUNC, random_state=int(N))
# train enet
enet = domseq.train(seq_df=seq_df, sample_size=3000, n_jobs=1)
# save enet
save_model(enet=enet, outfile=ENET_DIR + '/enet_' + N + '.joblib')