import os 
import sys
import numpy as np
import pandas as pd
from domseq import DomSeq, load_model

# HEMISPHERE: 'north', 'south'
# SUBTYPE: 'h1n1', 'h3n2'
# N: 0-20
X, HEMISPHERE, SUBTYPE, N = sys.argv
n = int(N)

# Construct seasons
YEARS = []
if HEMISPHERE == 'north':
    for i in np.arange(2, 24):
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
            
# Directory paths
OUT_DIR = HEMISPHERE + '_' + SUBTYPE
NAME = OUT_DIR + '_' + YEARS[n] # For accessing data and Enet
NAME1 = OUT_DIR + '_' + YEARS[n + 1] # For PRED_DIR, SAVE_DATA_DIR, and pred_df
PRED_DIR = 'results/enet_predictions/seasonal_predictions/' + OUT_DIR + '/' + NAME1 + '.csv'
os.makedirs('results/enet_predictions/seasonal_predictions/' + OUT_DIR, exist_ok=True)
# Save data
SAVE_DATA_DIR = 'results/enet_predictions/extra_data/' + OUT_DIR + '/' + NAME1 + '/'
os.makedirs(SAVE_DATA_DIR, exist_ok=True)
# Data and Enet
DATA_DIR = 'data/merged/' + OUT_DIR + '/' + NAME + '.csv'
PRED_DATA_DIR = 'data/merged/' + OUT_DIR + '/pred/' + NAME + '.csv'
ENET_DIR = 'enet_models/' + OUT_DIR + '/' + NAME + '.joblib'

# Length to truncate sequences
TRUNC = 565
    
# Initialize the DomSeq
domseq = DomSeq(seq_trunc_length=TRUNC, random_state=42)
# Load data
seq_df = pd.read_csv(DATA_DIR)
pred_seq_df = pd.read_csv(PRED_DATA_DIR)
# Load enet
enet = load_model(ENET_DIR)
# Compute prediction sequence
pred_df = domseq.predict_domseq(seq_df, pred_seq_df, enet, 3, 3000, SAVE_DATA_DIR)
    
# Dataframe to store predictions
pred_df.to_csv(PRED_DIR, index=False)