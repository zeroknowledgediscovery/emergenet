import os 
import sys
import pandas as pd
from domseq import DomSeq, save_model

# HEMISPHERE: 'north', 'south'
# SUBTYPE: 'h1n1', 'h3n2'
# YEAR: '02_03' to '22_23'
X, HEMISPHERE, SUBTYPE, YEAR = sys.argv
if HEMISPHERE == 'south':
    YEAR = YEAR[:2]
OUT_DIR = HEMISPHERE + '_' + SUBTYPE
NAME = OUT_DIR + '_' + YEAR 
os.makedirs('enet_models/' + OUT_DIR, exist_ok=True)

# Directory paths
DATA_DIR = 'data/merged/' + OUT_DIR + '/' + NAME + '.csv'
ENET_DIR = 'enet_models/' + OUT_DIR + '/' + NAME + '.joblib'

# Length to truncate sequences
TRUNC = 565

# Load data
seq_df = pd.read_csv(DATA_DIR)
# Initialize the DomSeq
domseq = DomSeq(seq_trunc_length=TRUNC, random_state=42)
# Train enet
enet = domseq.train(seq_df)
save_model(enet=enet, outfile=ENET_DIR)