import os 
import sys
import numpy as np
import pandas as pd
from domseq import DomSeq, save_model, load_model

# HEMISPHERE: 'north', 'south'
# SUBTYPE: 'h1n1', 'h3n2'
# SEGMENT: 'ha', 'na'
# N: 0-19
X, HEMISPHERE, SUBTYPE, SEGMENT, N = sys.argv
n = int(N)

# construct seasons
YEARS = []
if HEMISPHERE == 'north':
    for i in np.arange(3, 23):
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
    for i in np.arange(3, 23):
        if i < 10:
            YEARS.append('0' + str(i))
        else:
            YEARS.append(str(i))
            
# directory paths
OUT_DIR = HEMISPHERE + '_' + SUBTYPE
NAME = OUT_DIR + '_' + YEARS[n]
NAME_SEG = OUT_DIR + '_' + SEGMENT + '_' + YEARS[n]
# save data
SAVE_DATA_DIR = 'results/' + OUT_DIR + '_' + SEGMENT + '/' + NAME_SEG + '/'
os.makedirs(SAVE_DATA_DIR, exist_ok=True)
# data
DATA_DIR = '../raw_data/merged/' + OUT_DIR + '/' + NAME + '.csv'

# length to truncate sequences
TRUNC = 565
if SEGMENT == 'na':
    TRUNC = 468
    
# initialize the DomSeq
domseq = DomSeq(seq_trunc_length=TRUNC, random_state=42)
# load data
seq_df = pd.read_csv(DATA_DIR)
if SEGMENT == 'na':
    seq_df = seq_df.drop(columns=['acc','sequence'])
    seq_df.rename(columns={'acc_na':'acc','sequence_na':'sequence'}, inplace=True)
# compute dominant sequences
dom_seqs = domseq.compute_domseq(seq_df, save_data=SAVE_DATA_DIR)
dom_seqs.to_csv(SAVE_DATA_DIR + 'dom_seqs.csv', index=False)
