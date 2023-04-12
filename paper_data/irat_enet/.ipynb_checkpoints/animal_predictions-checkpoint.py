import os 
import sys
import numpy as np
import pandas as pd
import math
from Bio import SeqIO
from emergenet import Enet, save_model, load_model

# SUBTYPE: 'h1n1', 'h1n2', 'h3n2', 'h5n1', 'h5n2', 'h5n6', 'h5n8', 'h7', 'h9n2'
X, SUBTYPE = sys.argv
# X, SUBTYPE, N = sys.argv
# n = int(N)

# directory paths
DATA_DIR = 'raw_data/gisaid_animal/'
HUMAN_DATA_DIR = 'raw_data/gisaid_current/'
ENET_DIR = 'enet_models/current_enets/'
OUT_DIR = 'results/animal_predictions/'

# length to truncate sequences
NA_TRUNC = 449
HA_TRUNC = 550

# input: fasta file name, length to truncate each sequence, whether to represent strains as char arrays
# output: dataframe of sequences
def parse_fasta(file_name, trunc, seq_array = False):
    acc = []
    seq = []
    for record in SeqIO.parse(file_name, 'fasta'):
        if len(record.seq) < trunc:
            continue
        acc.append(record.id.split('|')[0])
        if seq_array:
            seq.append(np.array(record.seq[:trunc].upper()))
        else:
            seq.append(str(record.seq[:trunc].upper()))
    df = pd.DataFrame({'id':acc, 'sequence':seq})
    return df

# load animal data
ha_df = parse_fasta(DATA_DIR + SUBTYPE + '_ha.fasta', HA_TRUNC)
na_df = parse_fasta(DATA_DIR + SUBTYPE + '_na.fasta', NA_TRUNC)
df = ha_df.merge(na_df, how='inner', on='id').rename(columns={'sequence_x':'ha', 'sequence_y':'na'})
# split into 10 parts for H5N8
# chunk_size = int(df.shape[0]/10)
# start = chunk_size * n
# df = df.iloc[start:start + chunk_size]
# load human data
human_ha_df = parse_fasta(HUMAN_DATA_DIR + SUBTYPE + '_ha.fasta', HA_TRUNC, seq_array=True)
human_na_df = parse_fasta(HUMAN_DATA_DIR + SUBTYPE + '_na.fasta', NA_TRUNC, seq_array=True)

# load enet
ha_enet = load_model(ENET_DIR + SUBTYPE + '_ha.joblib')
na_enet = load_model(ENET_DIR + SUBTYPE + '_na.joblib')


ha_emergence_risk = []
na_emergence_risk = []
geometric_mean_risk = []

for i in range(len(df)):
    ha_risks = []
    na_risks = []
    geo_risks = []
    # repeat 10 times for variance computation
    for j in range(10):
        # access animal strain
        row = df.iloc[i]
        # initialize Enet
        enet_ha = Enet(seq=row['ha'], seq_trunc_length=HA_TRUNC, random_state=j)
        enet_na = Enet(seq=row['na'], seq_trunc_length=NA_TRUNC, random_state=j)
        # compute risk for ha and na
        emergence_risk_score_ha, variance_ha = enet_ha.emergence_risk(seq_df=human_ha_df, enet=ha_enet, sample_size=100)
        emergence_risk_score_na, variance_na = enet_na.emergence_risk(seq_df=human_na_df, enet=na_enet, sample_size=100)
        # add risk to list of risks for current animal strain
        ha_risks.append(emergence_risk_score_ha)
        na_risks.append(emergence_risk_score_na)
        geo_risks.append(math.sqrt(emergence_risk_score_ha * emergence_risk_score_na))
    
    # add risk lists to overall risk list
    ha_emergence_risk.append(ha_risks)
    na_emergence_risk.append(na_risks)
    geometric_mean_risk.append(geo_risks)

# add to dataframe
df['ha_risk'] = ha_emergence_risk
df['na_risk'] = na_emergence_risk
df['geometric_mean_risk'] = geometric_mean_risk

# save dataframe as csv
os.makedirs(OUT_DIR, exist_ok=True)
df.to_csv(OUT_DIR + SUBTYPE + '.csv', index=False)
# df.to_csv(OUT_DIR + SUBTYPE + '_' + N + '.csv', index=False)