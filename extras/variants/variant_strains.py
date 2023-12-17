import os 
import sys
import numpy as np
import pandas as pd
from tqdm.notebook import tqdm
from Bio import SeqIO
from emergenet import Enet, save_model, load_model, irat_risk

X, N = sys.argv
n = int(N)

# directory paths
DATA_DIR = 'data/variants/'
HUMAN_DIR = 'data/human/'
ENET_DIR = 'enet_models/'
OUT_DIR = 'results/'

# length to truncate sequences
NA_TRUNC = 445
HA_TRUNC = 565

names = []
ha_metadata = []
# na_metadata = []
ha_seq = []
# na_seq = []
risk_ha = []
var_ha = []
# risk_na = []
# var_na = []
# geom_mean = []
# irat_emergence = []
# irat_impact = []

# get name
STRAIN = 'A/Wisconsin/71/2016'
names.append(STRAIN)
NAME = STRAIN.replace('/', ':')

# initialize Enets
enet_ha = Enet(seq=DATA_DIR+NAME+'_ha.fasta', seq_trunc_length=HA_TRUNC, random_state=42)
# enet_na = Enet(seq=DATA_DIR+NAME+'_na.fasta', seq_trunc_length=NA_TRUNC, random_state=42)
ha_metadata.append(enet_ha.seq_metadata)
# na_metadata.append(enet_na.seq_metadata)
ha_seq.append(enet_ha.seq)
# na_seq.append(enet_na.seq)

# load human data
df_ha = enet_ha.load_data(filepath=HUMAN_DIR+'h1n2_ha_wisconsin.fasta')
# df_na = enet_na.load_data(filepath=HUMAN_DIR+'h1n2_na.fasta')

# # train Enets
enet_ha1 = enet_ha.train(seq_df=df_ha, sample_size=3000)
save_model(enet=enet_ha1, outfile=ENET_DIR+NAME+'_ha1.joblib')
# # enet_na1 = enet_na.train(seq_df=df_na, sample_size=3000)
# # save_model(enet=enet_na1, outfile=ENET_DIR+NAME+'_na.joblib')

# # compute risk scores
risk_score_ha, variance_ha, min_score, min_seq = enet_ha.emergence_risk(seq_df=df_ha, enet=enet_ha1, sample_size=1000)
# risk_score_na, variance_na = enet_na.emergence_risk(seq_df=df_na, enet=enet_na1, sample_size=1000)
# geom_mean_risk_score = np.sqrt(risk_score_ha * risk_score_na)
# irat_emergence_prediction, irat_impact_prediction = irat_risk(risk_score_ha, risk_score_na)

# # save data
risk_ha.append(risk_score_ha)
var_ha.append(variance_ha)
# risk_na.append(risk_score_na)
# var_na.append(variance_na)
# geom_mean.append(geom_mean_risk_score)
# irat_emergence.append(irat_emergence_prediction)
# irat_impact.append(irat_impact_prediction)

results = pd.DataFrame({'name':names,
                        'ha_metadata':ha_metadata,
#                         'na_metadata':na_metadata,
                        'ha_seq':ha_seq,
#                         'na_seq':na_seq,
                        'risk_ha':risk_ha,
                        'var_ha':var_ha,
                        'min_score':[min_score],
                        'min_seq':[min_seq]
#                         'risk_na':risk_na,
#                         'var_na':var_na,
#                         'geom_mean':geom_mean,
#                         'irat_emergence':irat_emergence,
#                         'irat_impact':irat_impact
                       })
results.to_csv(OUT_DIR+NAME+'1.csv')