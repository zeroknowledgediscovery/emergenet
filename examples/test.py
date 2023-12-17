from emergenet.emergenet import Enet, save_model, load_model, irat_risk
import numpy as np

target_seq='../extras/variants/data/variants/A:Alberta:01:2020_(H1N2)v_'
HA_ext='ha.fasta'
NA_ext='na.fasta'
HA_len=550
DATA_DIR = 'example_data/emergenet/'

enet_ha = Enet(seq=target_seq+HA_ext, seq_trunc_length=HA_len, random_state=22)

print('Initializing Enet with fasta file\n---------------------------------\n')
print(enet_ha.seq_metadata)
print(enet_ha.seq)
print('Length of target sequence:', len(enet_ha.seq))

# load fasta data
df_ha = enet_ha.load_data(filepath=DATA_DIR+'ha_sequences.fasta')

print('Number of sequences:', len(df_ha))
# load enet
enet_ha1 = load_model(filepath=DATA_DIR+'ha_enet.joblib')

# compute risk score
risk_score_ha, variance_ha = enet_ha.emergence_risk(seq_df=df_ha.head(1), enet=enet_ha1, sample_size=1000)

#avg_ha, min_ha, max_ha, var_ha = enet_ha.emergence_risk_qsampling(seq_df=df_ha, enet=enet_ha1, sample_size=1000, qsamples=10, steps=10)



print('Emergenet Risk Score:', risk_score_ha)
print('Variance:', variance_ha)

enet_na = Enet(seq='../extras/variants/data/variants/A:Alberta:01:2020_(H1N2)v_na.fasta', seq_trunc_length=449, random_state=42)

print('Initializing Enet with fasta file\n---------------------------------\n')
print(enet_na.seq_metadata)
print(enet_na.seq)
print('Length of target sequence:', len(enet_na.seq))

# load fasta data
df_na = enet_na.load_data(filepath=DATA_DIR+'na_sequences.fasta')

print('Number of sequences:', len(df_na))

enet_na1 = load_model(filepath=DATA_DIR+'na_enet.joblib.gz')

risk_score_na, variance_na = enet_na.emergence_risk(seq_df=df_na, enet=enet_na1, sample_size=1000)

print('Emergenet Risk Score:', risk_score_na)
print('Variance:', variance_na)



#avg_na, min_na, max_na, var_na = enet_na.emergence_risk_qsampling(seq_df=df_na, enet=enet_na1, sample_size=1000, qsamples=10, steps=10)

#print('Emergenet Risk Score:', avg_na)
#print('Bounds:', [min_na, max_na])
#print('Variance:', var_na)

# compute geometric mean emergence risk
geom_mean_risk_score = np.sqrt(risk_score_ha * risk_score_na)
irat_emergence_prediction, irat_impact_prediction = irat_risk(risk_score_ha, risk_score_na)

print('Geometric Mean of HA and NA risk scores:', round(geom_mean_risk_score, 6))
print('Emergenet prediction of IRAT emergence estimate:', round(irat_emergence_prediction, 1))
print('Emergenet prediction of IRAT impact estimate:', round(irat_impact_prediction, 1))
