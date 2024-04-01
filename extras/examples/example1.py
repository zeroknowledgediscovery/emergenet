from emergenet.emergenet import Enet, save_model, load_model, irat_risk

DATA_DIR = 'example_data/emergenet/'

enet_ha = Enet(seq='../extras/variants/data/variants/A:Alberta:01:2020_(H1N2)v_ha.fasta', seq_trunc_length=550, random_state=42)

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

print('Emergenet Risk Score:', risk_score_ha)
print('Variance:', variance_ha)
