import os, sys
import pandas as pd
from collections import Counter
from emergenet import Enet, save_model

X, IND = sys.argv
i = int(IND)

# Length to truncate sequences
NA_TRUNC = 449
HA_TRUNC = 560

# Load IRAT sequences
irat = pd.read_csv('data/animal/irat.csv')
row = irat.iloc[i]
virus_name = row['Influenza Virus'].replace('/',':')

# Directory paths
HUMAN_DIR = 'data/human/irat/'
ENET_DIR = 'enet_models/' + virus_name + '/'
os.makedirs(ENET_DIR, exist_ok=True)
RESULT_DIR = 'results/' + virus_name + '/'
os.makedirs(RESULT_DIR, exist_ok=True)

# Load human sequences
human = pd.read_csv(HUMAN_DIR + virus_name + '.csv', na_filter=False)

def evaluate(segment):
    TRUNC = NA_TRUNC
    if segment == 'HA':
        TRUNC = HA_TRUNC
    irat_seq = row[f'{segment} Sequence']
    if irat_seq != '-1':
        emergenet = Enet(seq=irat_seq, seq_trunc_length=TRUNC, random_state=42)
        human1 = human[human['segment'] == segment]
        human1 = human1[human1['sequence'].str.len() >= TRUNC]
        subtypes = Counter(human1[segment])
        for x in subtypes:
            if subtypes[x] < 15:
                continue
            # Use entire population for constructing Enet
            df = human1[human1[segment] == x]
            enet = emergenet.train(df, 10000, 8)
            save_model(enet, ENET_DIR + x + '.joblib')
            # Use only unique strains for inference
            df = df.drop_duplicates(subset=['sequence'])
            df = emergenet.risk(df, enet)
            df.to_csv(RESULT_DIR + x + '.csv', index=False)

evaluate('HA')
evaluate('NA')