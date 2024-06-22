import sys
import pandas as pd
from emergenet.emergenet import Enet


X, IND, RANDOM_STATE = sys.argv
i = int(IND)
random_state = int(RANDOM_STATE)

# Load IRAT sequence
irat = pd.read_csv('data/irat.csv')
row = irat.iloc[i]
virus_name = row['Influenza Virus'].replace('/',':')
analysis_date = row['Date of Risk Assessment']
ha_seq = row['HA Sequence']
na_seq = row['NA Sequence']
if ha_seq == '-1' or na_seq == '-1':
    sys.exit()

# Save directory
BASE_SAVE_DIR = 'results/irat_predictions/' + virus_name + '/sem/'
SAVE_DIR = BASE_SAVE_DIR + str(random_state) + '/'
enet = Enet(analysis_date, ha_seq, na_seq, save_data=SAVE_DIR, random_state=random_state)
ha_risk, na_risk = enet.risk(enet_sample_size=0.75)
