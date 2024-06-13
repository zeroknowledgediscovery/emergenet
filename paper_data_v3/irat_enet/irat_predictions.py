import sys
import pandas as pd
from emergenet.emergenet import Enet


X, IND = sys.argv
i = int(IND)

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
SAVE_DIR = 'results/irat_predictions/' + virus_name + '/'
enet = Enet(analysis_date, ha_seq, na_seq, save_data=SAVE_DIR, random_state=42)
ha_risk, na_risk = enet.risk(enet_sample_size=10000)
