import sys, os
import pandas as pd
from datetime import datetime, timedelta
from emergenet.emergenet import Enet


X, IND, MONTHS = sys.argv
i = int(IND)
if MONTHS == '0':
    d = 1
else:
    d = int(MONTHS) * 30

# Load variant sequence
variant = pd.read_csv('data/variant.csv')
row = variant.iloc[i]
virus_name = row['name'].replace('/',':')
analysis_date = row['date']
# Analyze d days before emergence
input_date = datetime.strptime(analysis_date, "%Y-%m-%d") - timedelta(days=d)
analysis_date = input_date.strftime("%Y-%m-%d")
ha_seq = row['ha_sequence']
na_seq = row['na_sequence']
if ha_seq == '-1' or na_seq == '-1':
    sys.exit()

# Save directory
SAVE_DIR = 'results/variant_predictions/' + virus_name + '/' + MONTHS + '/'
if not os.path.exists(SAVE_DIR):
    enet = Enet(analysis_date, ha_seq, na_seq, save_data=SAVE_DIR, random_state=42)
    ha_risk, na_risk = enet.risk(sample_size=10000)
