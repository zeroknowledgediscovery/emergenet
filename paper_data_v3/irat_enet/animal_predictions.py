import sys, math
import pandas as pd
from emergenet.emergenet import Enet


X, IND = sys.argv
i = int(IND)
num_partions = 100
SAVE_DIR = 'results/animal_predictions/'

# Load subset of animal sequences
animal = pd.read_csv('data/animal.csv')
partition_size = math.ceil(len(animal) / num_partions)
start = i * partition_size
end = None if i == num_partions - 1 else start + partition_size
animal = animal.iloc[start:end]
if len(animal) == 0:
    sys.exit()

# Compute risks for each sequence
ha_risks = []
na_risks = []
for j in range(len(animal)):
    row = animal.iloc[j]
    virus_name = row['name'].replace('/',':')
    ha_seq = row['ha_sequence']
    na_seq = row['na_sequence']
    enet = Enet('PRESENT', ha_seq, na_seq, random_state=42)
    ha_risk, na_risk = enet.risk()
    ha_risks.append(ha_risk)
    na_risks.append(na_risk)

# Save results
animal['ha_risk'] = ha_risks
animal['na_risk'] = na_risks
animal.to_csv(SAVE_DIR + 'animal' + IND + '.csv', index=False)
