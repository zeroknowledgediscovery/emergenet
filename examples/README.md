# Examples Using the Emergenet Package

## Quick Start

### Predicting Dominant Sequences


```python
from emergenet.domseq import DomSeq, save_model, load_model

# initialize DomSeq
domseq = DomSeq(seq_trunc_length=565, random_state=42)

# load data from current time period
df = domseq.load_data(filepath='sequences.fasta')

# compute dataframe of cluster-wise dominant sequences for current time period
dominant_sequences = domseq.compute_domseq(seq_df=df)

# train enet
enet_model = domseq.train(seq_df=df, sample_size=3000)
save_model(enet=enet_model, outfile='enet_modes/')

# load candidate sequences for recommendation
pred_df = domseq.load_data(filepath='pred_sequences.fasta')

# compute dataframe of cluster-wise recommendation sequences for next time period
prediction_sequences = domseq.predict_domseq(seq_df=df, pred_seq_df=pred_df, enet=enet_model, sample_size=3000)
```

### Evaluating Sequence Risk

```python
from emergenet.emergenet import Enet, predict_irat_emergence

# Initialize the Enet
enet = Enet(analysis_date, ha_seq, na_seq, save_data=SAVE_DIR, random_state=42)

# Estimate the Enet risk scores
ha_risk, na_risk = enet.risk(sample_size=10000)

# Map the Enet risk scores to the IRAT risk scale
irat, irat_low, irat_high = predict_irat_emergence(ha_risk, na_risk)
```
