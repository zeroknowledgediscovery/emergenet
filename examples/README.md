# Examples Using the Emergenet Package

## Quick Start

### Predicting Dominant Sequences

```
from emergenet.domseq import DomSeq, save_model, load_model

# initialize DomSeq
domseq = DomSeq(seq_trunc_length=566, random_state=42)

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

```
from emergenet.emergenet import Enet, save_model, load_model

# initialize Enet
enet = Enet(seq='target_sequence.fasta', seq_trunc_length=566, random_state=42)

# load data
df = enet.load_data(filepath='sequences.fasta')

# train enet
enet_model = enet.train(seq_df=df, sample_size=3000)
save_model(enet=enet_model, outfile='enet_modes/')

# compute emergence risk score
erisk, var = enet.emergence_risk(seq_df=df, enet=enet_model, sample_size=3000)
```
