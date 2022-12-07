# Emergenet

### Directory Structure

```
Emergenet
├── emergenet: files for Emergenet package
├── examples: examples using the Emergenet package
├── irat_qnet : IRAT Emergence Risk vs. Qnet q-distance comparison
├── qnet_predictions : compute dominant strains and Qnet predictions
└── tex : contains LaTeX and PDF files for paper
```

## Description
- Computing predicting dominant strains
- Superfast risk assessment of emerging pathogens


## Installation

To install with pip:

```
pip install emergenet
# to update
pip install emergenet --upgrade
```

### Dependencies

* [quasinet](https://github.com/zeroknowledgediscovery/quasinet/)
* numpy 
* pandas 
* Levenshtein 
* biopython

## Usage

### Predicting Dominant Sequences

```
from emergenet.domseq import DomSeq, save_model, load_model

# initialize DomSeq
domseq = DomSeq(seq_trunc_length=566, random_state=42)

# load data
df = domseq.load_data(filepath='sequences.fasta')

# compute dominant sequence for current time period
dom_id, dom_seq = domseq.compute_domseq(seq_df=df, sample_size=1000)

# train enet
enet = domseq.train(seq_df=df, sample_size=1000)

# compute prediction sequence for next time period
pred_id, pred_seq = domseq.predict_domseq(seq_df=df, enet=enet, sample_size=1000)
```

### Evaluating Sequence Risk

```
from emergenet.emergenet import Enet, save_model, load_model

# initialize enet
enet = Enet(seq='target_sequence.fasta', seq_trunc_length=550, random_state=42)

# load data
df = enet.load_data(filepath='sequences.fasta')

# train enet
enet = enet.train(seq_df=df, sample_size=1000)

# compute emergence risk score
erisk, var = enet.emergence_risk(seq_df=df, enet=enet, sample_size=1000)
```
 
### Examples

Examples are located [here](https://github.com/zeroknowledgediscovery/emergenet/tree/main/examples).

## Documentation

For more documentation, see [here](https://zeroknowledgediscovery.github.io/emergenet/).
