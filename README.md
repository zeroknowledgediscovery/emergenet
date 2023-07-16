# Emergenet

### Directory Structure

```
Emergenet
├── emergenet : files for Emergenet package
├── examples : examples using the Emergenet package
├── paper_data_v2 : results for current version of the paper
├── paper_data : results for older versions of the paper, not used in current version
├── paper_data_old : results for older versions of the paper, not used in current version
└── tex : contains LaTeX and PDF files for paper
```

## Description
- Computing and predicting dominant viral strains
- Superfast risk assessment of emerging pathogens


## Installation
[PyPI](https://pypi.org/project/emergenet/)

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
* scikit-learn
* statsmodels

## Usage

### Predicting Dominant Sequences

```
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

```
from emergenet.emergenet import Enet, save_model, load_model

# initialize Enet
enet = Enet(seq='target_sequence.fasta', seq_trunc_length=550, random_state=42)

# load data
df = enet.load_data(filepath='sequences.fasta')

# train enet
enet_model = enet.train(seq_df=df, sample_size=1000)
save_model(enet=enet_model, outfile='enet_modes/')

# compute emergence risk score
erisk, var = enet.emergence_risk(seq_df=df, enet=enet_model, sample_size=1000)
```
 
### Examples

Examples are located [here](https://github.com/zeroknowledgediscovery/emergenet/tree/main/examples).

## Documentation

For more documentation, see [here](https://zeroknowledgediscovery.github.io/emergenet/).
