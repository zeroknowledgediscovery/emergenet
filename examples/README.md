# Examples Using the Emergenet Package

Run `emergenet_examples.ipynb` for an example estimating emergence risk of an animal strain with `emergenet.emergenet` 

Run `domseq_examples.ipynb` for an example predicting a future dominant strain with `emergenet.domseq`. 

## Quick Start

### Estimating Emergence Risk with `emergenet.emergenet`

To evaluate a strain, use `estimate_risk.py`. You only need to provide the HA and NA sequences. Run `python estimate_risk.py -h` for more arguments.

Here is a detailed example using an IRAT strain, A/Indiana/08/2011, evaluated at present time.

```bash
$ HA=MKTIIAFSCILCLIFAQKLPGSDNSMATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGGICNSPHQILDGKNCTLIDALLGDPHCDDFQNKEWDLFVERSTAYSNCYPYYVPDYATLRSLVASSGNLEFTQESFNWTGVAQGGSSYACRRGSVNSFFSRLNWLYNLNYKYPEQNVTMPNNDKFDKLYIWGVHHPGTDKDQTNLYVQASGRVIVSTKRSQQTVIPNIGSRPWVRGVSSIISIYWTIVKPGDILLINSTGNLIAPRGYFKIQSGKSSIMRSDAHIDECNSECITPNGSIPNDKPFQNVNKITYGACPRYVKQNTLKLATGMRNVPEKQTRGIFGAIAGFIENGWEGMVDGWYGFRHQNSEGTGQAADLKSTQAAINQITGKLNRVIKKTNEKFHQIEKEFSEVEGRIQDLEKYVEDTKIDLWSYNAEILVALENQHTIDLTDSEMSKLFERTRRQLRENAEDMGNGCFKIYHKCDNACIGSIRNGTYDHDIYRNEALNNRFQIKGVQLKSGYKDWILWISFAISCFLLCVVLLGFIMWACQKGNIRCNICI
$ NA=MNPNQKIITIGSVSLIIATICFLMQIAILVTTVTLHFKQHDYNSPPNNQAMLCEPTIIERNTTEIVYLTNITIEKEICPKLAEYRNWSKPQCNITGFAPFSKDNSIRLSAGGDIWVTREPYVSCDPDKCYQFALGQGTTLNNGHSNNTVHDRTPYRTLLMNELGVPFHLGTRQVCMAWSSSSCHDGKAWLHVCITGNDNNATASFIYNGRLVDSIGSWSKNILRTQESECVCINGTCTVVMTDGSASGKADTKILFVEEGKIVHISTLSGSAQHVEECSCYPRFPGVRCVCRDNWKGSNRPIVDINVKNYSIVSSYVCSGLVGDTPRKSDSVSSSYCLDPNNEKGGHGVKGWAFDDGNDVWMGRTINETLRLGYETFKVIEGWSKANSKLQTNRQVIVEKGDRSGYSGIFSVEGKSCINRCFYVELIRGRKEETKVWWTSNSIVVFCGTSGTYGTGSWPDGADINLMPI
$ python estimate_risk.py $HA $NA --risk_sample_size=100

Estimated IRAT Emergence Score: 6.50
Time taken: 31.28 seconds
```

Here is a detailed example using the same IRAT strain, A/Indiana/08/2011, evaluated at the time of IRAT assessment.

```python
import pandas as pd
from emergenet.emergenet import Enet, predict_irat_emergence

DATA_DIR = 'data/emergenet/'

# Load IRAT sequence - A/Indiana/08/2011
irat_df = pd.read_csv(DATA_DIR+'irat.csv')
row = irat_df.iloc[20]

# We need the analysis date, and HA and NA sequences
# Optionally, we can proved a save_data directory
analysis_date = row['Date of Risk Assessment']
ha_seq = row['HA Sequence']
na_seq = row['NA Sequence']
SAVE_DIR = 'data/emergenet/example_results/'

# Initialize the Enet
enet = Enet(analysis_date=analysis_date, 
            ha_seq=ha_seq, 
            na_seq=na_seq, 
            save_data=SAVE_DIR, 
            random_state=42)

# Estimate the Enet risk scores
ha_risk, na_risk = enet.risk(risk_sample_size=100)

# Map the Enet risk scores to the IRAT risk scale
irat, irat_low, irat_high = predict_irat_emergence(ha_risk=ha_risk, 
                                                   na_risk=na_risk)
```

### Predicting Seasonal Vaccine Strains with `emergenet.domseq` 

```python
import pandas as pd
from emergenet.domseq import DomSeq
from emergenet.utils import save_model, load_model

DATA_DIR = 'data/domseq/'

# Initialize the DomSeq
domseq = DomSeq(seq_trunc_length=565, random_state=42)

# Load data from current time period (2021-2022 season)
df = pd.read_csv(DATA_DIR+'north_h1n1_21_22.csv')

# Train enet
enet = domseq.train(seq_df=df, sample_size=3000, n_jobs=1)

# Load candidate sequences for recommendation
# This includes all human H1N1 strains up until the 2021-2022 season
candidate_df = pd.read_csv(DATA_DIR+'north_h1n1_21_22_pred.csv')

# Compute prediction sequences (return predictions from top 3 largest clusters)
pred_df = domseq.predict_domseq(seq_df=df, 
                                pred_seq_df=candidate_df, 
                                enet=enet_model, 
                                n_clusters=3, 
                                sample_size=3000)

# Compute a single prediction for the dominant strain
single_pred_seq = domseq.predict_single_domseq(pred_seqs=pred_df, 
                                               pred_seq_df=candidate_df)
```
