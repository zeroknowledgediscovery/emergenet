# Qnet Predictions

```
EmergeNet
├── irat_qnet
├── qnet_predictions
│   ├── misc
│   │   ├── missing_prediction_2019_2020.ipynb : H3N2 NA Northern Hemisphere 2019-2020 prediction was missing, this notebook computes it
│   │   └── north_h3n2_na_18.fasta : use 2018-2019 sequence data to compute prediction
│   ├── qnet_models : zipped Qnets by hemisphere, sub-type (H1N1/H3N2), segment (HA/NA), year (ex. 20 = 2020-2021)
│   ├── raw_data
│   │   ├── gisaid : GISAID amino acid fasta data by hemisphere, sub-type (H1N1/H3N2), segment (HA/NA), year (ex. 20 = 2020-2021)
│   │   └── ncbi : all NCBI amino acid fasta data through June 2022, by sub-type (H1N1/H3N2), segment (HA/NA)
│   ├── results
│   │   ├── dominant_sequences_<YEAR>.csv
│   │   ├── influenza_qnet_predictions_<YEAR>.csv
│   │   └── influenza_qnet_predictions_3cluster_<YEAR>.csv
│   ├── tables : CSV tables for SI-Tab. 3-15
│   ├── dominant_sequences_<YEAR>.ipynb : finds dominant sequences for specific flu season using edit distance
│   ├── influenza_qnet_predictions_<YEAR>.ipynb : predicts dominant strain for next flu season using qdistance; single cluster and 3 cluster
│   ├── qnet_improvement.ipynb : computes values for Tab. 1 
│   └── update_tables.ipynb : updates CSV in "tables" directory using data from "results" directory
└── tex
```
