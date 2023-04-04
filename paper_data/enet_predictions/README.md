# Emergenet Seasonal Influenza Predictions

### Directory Structure

```
enet_predictions
├── enet_models: all HA Enet models trained
│   ├── north_h1n1_ha
│   ├── north_h3n2_ha
│   ├── south_h1n1_ha
│   └── south_h3n2_ha
├── raw_data
│   ├── gisaid: FASTA data from GISAID
│   ├── merged: CSV data combining GISAID and NCBI
│   ├── ncbi: FASTA data from NCBI
│   └── who: WHO seasonal vaccine recommendations
├── results
│   ├── dominant_sequences: dominant sequences for each season (edit distance centroid)
│   ├── enet_predictions: Emergenet recommended seqeunces for each season
│   ├── enet_who_comparison: compares Emergenet and WHO recommendations to dominant sequences
│   ├── num_seqs_north.csv: number of sequences per northern season
│   └── num_seqs_south.csv: number of sequences per southern season
├── tables: LaTeX tables for paper
├── dominant_sequences.ipynb: finds dominant strains using edit distance
├── domseq.py: local version of DomSeq module of the Emergnet package, which finds and predicts dominant sequences
├── emergenet_predictions.ipynb: predicts dominant strains using Emergenet
├── enet_predictions.py: script to predict dominant strains using Emergenet
├── enet_train.py: script to train Enet models
├── enet_who_comparison.ipynb: compares Emergenet and WHO recommendations to dominant sequences
├── make_tables.ipynb: make tables for paper (found in tables directory)
├── run_enet_predictions.sh: shell script to run enet_predictions.py
└── run_enet_train.sh: shell script to run enet_train.py
```
