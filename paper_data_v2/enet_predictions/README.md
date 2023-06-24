# Emergenet Seasonal Influenza Predictions

### Directory Structure

```
enet_predictions
├── dominant_sequences : compute cluster-wise dominant sequences for each season
│   ├── results : each subdirectory contains a dataframe of sequences; one dominant sequence per cluster
│   ├── compute_domseq.py : script to compute dominant sequences
│   ├── dominant_sequences.ipynb : example computation and visualization
│   └── run_compute_domseq.sh : shell script to run compute_domseq.py
├── enet_models : all HA Enet models trained
│   ├── north_h1n1_ha
│   ├── north_h3n2_ha
│   ├── south_h1n1_ha
│   └── south_h3n2_ha
├── raw_data
│   ├── gisaid : FASTA data from GISAID
│   ├── merged : CSV data combining GISAID and NCBI
│   ├── ncbi : FASTA data from NCBI
│   └── who : WHO seasonal vaccine recommendations
├── results
│   ├── enet_predictions : Emergenet recommended seqeunces for each season
│   └── enet_who_comparison : compares Emergenet and WHO recommendations to dominant sequences
├── tables : LaTeX tables for paper
├── data_preprocessing.ipynb : cleaning and merging raw FASTA data from GISAID and NCBI
├── domseq.py : local version of DomSeq module of the Emergnet package, which finds and predicts dominant sequences
├── emergenet_predictions.ipynb : aggregates predications from enet_predictions.py
├── enet_predictions.py : script to predict dominant strains using Emergenet
├── enet_train.py : script to train Enet models
├── enet_who_comparison.ipynb : compares Emergenet and WHO recommendations to dominant sequences
├── make_tables.ipynb : make tables for paper (found in tables directory)
├── run_enet_predictions.sh : shell script to run enet_predictions.py
└── run_enet_train.sh : shell script to run enet_train.py
```
