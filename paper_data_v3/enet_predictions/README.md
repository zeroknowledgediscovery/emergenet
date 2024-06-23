# Predicting Seasonal Vaccine Strains

Prediction of the seasonal vaccine strain as a close match to a historical strain allows out-of-sample validation against past World Health Organization (WHO) recommendations for the flu shot, which is reformulated February and September of each year for the northern and southern hemispheres, respectively, based on a cocktail of historical strains determined via global surveillance. For each year of the past two decades (2002-03 to 2022-23), we constructed a separate Emergenet for the southern and northern flu seasons using HA strains available from the previous season, for H1N1 and H3N2 subtypes. We use these models to predict seasonal vaccine strains, which we compare to the WHO's official vaccine recommendations.

## Replicating Paper Results

### 1) Aggregate Data

Run `data_collection.ipynb`. 

**NOTE: it is important that you do this before running any other code!**

### 2) Compute Seasonal Vaccine Strain Predictions

Follow the instructions in `emergenet_predictions.ipynb`.

### 3) Comparison with WHO Vaccine Recommendations

Run `enet_who_comparison.ipynb`.

### 4) Make Tables for Paper

Run `make_tables.ipynb`

### 5) Random Sampling Effects

Test if random sampling for training Emergenets affects E-distance between two strains.

Follow the instructions in `random_sampling/enet_random_sampling.ipynb`.

## File Tree

``` txt
enet_predictions
├── data : FASTA files containing sequences
│   ├── gisaid : https://platform.epicov.org/epi3/
│   ├── ncbi : https://www.ncbi.nlm.nih.gov/labs/virus/vssi/
│   ├── who : https://www.who.int/teams/global-influenza-programme/vaccines/who-recommendations
│   └── huddleston.csv : https://cdn.elifesciences.org/articles/60067/elife-60067-fig11-data1-v2.csv
├── random_sampling : Test the effect of random sampling in Emergenet training
│   ├── data
│   ├── results
│   ├── domseq.py
│   ├── enet_random_sampling.ipynb
│   ├── enet_random_sampling_test.py
│   ├── enet_random_sampling_train.py
│   ├── intercluster_distance_matrices.py
│   ├── run_enet_random_sampling_test.sh
│   ├── run_enet_random_sampling_train.sh
│   └── run_intercluster_distance_matrices.sh
├── results
│   ├── enet_predictions : Seasonal vaccine strain predictions
│   └── enet_who_comparison : Compare Emergenet predictions with WHO recommendations
├── tables : Tables and figures for paper
├── data_collection.ipynb : Creates CSV data from FASTA
├── domseq.py : Local version of emergenet.domseq module
├── emergenet_predictions.ipynb : Aggregate seasonal vaccine strain predictions
├── enet_predictions.py : Predict seasonal vaccine strains
├── enet_train.py : Train Enets for predicting seasonal vaccine strains
├── enet_who_comparison.ipynb : Compare Emergenet predictions with WHO recommendations
├── make_tables.ipynb : Make tables for paper
├── plot_predictions.ipynb : Example usage of plotseq.py
├── plotseq.py : Module for visualizing the strain space
├── run_enet_predictions.sh : Parallelize enet_predictions.py
└── run_enet_train.sh : Parallelize enet_train.py
```