# Estimating Emergence Risk of Non-human Strains

Our primary claim is the ability to estimate the emergence potential of animal strains that do not yet circulate in humans. We validate this claim by training Emergenet models for strains analyzed by the CDC's [Influenza Risk Assessment Tool (IRAT)](https://www.cdc.gov/flu/pandemic-resources/national-strategy/risk-assessment.htm) and variant strains (animal strains that emerged in humans).

## Replicating Paper Results

### 1) Aggregate Data

Run `data_collection.ipynb`. 

**NOTE: it is important that you do this before running any other code!**

### 2) Risk Assessment of IRAT Strains

We train Emergenet models for strains analyzed by the CDC's [Influenza Risk Assessment Tool (IRAT)](https://www.cdc.gov/flu/pandemic-resources/national-strategy/risk-assessment.htm). We construct Emergenet models for HA and NA sequences using subtype-specific human strains, collected within the year prior to the assessment date. We use these models to estimate the risk of strains previously assessed by IRAT.

Follow the instructions in `emergenet_predictions.ipynb`.

#### Final IRAT Model

The final IRAT results and model reported in the paper are generated with a slightly older version of the `emergenet` package. To run this, please go into the `final_irat_model` subdirectory.

Follow the instructions in `final_irat_model/emergenet_predictions.ipynb`.

#### Random Subsampling Effects

We explore sampling effects. We randomly sample 75% of the strain population for Enet training and risk evaluation, computing the mean and standard error over 20 random seeds.

Follow the instructions at the bottom `emergenet_predictions.ipynb`.

### 3) Risk Assessment of Variant Strains

A second demonstration of our approach is obtained by computing the E-risk scores for "variants" (animal Influenza A strains isolated in human hosts), which might be expected to pose high risk due to their successful replication in human cells, even if the possibility of HH transmission is not yet observed or guaranteed. We follow the same procedure as with the IRAT strains to estimate the risk of each variant at six, three, one, and zero months before their collection date.

Follow the instructions in `variant_predictions.ipynb`.

For the final figure paper, see `variant_new.ipynb`.

### 4) Risk Assessment of Animal Strains

Finally, we estimate the IRAT scores of all 6,354 wild Influenza A animal viruses collected globally between January
2020 and January 2024. 

Follow the instructions in `animal_predictions.ipynb`.

## File Tree

``` txt
irat_enet
├── data
│   ├── animal : All animal strains from 1/1/2020 - 1/1/2024
│   ├── human : All human strains from 1/1/2010 - 1/1/2024
│   ├── irat : Strains previously analyzed by IRAT
│   └── variant : Variant strains (animal strains that have emerged in humans)
├── emergenet : Local version of emergenet.emergenet module
├── final_irat_model : Generates IRAT results reported in the paper
│   ├── data : Data for IRAT risk estimation
│   ├── results : Results for each IRAT strain
│   ├── data_collection.ipynb : Creates CSV data from FASTA
│   ├── emergenet.py : Older Local version of emergenet.emergenet module
│   ├── irat_predictions.ipynb : Analyze risk estimates of IRAT strains
│   ├── irat_predictions.py : Estimate risk of IRAT strains
│   └── run_irat_predictions.sh : Parallelize irat_predictions.py
├── results
│   ├── animal_predictions : Risk estimates of animal strains
│   ├── irat_predictions : Risk estimates of IRAT strains
│   └── variant_predictions : Risk estimates of variant strains
├── tables : Tables and figures for paper
├── animal_predictions.ipynb : Analyze risk estimates of animal strains
├── animal_predictions.py : Estimate risk of animal strains
├── data_collection.ipynb : Creates CSV data from FASTA
├── irat_predictions.ipynb : Analyze risk estimates of IRAT strains
├── irat_predictions.py : Estimate risk of IRAT strains
├── irat_predictions_sem.py : Estimate risk of IRAT strains with subsampling
├── run_animal_predictions.sh : Parallelize animal_predictions.py
├── run_irat_predictions_sem.sh : Parallelize irat_predictions_sem.py
├── run_irat_predictions.sh : Parallelize irat_predictions.py
├── run_variant_predictions.sh : Parallelize variant_predictions.py
├── variant_new.ipynb : Final variant figure reported in the paper
├── variant_predictions.ipynb : Analyze risk estimates of variant strains
└── variant_predictions.py : Estimate risk of variant strains
```