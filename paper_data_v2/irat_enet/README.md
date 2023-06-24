# CDC IRAT Risk vs. Emergenet

### Directory Structure

```
irat_qnet
├── enet_models : all Enet models trained
│   ├── current_enets : trained on current data for each subtype
│   └── irat_enets : trained on data available at time of IRAT analysis
├── figures
├── raw_data
│   ├── gisaid : GISAID amino acid fasta data to compute qnet, named each IRAT strain
│   ├── gisaid_animal : GISAID current animal sequence data
│   ├── gisaid_current : GISAID current human sequence data
│   └── irat_sequences : sequence of IRAT strain
├── results
│   ├── animal_predictions: uses GLM models to predict risk of current circulating animal Influenza strains
│   ├── glm_models: GLM models fit on Emergenet E-distance
│   ├── predict_linear_bounds: linear bounds for Fig. 3 in paper
│   ├── irat_average_qdistances.csv : Emergenet E-distance scores added to IRAT table, using data available at time of IRAT analysis
│   ├── irat_average_qdistances_current.csv : Emergenet E-distance scores added to IRAT table, using data available currently
│   ├── irat_data.csv : replicated CDC IRAT table
│   ├── irat_predictions.csv : predicted IRAT scores from irat_average_qdistances.csv using GLM models
│   └── irat_predictions_current.csv : predicted IRAT scores from irat_average_qdistances_current.csv using GLM models
├── animal_predictions.ipynb : evaluate results of running animal_predictions.py
├── animal_predictions.py : script to predict risk scores of current circulating animal Influenza strains
├── emergenet.py : local version of Emergenet module of the Emergnet package, which predicts risk of a strain
├── irat_data_collection.ipynb : create irat_data.csv
├── irat_enet_comparison.ipynb : compares risk assesment from IRAT and Emergenet E-distance
├── irat_enet_current.ipynb : same as irat_enet_comparison.ipynb but using current data
├── irat_figures.ipynb : create figures for paper
├── linear_models.ipynb : create GLM models and makes irat_predictions.csv and irat_predictions_current.csv
├── predict_linear_bounds.ipynb : linear bounds for Fig. 3 in paper
└── run_animal_predictions.sh : shell script to run animal_predictions.py
```
