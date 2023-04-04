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
│   ├── irat_average_qdistances.csv : compares risk assesment from IRAT and Qnet q-distance
│   ├── irat_average_qdistances.tex : SI Tab. 16
│   ├── irat_average_qdistances_all_sequences.csv : comparison for H1- and H3- using all data of same subtype
│   ├── irat_average_qdistances_all_sequences.tex : LaTeX file for irat_average_qdistances_all_sequences.csv
│   ├── irat_data.csv : replicated CDC IRAT table
│   └── irat_data_all_sequences.csv : CDC IRAT table with only H1- and H3- strains
├── irat_data_collection.ipynb : create irat_data.csv
├── irat_data_collection_all_sequences.ipynb : create irat_data_all_sequences.csv
├── irat_figures.ipynb : create figures with irat_data.csv
├── irat_figures_all_sequences.ipynb : create figures with irat_data_all_sequences.csv
├── irat_qnet_comparison.ipynb : compares risk assesment from IRAT and Qnet q-distance
└── irat_qnet_comparison_all_sequences.ipynb : comparison for H1- and H3- using all data of same subtype
```
