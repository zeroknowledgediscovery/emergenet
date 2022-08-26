# CDC IRAT Emergence Risk vs. Qnet Q-distance Comparison

## IRAT vs. Qnet Comparison
1. Update IRAT strain data in `raw_data/irat_sequences` from [IRAT](https://www.cdc.gov/flu/pandemic-resources/monitoring/irat-virus-summaries.htm#H1N2variant)
    - File name is "VIRUS_NAME_na.fasta" or "VIRUS_NAME_ha.fasta"
    - Update `irat_data.csv` using `irat_data_collection.ipynb`
2. For each strain previously analyzed by IRAT:
    - Collect strains from GISAID one year leading up to month of analysis
    - For example, the "A/swine/Shandong/1207/2016" strain was assessed by IRAT in July 2020, so we will use human H1N1 strains circulating between July 1, 2019 through June 30, 2020
    - Note: had difficulty finding 'A/duck/New York/1996', only NA available, do not use in final results
    - For the following strains, only use upper bound of date due to small sample size
        - H1N2, H5N1, H5N6, H7N7, H9N2
    - The following strains have no human strains available
        - H5N2, H5N8, H7N8, H10N8
    - Strains with 'Qnet Sample' = -1 have no available human strains
3. See `irat_qnet_comparison.ipynb`
    - Construct a Qnet using these strains **if there are more than 30 strains in the population for both NA and HA**
    - Compute the average q-distance among the strain in question and the circulating human strains for both NA and HA
        - Average the NA and HA averages for the "Both Average Qdistance" column, which is used in the final graph

## IRAT vs. Qnet Comparison - All Sequences
**H1N1, H1N2, and H3N2 only**
1. Update `irat_data_all_sequences.csv` using `irat_data_collection_all_sequences.ipynb`
2. For each strain previously analyzed by IRAT
    - Collect all strains of that variety from NCBI
    - Note: had difficulty finding "A/duck/New York/1996", only NA available
3. See `irat_qnet_comparison_all_sequences.ipynb`
    - Construct a Qnet using these strains
        - For H1N2, the HA is similar to that of H1N1 and the NA to that of H3N2, so use H1N1 HA for HA Qnet and H3N2 NA for NA Qnet
    - Compute the average q-distance among the strain in question and the circulating human strains for both NA and HA
        - Average the NA and HA averages for the "Both Average Qdistance" column, which is used in the final graph

## File Tree
```
Emergenet
├── irat_qnet
│   ├── figures
│   │   ├── irat_combined.png : Fig. 3
│   │   └── irat_split.png : SI Fig. 3
│   ├── qnet_models : zipped Qnets named after IRAT strains
│   │   └── zipped Qnets using all NCBI strains from raw_data/ncbi
│   ├── raw_data
│   │   ├── gisaid : GISAID amino acid fasta data to compute qnet, named each IRAT strain, NA and HA
│   │   ├── irat_sequences : sequence of IRAT strain, NA and HA
│   │   └── ncbi : all NCBI amino acid fasta data through May, 1 2022, by sub-type (H1N1/H3N2), segment (HA/NA)
│   ├── results
│   │   ├── irat_average_qdistances.csv : compares risk assesment from IRAT and Qnet q-distance
│   │   ├── irat_average_qdistances.tex : SI Tab. 16
│   │   ├── irat_average_qdistances_all_sequences.csv : comparison for H1- and H3- using all data of same subtype
│   │   ├── irat_average_qdistances_all_sequences.tex : LaTeX file for irat_average_qdistances_all_sequences.csv
│   │   ├── irat_data.csv : replicated CDC IRAT table
│   │   └── irat_data_all_sequences.csv : CDC IRAT table with only H1- and H3- strains
│   ├── irat_data_collection.ipynb : create irat_data.csv
│   ├── irat_data_collection_all_sequences.ipynb : create irat_data_all_sequences.csv
│   ├── irat_figures.ipynb : create figures with irat_data.csv
│   ├── irat_figures_all_sequences.ipynb : create figures with irat_data_all_sequences.csv
│   ├── irat_qnet_comparison.ipynb : compares risk assesment from IRAT and Qnet q-distance
│   ├── irat_qnet_comparison_all_sequences.ipynb : comparison for H1- and H3- using all data of same subtype
├── qnet_predictions
└── tex
```
