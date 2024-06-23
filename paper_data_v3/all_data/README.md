# Emergenet Data

This directory contains all sequence metadata used in the paper *Emergenet: Fast Scalable Emergence Risk Assessment of Influenza A Strains Circulating In Non-human Hosts*. Data is obtained from [GISAID](https://platform.epicov.org/epi3/cfrontend#586f5f) and [NCBI](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus).

The Excel file `seq_metadata.xlsx` contains fours sheets (HA/NA, GISAID/NCBI) containing the metadata (name, subtype, accession) of each sequence.

## File Tree

``` txt
all_data
├── gisaid_metadata_ha.csv : (name, subtype, accession) for HA sequences from GISAID
├── gisaid_metadata_na.csv : (name, subtype, accession) for NA sequences from GISAID
├── ncbi_metadata_ha.csv : (name, subtype, accession) for HA sequences from NCBI
├── ncbi_metadata_na.csv : (name, subtype, accession) for NA sequences from NCBI
├── seq_metadata.xlsx : Contains sheets for each of the 4 CSVs above
├── total_sequences.csv : Summarizes total sequences of each subtype
├── total_sequences.ipynb : Compiles metadata and total_sequences.csv  
└── total_sequences.tex : Table of total_sequences.csv for paper
```
