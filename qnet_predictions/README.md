# Seasonal Dominant Strain and Qnet Predictions

## Dominant Strain
To compute dominant strain for the **YEAR1 - YEAR2** flu season:
1. Download GISAID amino acid data with the following filters:
    - Host: Human
    - Flu Season: 
        - Northern strains from 10/01/**YEAR2** - 5/01/**YEAR2**
        - Southern strains from 04/01/**YEAR1** - 10/1/**YEAR1**
        - Flu season dates from [CDC](https://www.cdc.gov/flu/school-business/travelersfacts.htm)
    - Segment: HA (4) and NA (6)
    - Add data to `raw_data/gisaid` with the appropriate naming convention
        - HEMISPHERE_SEQUENCE_SEGMENT_SEASON
        - ex. north_h1n1_ha_20.fasta
2. Duplicate any `dominant_sequences_<YEAR>.ipynb` and modify it for current season (more instructions in notebook)
    - Levenshtein Centroid: $\widehat{x}^{dom} = argmin_{x\in P^t} \sum_{y \in P^t} \theta(x,y)$
        - Where $P^t$ is the sequence population at time $t$.
        - $\theta(x,y)$ is the edit distance between x and y

## Qnet Predictions
To compute Qnet predictions strain for the **YEAR1 - YEAR2** flu season:
1. Download GISAID amino acid data with the following filters:
    - Host: Human
    - Flu Season: 
        - Northern strains from 10/01/**YEAR2** - 5/01/**YEAR2**
        - Southern strains from 04/01/**YEAR1** - 10/1/**YEAR1**
        - Flu season dates from [CDC](https://www.cdc.gov/flu/school-business/travelersfacts.htm)
    - Segment: HA (4) and NA (6)
    - Add data to `raw_data/gisaid` with the appropriate naming convention
        - HEMISPHERE_SEQUENCE_SEGMENT_SEASON
        - ex. north_h1n1_ha_20.fasta
2. Download NCBI amino acid data, change upper bound of year to May 1, **YEAR1**
    - H1N1 HA: [Link](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein&VirusLineage_ss=H1N1%20subtype,%20taxid:114727&HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&ProtNames_ss=hemagglutinin&LabHost_s=include&SLen_i=550%20TO%20600&QualNum_i=0&CollectionDate_dr=2000-01-01T00:00:00.00Z%20TO%202022-05-01T23:59:59.00Z), H1N1 NA: [Link](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein&VirusLineage_ss=H1N1%20subtype,%20taxid:114727&HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&LabHost_s=include&QualNum_i=0&CollectionDate_dr=2000-01-01T00:00:00.00Z%20TO%202022-05-01T23:59:59.00Z&SLen_i=450%20TO%20500&ProtNames_ss=neuraminidase), H3N2 HA: [Link](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein&HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&LabHost_s=include&QualNum_i=0&CollectionDate_dr=2000-01-01T00:00:00.00Z%20TO%202022-05-01T23:59:59.00Z&VirusLineage_ss=H3N2%20subtype,%20taxid:119210&SLen_i=550%20TO%20650&ProtNames_ss=hemagglutinin), H3N2 NA: [Link](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/virus?SeqType_s=Protein&HostLineage_ss=Homo%20sapiens%20(human),%20taxid:9606&LabHost_s=include&QualNum_i=0&CollectionDate_dr=2000-01-01T00:00:00.00Z%20TO%202022-05-01T23:59:59.00Z&SLen_i=450%20TO%20500&ProtNames_ss=neuraminidase&VirusLineage_ss=H3N2%20subtype,%20taxid:119210)
    - Add data to `raw_data/ncbi` with the appropriate naming convention
        - EQUENCE_SEGMENT
        - ex. h1n1_ha.fasta
3. Duplicate any `influenza_qnet_predictions_<YEAR>.ipynb` and modify it for current season (more instructions in notebook)
    - Create new Qnets for each hemisphere, sub-type, and segment
    - Q-Centroid: $\widehat{x}^{t+1} = argmin_{x\in P} \sum_{y \in P^t} \theta(x,y)$
        - Where $P^t$ is the sequence population at time $t$ and $P = P^t \cup P^{t-1} \cup P^{t-2} \cup \dots \cup P^1$.
        - $\theta(x,y)$ is the qdistance between x and y in their respective Qnets
    - Multi-Cluster Predictions
        - Compute distance matrix between sequences in $P^t$
        - Create three clusters, then find the dominant strain of each cluster

## WHO vs. Qnet Predictions
        - See this [link](https://www.fludb.org/brc/vaccineRecommend.spg?decorator=influenza#:~:text=From%20these%20data%2C%20the%20WHO,recommendation%20are%20also%20usually%20suggested.) for WHO recommendations and this [link](https://platform.epicov.org/epi3/frontend#507f8c) to search the sequences
        - Compute edit distance between dominant strain and WHO predicted strain, Qnet predicted strains
        - Tables 4 - 15 in the paper


## File Tree
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
│   │   └── ncbi : all NCBI amino acid fasta data through May, 1 2022, by sub-type (H1N1/H3N2), segment (HA/NA)
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
