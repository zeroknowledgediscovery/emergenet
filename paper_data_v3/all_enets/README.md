# Emergenet Models

This directory contains gzipped versions of the Emergenet models trained for the paper *Emergenet: Fast Scalable Emergence Risk Assessment of Influenza A Strains Circulating In Non-human Hosts*.

## `enet_predictions`: Predicting Seasonal Vaccine Strains

Prediction of the seasonal vaccine strain as a close match to a historical strain allows out-of-sample validation against past World Health Organization (WHO) recommendations for the flu shot, which is reformulated February and September of each year for the northern and southern hemispheres, respectively, based on a cocktail of historical strains determined via global surveillance. For each year of the past two decades (2002-03 to 2022-23), we constructed a separate Emergenet for the southern and northern flu seasons using HA strains available from the previous season, for H1N1 and H3N2 subtypes. For seasons with > 3000 strains available, we randomly sampled 3000 strains which provide an accurate representation of the population.


## `irat_enet`: Estimating Emergence Risk of Non-human Strains

Our primary claim is the ability to estimate the emergence potential of animal strains that do not yet circulate in humans. We validate this claim by training Emergenet models for strains analyzed by the CDC's [Influenza Risk Assessment Tool (IRAT)](https://www.cdc.gov/flu/pandemic-resources/national-strategy/risk-assessment.htm) and variant strains (animal strains that emerged in humans).

### `irat_enet/irat`

We constructed Emergenet models for HA and NA sequences using subtype-specific human strains, collected within the year prior to the assessment date, e.g., the assessment date for A/swine/Shandong/1207/2016 is 06/2020, and we use all human strains collected between 01/07/2019 and 06/30/2020 for the Emergenet inference. More specifically, for each HA subtype (H1Nx, H3Nx, etc) with > 15 strains collected in the past year, we construct an Emergenet model. 

### `irat_enet/variant`

A second demonstration of our approach is obtained by computing the E-risk scores for "variants" (animal Influenza A strains isolated in human hosts), which might be expected to pose high risk due to their successful replication in human cells, even if the possibility of HH transmission is not yet observed or guaranteed. We analyze 15 H1N1 and H1N2 variants isolated post-2015, 14 of which are listed as antigenic prototypes for canidate vaccine viruses by the CDC. We construct Enet models for each variant at six, three, one, and zero months before their collection date.

## File Tree

``` txt
all_enets
├── enet_predictions
│   ├── north_h1n1
│   ├── north_h3n2
│   ├── south_h1n1
│   └── south_h3n2
└── irat_enet
    ├── irat
    └── variant
```
