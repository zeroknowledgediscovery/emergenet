# Publication Draft

## Learning Mutational Patterns at Scale to Analyze Sequence Divergence in Novel Pathogens
Influenza viruses constantly evolve, and mismatches between predicted and circulating strains impact vaccine effectiveness. A barrier to predicting the season-specific dominant strains is the limited ability to predict future mutations, or estimate the numerical likelihood of specific future strains. In this study, we introduce a biology-aware sequence similarity metric based on deep pattern recognition of emergent evolutionary constraints. We use our model in two applications. One, we calculate the odds of future mutations, outperforming WHO recommended flu vaccine compositions almost consistently over the past two decades. Two, we compute emergence risk of strains previously analyzed by the CDC's Influenza Risk Assessment Tool, showing a moderately strong linear correlation between our predictions and the CDC's, though our predictions require much less time and resources.

## File Tree
```
EmergeNet
├── irat_qnet
├── qnet_predictions
└── tex
    ├── Figures
    │   ├── External : Fig. 1-3, 5, SI-Fig. 1-3
    │   ├── plotdata : CSV data for Fig. 4
    │   ├── seasonalpred_both.tex : LaTeX file for Fig. 4
    │   └── tabdata : LaTeX files for Tab. 1, SI-Tab. 1, 3-18
    ├── SI.pdf : supplementary text PDF
    ├── SI.tex : supplementary text body
    ├── SIfig.tex : supplementary text figures
    ├── allbib.bib : citations
    ├── authorpdf.pdf : main text + supplementary text PDF
    ├── customdef.tex : header definitions
    ├── main.pdf : main text PDF
    ├── main.tex : main text body
    ├── method.tex : main text computational framework
    └── preamble.tex : packages and styles
```
