# Analysis Scripts

This directory contains the R code used for the transcriptomic analysis of
*Kalanchoë fedtschenkoi*.

The analysis includes:

- RNA-seq count-data preparation and sample metadata construction
- DESeq2 normalisation and differential-expression analysis
- LP6 versus LP1 comparisons across six diel timepoints
- Principal component analysis (PCA)
- Developmental stage × time interaction analysis
- EggNOG functional annotation
- KEGG pathway enrichment analysis
- KEGG metabolic pathway mapping
- Analysis of CAM-associated genes
- Candidate-gene prioritisation
- Generation of dissertation figures

The primary analysis was conducted in R. Positive log2 fold-change values
represent higher expression in mature LP6 leaves, while negative values
represent higher expression in young LP1 leaves.
