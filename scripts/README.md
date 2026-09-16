# Analysis Scripts

This directory contains the R code used for the transcriptomic analysis of
*Kalanchoë fedtschenkoi*.

## Main analysis script

`01_complete_RNAseq_analysis.R`

The analysis includes:

- RNA-seq count-data preparation
- Sample metadata construction
- DESeq2 normalisation and differential-expression analysis
- LP6 versus LP1 comparisons across six diel timepoints
- Principal component analysis (PCA)
- Developmental stage × time interaction analysis
- EggNOG functional annotation
- KEGG pathway enrichment analysis
- KEGG metabolic pathway analysis
- Analysis of CAM-associated genes
- Candidate-gene prioritisation
- Generation of dissertation figures

Positive log2 fold-change values indicate higher expression in mature LP6
(CAM) leaves, while negative values indicate higher expression in young LP1
(predominantly C3) leaves.

The complete analysis workflow is provided in
`01_complete_RNAseq_analysis.R`.
