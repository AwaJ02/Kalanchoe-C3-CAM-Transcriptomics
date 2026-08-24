# Kalanchoe-C3-CAM-Transcriptomics
R-based transcriptomic analysis of the developmental C3-to-CAM transition in Kalanchoë fedtschenkoi across the diel cycle. MSc Bioinformatics dissertation project
# Diel transcriptomic analysis of the developmental C3-to-CAM transition in Kalanchoë fedtschenkoi

## Overview

This repository contains the R scripts, processed results and figures
generated for an MSc Bioinformatics dissertation investigating
transcriptional changes associated with the developmental transition
from predominantly C3 photosynthesis in young LP1 leaves to established
CAM in mature LP6 leaves of Kalanchoë fedtschenkoi.

Gene expression was compared between LP1 and LP6 leaves at six diel
sampling times (10, 14, 18, 22, 02 and 06 h).

## Analysis

The analysis includes:

1. RNA-seq count-data processing
2. DESeq2 differential expression analysis
3. LP6 vs LP1 contrasts at each diel timepoint
4. DEG summarisation
5. KEGG pathway enrichment
6. Functional annotation of CAM-associated genes
7. PEPC, PPCK and PPDK analysis
8. Developmental stage × time interaction testing
9. Priority CAM candidate identification
10. Generation of dissertation figures

## Repository structure

- `scripts/` – R analysis scripts
- `data/` – input and processed datasets
- `results/` – statistical analysis outputs
- `figures/` – final dissertation figures
- `supplementary/` – additional supporting results

## Software

Analysis was performed in R using packages including DESeq2,
clusterProfiler, KEGGREST, dplyr, tidyr and ggplot2.
