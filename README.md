# Developmental and Diel Transcriptomic Reprogramming Associated with the C3-to-CAM Transition in *Kalanchoë fedtschenkoi*

## Overview

This repository contains the R-based transcriptomic analysis performed for an
MSc Bioinformatics dissertation investigating developmental and diel
transcriptional reprogramming associated with the transition from C3-like
photosynthesis to Crassulacean acid metabolism (CAM) in
*Kalanchoë fedtschenkoi*.

The analysis compares young leaf pair 1 (LP1), which is predominantly C3,
with mature leaf pair 6 (LP6), which exhibits established CAM, across a
24-hour diel cycle.

## Experimental Design

The RNA-seq dataset consists of 36 samples:

- 2 developmental leaf stages: LP1 and LP6
- 6 sampling times: 10 h, 14 h, 18 h, 22 h, 02 h and 06 h
- 3 biological replicates for each leaf stage × timepoint combination

Plants were sampled across a 12 h light / 12 h dark cycle.

## Analysis Workflow

The analysis was performed primarily in R and included:

1. RNA-seq count-data preparation and sample metadata construction
2. DESeq2 normalisation and differential expression analysis
3. Principal component analysis (PCA)
4. LP6 versus LP1 differential expression at each diel timepoint
5. Developmental stage × time interaction analysis
6. Functional annotation using EggNOG-derived annotations
7. Directional KEGG pathway enrichment
8. KEGG Mapper visualisation of metabolic pathways
9. Analysis of CAM-associated genes
10. Candidate-gene prioritisation

## Differential Expression

Differentially expressed genes (DEGs) were defined using:

- adjusted p-value < 0.05
- |log2 fold change| > 1

For LP6 versus LP1 contrasts:

- positive log2 fold change = higher expression in LP6 (mature CAM)
- negative log2 fold change = higher expression in LP1 (young, predominantly C3)

Substantial transcriptional differences were detected throughout the diel
cycle, with the largest number of DEGs observed at 22 h.

A likelihood-ratio test was also used to identify genes whose developmental
expression differences varied across the diel cycle.

## Functional and Pathway Analysis

Functional annotation and KEGG Orthology (KO) assignments were used to
investigate the biological processes associated with the developmental
C3-to-CAM transition.

Directional KEGG enrichment was performed separately for genes showing
higher expression in LP6 and genes showing higher expression in LP1.

The analysis identified coordinated differences in pathways including carbon
fixation, carbon metabolism, glycolysis/gluconeogenesis, starch and sucrose
metabolism, oxidative phosphorylation and plant hormone signalling.

KEGG Mapper was additionally used to visualise transcriptional differences
within metabolic pathways.

### KEGG Mapper Colour Convention

For KEGG Mapper pathway images:

- blue = higher expression in LP6 (CAM)
- red = higher expression in LP1 (C3-like)
- uncoloured = no significant differentially expressed mapped KO

KEGG Mapper visualisations represent transcript-level differences and should
not be interpreted as direct measurements of enzyme activity or metabolic flux.

## CAM-Associated Gene Analysis

CAM-associated genes were examined across the diel cycle to identify
developmental and temporal expression patterns.

Particular attention was given to genes associated with PEPC, PPCK and PPDK.
The analysis demonstrated that CAM-associated transcriptional changes were
not uniformly directed towards LP6, with different paralogues displaying
distinct developmental and diel expression patterns.

These findings support a network-level and temporally regulated model of CAM
development rather than a simple uniform increase in individual CAM genes.

## Repository Structure

```text
Kalanchoe-C3-CAM-Transcriptomics/
├── README.md
├── .gitignore
├── data/
│   ├── README.md
│   └── sample metadata
├── scripts/
│   ├── README.md
│   └── 01_complete_RNAseq_analysis.R
├── results/
│   ├── README.md
│   ├── differential expression summaries
│   ├── interaction results
│   ├── annotation summaries
│   ├── CAM candidate results
│   └── KEGG_enrichment/
├── figures/
│   ├── README.md
│   └── final dissertation figures
└── supplementary/
    ├── CAM_analysis/
    └── KEGG_Mapper/
