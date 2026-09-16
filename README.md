# Developmental and Diel Transcriptomic Reprogramming Associated
with the C3-to-CAM Transition in Kalanchoë fedtschenkoi

## Overview

This repository contains the R scripts, processed results and figures
used for my MSc Bioinformatics dissertation.

The project investigates developmental and diel transcriptional changes
associated with the transition from predominantly C3 photosynthesis in
young LP1 leaves to established Crassulacean Acid Metabolism (CAM) in
mature LP6 leaves of Kalanchoë fedtschenkoi.

## Experimental Design

RNA-seq data consisted of 36 samples:

- Two developmental stages: LP1 and LP6
- Six diel timepoints: 10, 14, 18, 22, 02 and 06 h
- Three biological replicates per developmental stage × timepoint

## Analysis

The bioinformatics workflow included:

1. RNA-seq count-data preparation
2. DESeq2 normalisation and differential-expression analysis
3. LP6 versus LP1 comparisons at each diel timepoint
4. Developmental stage × time interaction analysis
5. Functional gene annotation
6. KEGG pathway enrichment
7. KEGG metabolic pathway mapping
8. Analysis of CAM-associated genes
9. Candidate-gene prioritisation

Differentially expressed genes were defined using:

adjusted p-value < 0.05 and |log2 fold change| > 1.

## Repository Structure

`scripts/` – R scripts used for the analyses  
`data/` – metadata and information describing input datasets  
`results/` – processed analysis outputs and summary tables  
`figures/` – figures used in the dissertation  
`supplementary/` – additional pathway maps and supporting analyses

## Key Findings

- 7,033–10,743 genes were differentially expressed between LP1 and LP6,
  depending on diel time.
- 12,588 of 29,415 genes (42.79%) showed a significant developmental
  stage × time interaction.
- LP6-higher genes were enriched in pathways associated with carbon
  fixation, carbon metabolism, carbohydrate metabolism and energy
  metabolism.
- KfGene007051-RA (PPDK) showed consistently higher expression in LP6
  across all six diel timepoints.
- PEPC and PPCK paralogues displayed contrasting developmental and
  temporal expression patterns.

## Software

Analyses were performed in R using packages including:

- DESeq2
- clusterProfiler
- ggplot2

KEGG Mapper was used to visualise differential expression within
metabolic pathways.

## Author

Awa Jammeh  
MSc Bioinformatics  
University of Liverpool
