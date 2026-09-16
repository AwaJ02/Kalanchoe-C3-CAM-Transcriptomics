
# LIFE703 MSc Bioinformatics Project
# Corrected RNA-seq analysis of Kalanchoe fedtschenkoi
# LP1 (young C3 leaves) vs LP6 (mature full-CAM leaves)
# across six 4-hourly timepoints


# GAI DECLARATION
# I, [201945294], used ChatGPT (OpenAI) to assist with structuring and
# refining R code, debugging errors, and explaining analytical steps.
# All analyses, outputs and interpretations were reviewed, checked and
# modified where necessary for correctness and suitability for this project.



# ANALYSIS OVERVIEW

# Experimental design confirmed by the data provider:
#
#   LP1 (L1) = youngest leaf pair flanking the shoot apical meristem;
#              predominantly C3 photosynthesis
#
#   LP6 (L6) = fully developed/expanded leaf pair;
#              fully developed CAM
#
#   6 timepoints sampled every 4 hours across a 24-hour
#   12 h light / 12 h dark cycle
#
#   3 biological replicates per leaf pair per timepoint
#
#   2 leaf ages x 6 timepoints x 3 replicates = 36 samples
#
# IMPORTANT:
# The columns in readscount.tab.txt follow the original CGR sample order:
#
# samples 01-03 = LP1 at 10:00
# samples 04-06 = LP6 at 10:00
# samples 07-09 = LP1 at 14:00
# samples 10-12 = LP6 at 14:00
# samples 13-15 = LP1 at 18:00
# samples 16-18 = LP6 at 18:00
# samples 19-21 = LP1 at 22:00
# samples 22-24 = LP6 at 22:00
# samples 25-27 = LP1 at 02:00
# samples 28-30 = LP6 at 02:00
# samples 31-33 = LP1 at 06:00
# samples 34-36 = LP6 at 06:00
#
# The primary biological comparison is LP6 vs LP1 WITHIN each timepoint.
#
# For all per-timepoint contrasts below:
#   positive log2FoldChange = higher expression in LP6 (mature CAM)
#   negative log2FoldChange = higher expression in LP1 (young C3)


# 1. LOAD PACKAGES
# ========================================
library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidyr)

library(clusterProfiler)
library(enrichplot)


# ============================================================
# 2. IMPORT RNA-SEQ COUNT MATRIX
# ============================================================

counts <- read.table(
  "readscount.tab.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

cat("Dimensions of full supplied table:\n")
print(dim(counts))

# The first 36 columns are raw read counts.
# Columns after these are FPKM values and must NOT be used as DESeq2 input.

counts_raw <- counts[, 1:36]

cat("\nDimensions of raw count matrix:\n")
print(dim(counts_raw))

cat("\nRaw-count sample names:\n")
print(colnames(counts_raw))

# DESeq2 requires integer counts.
counts_raw <- round(as.matrix(counts_raw))
storage.mode(counts_raw) <- "integer"


# ============================================================
# 3. BUILD CORRECT SAMPLE METADATA
# ============================================================

time_order <- c("10", "14", "18", "22", "02", "06")

# Within EVERY timepoint:
# first 3 samples = LP1
# next 3 samples  = LP6

leaf <- rep(
  c(rep("LP1", 3), rep("LP6", 3)),
  times = 6
)

timepoint <- rep(
  time_order,
  each = 6
)

replicate <- rep(
  1:3,
  times = 12
)

sample_info <- data.frame(
  leaf = factor(leaf, levels = c("LP1", "LP6")),
  timepoint = factor(timepoint, levels = time_order),
  replicate = factor(replicate),
  stringsAsFactors = FALSE
)

rownames(sample_info) <- colnames(counts_raw)

# One 12-level group factor 
# This makes LP6-vs-LP1 contrasts at each timepoint direct and unambiguous.

group_levels <- as.vector(
  rbind(
    paste0("LP1_", time_order),
    paste0("LP6_", time_order)
  )
)

sample_info$group <- factor(
  paste(sample_info$leaf, sample_info$timepoint, sep = "_"),
  levels = group_levels
)

cat("\nCorrected sample metadata:\n")
print(sample_info)

cat("\nSamples per leaf/timepoint group:\n")
print(table(sample_info$leaf, sample_info$timepoint))

# Sanity checks
stopifnot(nrow(sample_info) == 36)
stopifnot(ncol(counts_raw) == 36)
stopifnot(identical(rownames(sample_info), colnames(counts_raw)))
stopifnot(all(table(sample_info$group) == 3))


# Save metadata for dissertation/reproducibility

write.csv(
  sample_info,
  "Corrected_sample_metadata_LP1_LP6.csv",
  row.names = TRUE
)


# ============================================================
# 4.DESEQ2 DATASET
# ============================================================

# A group-based design is used because the main question is:
# LP6 vs LP1 at EACH of the six timepoints.
#
# This design fits one mean for each of the 12 leaf-age/time combinations
# and avoids confusing an interaction coefficient with a within-timepoint
# LP6-vs-LP1 comparison.

dds <- DESeqDataSetFromMatrix(
  countData = counts_raw,
  colData = sample_info,
  design = ~ group
)


# ============================================================
# 5. FILTER LOW-COUNT GENES
# ============================================================

# Retain genes with >10 reads in total across all 36 libraries.

keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]

cat("\nGenes retained after filtering:", nrow(dds), "\n")


# ============================================================
# 6. RUN DESEQ2
# ============================================================

dds <- DESeq(dds)

cat("\nDESeq2 coefficient names:\n")
print(resultsNames(dds))


# ============================================================
# 7. PCA: QUALITY CONTROL AND GLOBAL EXPRESSION STRUCTURE
# ============================================================

vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(
  vsd,
  intgroup = c("leaf", "timepoint"),
  returnData = TRUE
)

percent_var <- round(
  100 * attr(pca_data, "percentVar"),
  1
)

theme_dissertation <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    panel.grid.minor = element_blank(),
    legend.title = element_text(face = "bold")
  )

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    colour = leaf,
    shape = timepoint
  )
) +
  geom_point(size = 4, alpha = 0.9) +
  labs(
    title = "PCA of Kalanchoe RNA-seq samples",
    subtitle = "LP1 (young C3) vs LP6 (mature CAM) across the diel cycle",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    colour = "Leaf pair",
    shape = "Time"
  ) +
  theme_dissertation

print(pca_plot)

ggsave(
  "PCA_LP1_LP6_timecourse.png",
  pca_plot,
  width = 9,
  height = 7,
  dpi = 300
)

ggsave(
  "PCA_LP1_LP6_timecourse.pdf",
  pca_plot,
  width = 9,
  height = 7
)


# ============================================================
# 8. FUNCTION TO EXTRACT LP6 VS LP1 AT ONE TIMEPOINT
# ============================================================

get_leaf_contrast <- function(dds_object, tp) {
  
  lp6_group <- paste0("LP6_", tp)
  lp1_group <- paste0("LP1_", tp)
  
  res <- results(
    dds_object,
    contrast = c(
      "group",
      lp6_group,
      lp1_group
    ),
    alpha = 0.05
  )
  
  res_df <- as.data.frame(res)
  res_df$GeneID <- rownames(res_df)
  res_df$Timepoint <- tp
  
  # Biological direction:
  # positive = LP6 > LP1
  # negative = LP1 > LP6
  
  res_df$Direction <- case_when(
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      res_df$log2FoldChange > 1 ~ "Higher_in_LP6_CAM",
    
    !is.na(res_df$padj) &
      res_df$padj < 0.05 &
      res_df$log2FoldChange < -1 ~ "Higher_in_LP1_C3",
    
    TRUE ~ "Not_significant"
  )
  
  res_df
}


# ============================================================
# 9. RUN ALL SIX LP6 VS LP1 CONTRASTS
# ============================================================

dir.create(
  "DESeq2_corrected_results",
  showWarnings = FALSE
)

all_de_results <- list()
de_summary <- data.frame()

for (tp in time_order) {
  
  cat(
    "\n====================================\n",
    "LP6 vs LP1 at ", tp, ":00\n",
    "====================================\n",
    sep = ""
  )
  
  res_df <- get_leaf_contrast(dds, tp)
  
  all_de_results[[tp]] <- res_df
  
  # Save ALL tested genes
  write.csv(
    res_df,
    file.path(
      "DESeq2_corrected_results",
      paste0("DESeq2_LP6_vs_LP1_", tp, "h_all_genes.csv")
    ),
    row.names = FALSE
  )
  
  # Significant genes:
  # adjusted p-value < 0.05
  # absolute log2 fold change > 1
  
  sig <- res_df %>%
    filter(
      !is.na(padj),
      padj < 0.05,
      abs(log2FoldChange) > 1
    )
  
  write.csv(
    sig,
    file.path(
      "DESeq2_corrected_results",
      paste0("Significant_LP6_vs_LP1_", tp, "h.csv")
    ),
    row.names = FALSE
  )
  
  n_up_lp6 <- sum(sig$log2FoldChange > 0)
  n_up_lp1 <- sum(sig$log2FoldChange < 0)
  
  de_summary <- bind_rows(
    de_summary,
    data.frame(
      Timepoint = tp,
      Significant_total = nrow(sig),
      Higher_in_LP6_CAM = n_up_lp6,
      Higher_in_LP1_C3 = n_up_lp1
    )
  )
  
  cat("Significant genes:", nrow(sig), "\n")
  cat("Higher in LP6 (mature CAM):", n_up_lp6, "\n")
  cat("Higher in LP1 (young C3):", n_up_lp1, "\n")
}

write.csv(
  de_summary,
  "DESeq2_corrected_results/DE_gene_counts_by_timepoint.csv",
  row.names = FALSE
)

print(de_summary)


# ============================================================
# 10. PLOT NUMBER OF DIFFERENTIALLY EXPRESSED GENES OVER TIME
# ============================================================

de_summary_long <- de_summary %>%
  select(
    Timepoint,
    Higher_in_LP6_CAM,
    Higher_in_LP1_C3
  ) %>%
  pivot_longer(
    cols = -Timepoint,
    names_to = "Direction",
    values_to = "Genes"
  )

de_summary_long$Timepoint <- factor(
  de_summary_long$Timepoint,
  levels = time_order
)

de_count_plot <- ggplot(
  de_summary_long,
  aes(
    x = Timepoint,
    y = Genes,
    fill = Direction
  )
) +
  geom_col(position = "dodge") +
  labs(
    title = "Differentially expressed genes across the diel cycle",
    subtitle = "LP6 (mature CAM) vs LP1 (young C3)",
    x = "Sampling time",
    y = "Number of significant genes",
    fill = "Expression direction"
  ) +
  theme_dissertation

print(de_count_plot)

ggsave(
  "DESeq2_corrected_results/DE_gene_counts_by_timepoint.png",
  de_count_plot,
  width = 9,
  height = 6,
  dpi = 300
)


# ============================================================
# 11. MA AND VOLCANO PLOTS FOR EACH TIMEPOINT
# ============================================================

dir.create(
  "DESeq2_corrected_results/plots",
  showWarnings = FALSE
)

for (tp in time_order) {
  
  res_df <- all_de_results[[tp]]
  
  # ---------------- MA plot ----------------
  
  ma_plot <- ggplot(
    res_df,
    aes(
      x = baseMean,
      y = log2FoldChange,
      colour = Direction
    )
  ) +
    geom_point(alpha = 0.55, size = 1.1) +
    scale_x_log10() +
    geom_hline(
      yintercept = 0,
      linewidth = 0.7
    ) +
    labs(
      title = paste0("MA plot: LP6 vs LP1 at ", tp, ":00"),
      x = "Mean normalised count (log10 scale)",
      y = "Log2 fold change (LP6 / LP1)",
      colour = NULL
    ) +
    theme_dissertation
  
  ggsave(
    file.path(
      "DESeq2_corrected_results/plots",
      paste0("MA_LP6_vs_LP1_", tp, "h.png")
    ),
    ma_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
  
  # ---------------- Volcano plot ----------------
  
  # Use adjusted p-value for both significance classification
  # and the plotted significance axis for consistency.
  
  volcano_df <- res_df %>%
    mutate(
      negLog10Padj = -log10(
        pmax(padj, .Machine$double.xmin)
      )
    )
  
  volcano_plot <- ggplot(
    volcano_df,
    aes(
      x = log2FoldChange,
      y = negLog10Padj,
      colour = Direction
    )
  ) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(
      xintercept = c(-1, 1),
      linetype = "dashed"
    ) +
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dashed"
    ) +
    labs(
      title = paste0("Volcano plot: LP6 vs LP1 at ", tp, ":00"),
      x = "Log2 fold change (LP6 / LP1)",
      y = expression(-log[10]("adjusted p-value")),
      colour = NULL
    ) +
    theme_dissertation
  
  ggsave(
    file.path(
      "DESeq2_corrected_results/plots",
      paste0("Volcano_LP6_vs_LP1_", tp, "h.png")
    ),
    volcano_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}


# ============================================================
# 12. OPTIONAL: GLOBAL TEST FOR LEAF-AGE x TIME INTERACTION
# ============================================================

# This asks:
# "Which genes show a DIFFERENT LP6-vs-LP1 effect depending on time?"
#
# It is separate from the six direct LP6-vs-LP1 contrasts above.

dds_interaction <- DESeqDataSetFromMatrix(
  countData = counts_raw,
  colData = sample_info,
  design = ~ leaf + timepoint + leaf:timepoint
)

dds_interaction <- dds_interaction[
  rowSums(counts(dds_interaction)) > 10,
]

dds_interaction <- DESeq(
  dds_interaction,
  test = "LRT",
  reduced = ~ leaf + timepoint
)

interaction_res <- results(
  dds_interaction,
  alpha = 0.05
)

interaction_df <- as.data.frame(interaction_res)
interaction_df$GeneID <- rownames(interaction_df)

interaction_sig <- interaction_df %>%
  filter(
    !is.na(padj),
    padj < 0.05
  )

write.csv(
  interaction_df,
  "DESeq2_corrected_results/Leaf_by_time_interaction_all_genes.csv",
  row.names = FALSE
)

write.csv(
  interaction_sig,
  "DESeq2_corrected_results/Leaf_by_time_interaction_significant.csv",
  row.names = FALSE
)

cat(
  "\nGenes with significant leaf-age x time interaction:",
  nrow(interaction_sig),
  "\n"
)


# ============================================================
# 13. READING AND CLEANING EGGNOG ANNOTATIONS
# ============================================================

# The .xls file is actually tab-delimited text.
#
# Use whichever copy exists.

annotation_candidates <- c(
  "Kalanchoe.fedtschenkoi.eggnog.annotations(1).xls",
  "Kalanchoe.fedtschenkoi.eggnog.annotations.xls"
)

annotation_file <- annotation_candidates[
  file.exists(annotation_candidates)
][1]

if (is.na(annotation_file)) {
  stop("EggNOG annotation file was not found.")
}

anno <- read.delim(
  annotation_file,
  header = TRUE,
  sep = "\t",
  skip = 1,
  comment.char = "",
  quote = "",
  fill = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

colnames(anno)[1] <- "GeneID"

required_annotation_cols <- c(
  "GeneID",
  "predicted_gene_name",
  "KEGG_KOs",
  "eggNOG annot"
)

missing_annotation_cols <- setdiff(
  required_annotation_cols,
  colnames(anno)
)

if (length(missing_annotation_cols) > 0) {
  stop(
    "Missing annotation columns: ",
    paste(missing_annotation_cols, collapse = ", ")
  )
}

anno_clean <- anno %>%
  transmute(
    GeneID = trimws(as.character(GeneID)),
    Gene_Name = predicted_gene_name,
    KO = KEGG_KOs,
    Description = `eggNOG annot`
  ) %>%
  distinct()

# Converting blank/text missing values to NA

for (column_name in c("Gene_Name", "KO", "Description")) {
  anno_clean[[column_name]][
    anno_clean[[column_name]] %in%
      c("", "NA", "N/A", "-", "nan")
  ] <- NA
}

cat(
  "\nAnnotation rows:",
  nrow(anno_clean),
  "\n"
)


# ============================================================
# 14. ANNOTATING SIGNIFICANT GENES + TOP 15 PER TIMEPOINT
# ============================================================

dir.create(
  "Annotated_corrected_results",
  showWarnings = FALSE
)

all_annotated_sig <- data.frame()
candidate_summary <- data.frame()

for (tp in time_order) {
  
  sig <- all_de_results[[tp]] %>%
    filter(
      !is.na(padj),
      padj < 0.05,
      abs(log2FoldChange) > 1
    )
  
  sig_annotated <- sig %>%
    left_join(
      anno_clean,
      by = "GeneID"
    )
  
  write.csv(
    sig_annotated,
    file.path(
      "Annotated_corrected_results",
      paste0("Significant_LP6_vs_LP1_", tp, "h_annotated.csv")
    ),
    row.names = FALSE
  )
  
  # Top 15 higher in LP6/CAM
  top_lp6 <- sig_annotated %>%
    filter(log2FoldChange > 0) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 15) %>%
    mutate(
      Biological_direction = "Higher in LP6 (mature CAM)"
    )
  
  # Top 15 higher in LP1/C3
  top_lp1 <- sig_annotated %>%
    filter(log2FoldChange < 0) %>%
    arrange(padj, desc(abs(log2FoldChange))) %>%
    slice_head(n = 15) %>%
    mutate(
      Biological_direction = "Higher in LP1 (young C3)"
    )
  
  write.csv(
    top_lp6,
    file.path(
      "Annotated_corrected_results",
      paste0("Top15_higher_in_LP6_CAM_", tp, "h.csv")
    ),
    row.names = FALSE
  )
  
  write.csv(
    top_lp1,
    file.path(
      "Annotated_corrected_results",
      paste0("Top15_higher_in_LP1_C3_", tp, "h.csv")
    ),
    row.names = FALSE
  )
  
  candidate_summary <- bind_rows(
    candidate_summary,
    top_lp6,
    top_lp1
  )
  
  all_annotated_sig <- bind_rows(
    all_annotated_sig,
    sig_annotated
  )
}

write.csv(
  candidate_summary,
  "Annotated_corrected_results/Top_candidate_gene_summary_all_timepoints.csv",
  row.names = FALSE
)

write.csv(
  all_annotated_sig,
  "Annotated_corrected_results/All_significant_annotated_all_timepoints.csv",
  row.names = FALSE
)


# ============================================================
# 15. ANNOTATION COVERAGE SUMMARY
# ============================================================

annotation_summary <- all_annotated_sig %>%
  group_by(Timepoint) %>%
  summarise(
    Significant_genes = n(),
    Genes_with_description = sum(
      !is.na(Description) & Description != ""
    ),
    Genes_with_KO = sum(
      !is.na(KO) & KO != ""
    ),
    Percentage_with_description = round(
      100 * Genes_with_description / Significant_genes,
      2
    ),
    Percentage_with_KO = round(
      100 * Genes_with_KO / Significant_genes,
      2
    ),
    .groups = "drop"
  )

write.csv(
  annotation_summary,
  "Annotated_corrected_results/Annotation_coverage_summary.csv",
  row.names = FALSE
)

print(annotation_summary)


# ============================================================
# 16. SEARCHING FOR CAM-RELATED CANDIDATE GENES
# ============================================================

# This is a hypothesis-guided screen, not proof that a gene is CAM-specific.

cam_keywords <- c(
  "phosphoenolpyruvate",
  "PEPC",
  "PEP carboxylase",
  "PPCK",
  "malate",
  "malic enzyme",
  "malate dehydrogenase",
  "MDH",
  "pyruvate",
  "PPDK",
  "carbonic anhydrase",
  "starch",
  "sucrose",
  "vacuolar",
  "tonoplast",
  "aquaporin",
  "circadian",
  "clock",
  "stomata",
  "guard cell"
)

cam_pattern <- paste(
  cam_keywords,
  collapse = "|"
)

cam_candidates <- all_annotated_sig %>%
  filter(
    grepl(
      cam_pattern,
      paste(
        Gene_Name,
        Description
      ),
      ignore.case = TRUE
    )
  ) %>%
  arrange(
    Timepoint,
    padj
  )

write.csv(
  cam_candidates,
  "Annotated_corrected_results/CAM_candidate_genes_corrected.csv",
  row.names = FALSE
)


# ============================================================
# 17. PREPARE KALANCHOE GENE -> KO -> KEGG PATHWAY MAPPING
# ============================================================
# ============================================================
# 17.1 CHECK REQUIRED OBJECTS AND FILES
# ============================================================

if (!exists("anno_clean")) {
  stop(
    "anno_clean does not exist. Run the EggNOG annotation section first."
  )
}

if (!exists("dds")) {
  stop(
    "dds does not exist. Run the DESeq2 analysis first."
  )
}

if (!exists("all_de_results")) {
  stop(
    "all_de_results does not exist. Run the six LP6 vs LP1 contrasts first."
  )
}

if (!exists("time_order")) {
  stop(
    "time_order does not exist. Run the sample metadata section first."
  )
}

if (!file.exists("gene2kegg.tab.txt")) {
  stop(
    "gene2kegg.tab.txt was not found in the working directory."
  )
}

# ============================================================
# 17.2 CLEAN EGGNOG KO ANNOTATIONS
# ============================================================

# Check that KO annotations exist

cat("\nExample KO annotations from EggNOG:\n")

print(
  head(
    anno_clean$KO[
      !is.na(anno_clean$KO) &
        anno_clean$KO != ""
    ],
    20
  )
)


cat(
  "\nGenes with a non-missing KO field:",
  sum(
    !is.na(anno_clean$KO) &
      anno_clean$KO != ""
  ),
  "\n"
)


# Some genes contain multiple KO identifiers, for example:
# K01080,K18693
# These need to be split into separate rows.

gene_to_ko <- anno_clean %>%
  
  select(
    GeneID,
    KO
  ) %>%
  
  filter(
    !is.na(KO),
    KO != ""
  ) %>%
  
  tidyr::separate_rows(
    KO,
    sep = "[,;]"
  ) %>%
  
  mutate(
    
    GeneID = trimws(
      as.character(GeneID)
    ),
    
    KO = trimws(
      as.character(KO)
    ),
    
    # Remove possible prefixes
    KO = sub(
      "^ko:",
      "",
      KO
    )
  ) %>%
  
  filter(
    grepl(
      "^K[0-9]{5}$",
      KO
    )
  ) %>%
  
  distinct()


cat(
  "\nKalanchoe genes with valid KO assignments:",
  length(
    unique(gene_to_ko$GeneID)
  ),
  "\n"
)


cat(
  "Unique KO identifiers in EggNOG:",
  length(
    unique(gene_to_ko$KO)
  ),
  "\n"
)


cat(
  "\nExample Kalanchoe gene -> KO mappings:\n"
)

print(
  head(
    gene_to_ko,
    10
  )
)


# ============================================================
# 17.3 READ KO -> KEGG PATHWAY FILE
# ============================================================

# IMPORTANT:
#
# gene2kegg.tab.txt contains values such as:
#
# K00001    00010
#
# The pathway ID MUST remain "00010".
#
# Therefore colClasses = "character" is essential.
# Otherwise R may convert 00010 -> 10.

ko_to_pathway <- read.delim(
  
  "gene2kegg.tab.txt",
  
  header = FALSE,
  
  sep = "\t",
  
  quote = "\"",
  
  stringsAsFactors = FALSE,
  
  colClasses = "character",
  
  fill = TRUE,
  
  check.names = FALSE
)


# We only need the first two columns

if (ncol(ko_to_pathway) < 2) {
  
  stop(
    "gene2kegg.tab.txt contains fewer than two columns."
  )
}


ko_to_pathway <- ko_to_pathway[, 1:2]


colnames(
  ko_to_pathway
) <- c(
  "KO",
  "Pathway"
)


# ============================================================
# 17.4 CLEAN KO -> PATHWAY MAPPING
# ============================================================

ko_to_pathway <- ko_to_pathway %>%
  
  mutate(
    
    KO = trimws(
      as.character(KO)
    ),
    
    Pathway = trimws(
      as.character(Pathway)
    ),
    
    # Remove possible prefixes
    KO = sub(
      "^ko:",
      "",
      KO
    ),
    
    Pathway = sub(
      "^path:",
      "",
      Pathway
    ),
    
    Pathway = sub(
      "^ko",
      "",
      Pathway
    )
  ) %>%
  
  filter(
    
    # KO format e.g. K00001
    grepl(
      "^K[0-9]{5}$",
      KO
    ),
    
    # KEGG pathway format e.g. 00010
    grepl(
      "^[0-9]{5}$",
      Pathway
    )
  ) %>%
  
  distinct()


cat(
  "\nKO identifiers in gene2kegg file:",
  length(
    unique(ko_to_pathway$KO)
  ),
  "\n"
)


cat(
  "KEGG pathways represented in mapping file:",
  length(
    unique(ko_to_pathway$Pathway)
  ),
  "\n"
)


cat(
  "\nExample KO -> pathway mappings:\n"
)

print(
  head(
    ko_to_pathway,
    10
  )
)


# ============================================================
# 17.5 CHECK KO OVERLAP
# ============================================================

shared_kos <- base::intersect(
  
  unique(
    gene_to_ko$KO
  ),
  
  unique(
    ko_to_pathway$KO
  )
)


cat(
  "\nShared KO identifiers between EggNOG and gene2kegg:",
  length(shared_kos),
  "\n"
)


cat(
  "\nExamples of shared KO identifiers:\n"
)

print(
  head(
    shared_kos,
    20
  )
)


# This should NOT be zero.

if (length(shared_kos) == 0) {
  
  cat(
    "\nFirst EggNOG KO values:\n"
  )
  
  print(
    head(
      unique(gene_to_ko$KO),
      20
    )
  )
  
  
  cat(
    "\nFirst gene2kegg KO values:\n"
  )
  
  print(
    head(
      unique(ko_to_pathway$KO),
      20
    )
  )
  
  
  stop(
    paste(
      "No KO overlap was found.",
      "Do not continue to enrichment.",
      "Inspect the KO formats printed above."
    )
  )
}


# ============================================================
# 17.6 MAP KALANCHOE GENES TO KEGG PATHWAYS
# ============================================================

# Join:
#
# Kalanchoe GeneID
#       ↓
#      KO
#       ↓
# KEGG pathway

gene_to_pathway <- gene_to_ko %>%
  
  inner_join(
    ko_to_pathway,
    by = "KO",
    relationship = "many-to-many"
  ) %>%
  
  select(
    GeneID,
    KO,
    Pathway
  ) %>%
  
  distinct()


cat(
  "\nKalanchoe genes mapped to at least one KEGG pathway:",
  length(
    unique(
      gene_to_pathway$GeneID
    )
  ),
  "\n"
)


cat(
  "Unique KEGG pathways represented:",
  length(
    unique(
      gene_to_pathway$Pathway
    )
  ),
  "\n"
)


cat(
  "Total gene-pathway relationships:",
  nrow(
    gene_to_pathway
  ),
  "\n"
)


cat(
  "\nExample Kalanchoe gene -> KO -> pathway mappings:\n"
)

print(
  head(
    gene_to_pathway,
    20
  )
)


# ============================================================
# 17.7 CREATE TERM2GENE TABLE
# ============================================================

# clusterProfiler::enricher() requires:
#
# column 1 = pathway/term
# column 2 = gene

term2gene <- gene_to_pathway %>%
  
  select(
    Pathway,
    GeneID
  ) %>%
  
  distinct()


colnames(
  term2gene
) <- c(
  "term",
  "gene"
)


cat(
  "\nTERM2GENE dimensions:\n"
)

print(
  dim(term2gene)
)


# ============================================================
# 17.8 CREATE KEGG RESULTS DIRECTORY
# ============================================================

dir.create(
  "KEGG_corrected_results",
  showWarnings = FALSE
)


# Save mappings for reproducibility

write.csv(
  
  gene_to_ko,
  
  "KEGG_corrected_results/Kalanchoe_gene_to_KO_mapping.csv",
  
  row.names = FALSE
)


write.csv(
  
  gene_to_pathway,
  
  "KEGG_corrected_results/Kalanchoe_gene_KO_pathway_mapping.csv",
  
  row.names = FALSE
)


# ============================================================
# 18. CREATEING KEGG ENRICHMENT BACKGROUND
# ============================================================

# Background/universe:
#
# All genes that:
#
# 1. survived DESeq2 filtering
# 2. were therefore tested
# 3. have a valid KEGG pathway mapping
#
# This is more appropriate than using every gene in KEGG.

tested_gene_ids <- rownames(dds)


enrichment_universe <- base::intersect(
  
  tested_gene_ids,
  
  unique(
    term2gene$gene
  )
)


cat(
  "\nGenes tested by DESeq2:",
  length(
    tested_gene_ids
  ),
  "\n"
)


cat(
  "Genes tested by DESeq2 with KEGG pathway mappings:",
  length(
    enrichment_universe
  ),
  "\n"
)


if (length(enrichment_universe) == 0) {
  
  stop(
    paste(
      "The KEGG enrichment universe contains zero genes.",
      "Do not continue until the GeneID mapping is checked."
    )
  )
}


# ============================================================
# 18.1 KEGG ENRICHMENT FUNCTION
# ============================================================

run_pathway_enrichment <- function(
    gene_ids,
    tp,
    direction_label
) {
  
  
  # ----------------------------------------------------------
  # Keep genes represented in the KEGG background
  # ----------------------------------------------------------
  
  gene_ids <- unique(
    
    base::intersect(
      
      gene_ids,
      
      enrichment_universe
    )
  )
  
  
  cat(
    "\n",
    tp,
    ":00 | ",
    direction_label,
    " | ",
    length(gene_ids),
    " pathway-mapped genes",
    "\n",
    sep = ""
  )
  
  
  # ----------------------------------------------------------
  # Require enough genes for meaningful enrichment
  # ----------------------------------------------------------
  
  if (length(gene_ids) < 10) {
    
    warning(
      
      paste0(
        
        "Fewer than 10 pathway-mapped genes for ",
        
        tp,
        
        ":00 ",
        
        direction_label,
        
        ". Enrichment skipped."
      )
    )
    
    
    return(NULL)
  }
  
  
  # ----------------------------------------------------------
  # Over-representation analysis
  # ----------------------------------------------------------
  
  enrichment <- clusterProfiler::enricher(
    
    gene = gene_ids,
    
    universe = enrichment_universe,
    
    TERM2GENE = term2gene,
    
    pvalueCutoff = 0.05,
    
    pAdjustMethod = "BH",
    
    qvalueCutoff = 0.20,
    
    minGSSize = 10,
    
    maxGSSize = 1000
  )
  
  
  enrichment_df <- as.data.frame(
    enrichment
  )
  
  
  # ----------------------------------------------------------
  # Add biological information
  # ----------------------------------------------------------
  
  if (nrow(enrichment_df) > 0) {
    
    enrichment_df$Timepoint <- tp
    
    enrichment_df$Direction <- direction_label
  }
  
  
  # ----------------------------------------------------------
  # Save table
  # ----------------------------------------------------------
  
  write.csv(
    
    enrichment_df,
    
    file.path(
      
      "KEGG_corrected_results",
      
      paste0(
        
        "KEGG_",
        
        tp,
        
        "h_",
        
        direction_label,
        
        ".csv"
      )
    ),
    
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Produce dotplot
  # ----------------------------------------------------------
  
  if (nrow(enrichment_df) > 0) {
    
    
    pathway_plot <- enrichplot::dotplot(
      
      enrichment,
      
      showCategory = min(
        15,
        nrow(enrichment_df)
      )
    ) +
      
      labs(
        
        title = paste0(
          
          "KEGG enrichment: ",
          
          tp,
          
          ":00 – ",
          
          direction_label
        ),
        
        subtitle = "LP6 mature CAM leaves vs LP1 young C3 leaves"
      ) +
      
      theme_minimal(
        base_size = 12
      )
    
    
    ggsave(
      
      filename = file.path(
        
        "KEGG_corrected_results",
        
        paste0(
          
          "KEGG_",
          
          tp,
          
          "h_",
          
          direction_label,
          
          "_dotplot.png"
        )
      ),
      
      plot = pathway_plot,
      
      width = 10,
      
      height = 7,
      
      dpi = 300
    )
  }
  
  
  return(
    enrichment
  )
}


# ============================================================
# 18.2 RUN KEGG FOR ALL SIX TIMEPOINTS
# ============================================================

combined_kegg <- data.frame()


for (tp in time_order) {
  
  
  cat(
    
    "\n====================================\n",
    
    "KEGG enrichment for ",
    
    tp,
    
    ":00\n",
    
    "====================================\n",
    
    sep = ""
  )
  
  
  # ----------------------------------------------------------
  # Retrieve corrected LP6 vs LP1 result
  # ----------------------------------------------------------
  
  res_df <- all_de_results[[tp]]
  
  
  # ----------------------------------------------------------
  # Select significant DE genes
  #
  # padj < 0.05
  # |log2FC| > 1
  # ----------------------------------------------------------
  
  sig <- res_df %>%
    
    filter(
      
      !is.na(padj),
      
      padj < 0.05,
      
      abs(log2FoldChange) > 1
    )
  
  
  cat(
    "Total significant DE genes:",
    nrow(sig),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # ALL significant genes
  # ----------------------------------------------------------
  
  all_ids <- unique(
    sig$GeneID
  )
  
  
  # ----------------------------------------------------------
  # Higher in LP6
  #
  # Positive log2FC
  # Mature leaves
  # Full CAM
  # ----------------------------------------------------------
  
  lp6_ids <- sig %>%
    
    filter(
      log2FoldChange > 0
    ) %>%
    
    pull(
      GeneID
    ) %>%
    
    unique()
  
  
  # ----------------------------------------------------------
  # Higher in LP1
  #
  # Negative log2FC
  # Young leaves
  # Predominantly C3
  # ----------------------------------------------------------
  
  lp1_ids <- sig %>%
    
    filter(
      log2FoldChange < 0
    ) %>%
    
    pull(
      GeneID
    ) %>%
    
    unique()
  
  
  cat(
    "Higher in LP6 / mature CAM:",
    length(lp6_ids),
    "\n"
  )
  
  
  cat(
    "Higher in LP1 / young C3:",
    length(lp1_ids),
    "\n"
  )
  
  
  # ==========================================================
  # Run three enrichment analyses
  # ==========================================================
  
  
  result_all <- run_pathway_enrichment(
    
    all_ids,
    
    tp,
    
    "All_significant"
  )
  
  
  result_lp6 <- run_pathway_enrichment(
    
    lp6_ids,
    
    tp,
    
    "Higher_in_LP6_CAM"
  )
  
  
  result_lp1 <- run_pathway_enrichment(
    
    lp1_ids,
    
    tp,
    
    "Higher_in_LP1_C3"
  )
  
  
  # ==========================================================
  # Store results
  # ==========================================================
  
  result_list <- list(
    
    All_significant = result_all,
    
    Higher_in_LP6_CAM = result_lp6,
    
    Higher_in_LP1_C3 = result_lp1
  )
  
  
  for (
    direction_name in names(result_list)
  ) {
    
    
    current_result <- result_list[[direction_name]]
    
    
    if (!is.null(current_result)) {
      
      
      current_df <- as.data.frame(
        current_result
      )
      
      
      if (nrow(current_df) > 0) {
        
        
        current_df$Timepoint <- tp
        
        current_df$Direction <- direction_name
        
        
        combined_kegg <- bind_rows(
          
          combined_kegg,
          
          current_df
        )
      }
    }
  }
}


# ============================================================
# 19. SAVE COMBINED KEGG RESULTS
# ============================================================

if (nrow(combined_kegg) > 0) {
  
  
  combined_kegg <- combined_kegg %>%
    
    arrange(
      
      factor(
        Timepoint,
        levels = time_order
      ),
      
      Direction,
      
      p.adjust
    )
  
  
  # ----------------------------------------------------------
  # Save ALL significant enrichment results
  # ----------------------------------------------------------
  
  write.csv(
    
    combined_kegg,
    
    "KEGG_corrected_results/KEGG_all_corrected_enrichment_results.csv",
    
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 19.1 TOP 10 PATHWAYS PER TIMEPOINT AND DIRECTION
  # ==========================================================
  
  top10_kegg <- combined_kegg %>%
    
    group_by(
      
      Timepoint,
      
      Direction
    ) %>%
    
    slice_min(
      
      order_by = p.adjust,
      
      n = 10,
      
      with_ties = FALSE
    ) %>%
    
    ungroup()
  
  
  write.csv(
    
    top10_kegg,
    
    "KEGG_corrected_results/KEGG_top10_corrected_pathways.csv",
    
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 19.2 PRINT SUMMARY
  # ==========================================================
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    "TOP KEGG ENRICHMENT RESULTS\n"
  )
  
  cat(
    "====================================\n\n"
  )
  
  
  print(
    
    top10_kegg %>%
      
      select(
        
        Timepoint,
        
        Direction,
        
        ID,
        
        GeneRatio,
        
        BgRatio,
        
        pvalue,
        
        p.adjust,
        
        qvalue,
        
        Count
      )
  )
  
  
  # ==========================================================
  # 19.3 NUMBER OF ENRICHED PATHWAYS
  # ==========================================================
  
  pathway_summary <- combined_kegg %>%
    
    group_by(
      
      Timepoint,
      
      Direction
    ) %>%
    
    summarise(
      
      Enriched_pathways = n(),
      
      .groups = "drop"
    )
  
  
  cat(
    "\nNumber of enriched pathways:\n"
  )
  
  print(
    pathway_summary
  )
  
  
  write.csv(
    
    pathway_summary,
    
    "KEGG_corrected_results/KEGG_enriched_pathway_summary.csv",
    
    row.names = FALSE
  )
  
  
} else {
  
  
  message(
    
    paste(
      
      "No KEGG pathways passed the selected",
      
      "enrichment thresholds."
    )
  )
}


# ============================================================
# 19.4 FINAL KEGG ANALYSIS SUMMARY
# ============================================================

cat(
  "\n====================================\n"
)

cat(
  "KEGG ANALYSIS COMPLETE\n"
)

cat(
  "====================================\n"
)


cat(
  
  "\nBiological interpretation of direction:\n",
  
  "\nPositive log2FC:",
  
  "\n  Higher expression in LP6",
  
  "\n  Mature leaf",
  
  "\n  Fully developed CAM\n",
  
  "\nNegative log2FC:",
  
  "\n  Higher expression in LP1",
  
  "\n  Young leaf",
  
  "\n  Predominantly C3 photosynthesis\n",
  
  sep = ""
)


# ============================================================
# 20. SAVE NORMALISED EXPRESSION MATRIX
# ============================================================

normalised_counts <- counts(
  dds,
  normalized = TRUE
)

write.csv(
  normalised_counts,
  "DESeq2_corrected_results/DESeq2_normalised_counts.csv"
)


# ============================================================
# 21. SAVE SESSION INFORMATION
# ============================================================

writeLines(
  capture.output(sessionInfo()),
  "DESeq2_corrected_results/R_sessionInfo.txt"
)


# ============================================================
# ============================================================

cat(
  "\nAnalysis completed.\n",
  "Primary interpretation:\n",
  "  log2FC > 0 = higher in LP6 / mature CAM leaf\n",
  "  log2FC < 0 = higher in LP1 / young C3 leaf\n",
  "\nDo not use the previous mesophyll/epidermal outputs as final results.\n",
  sep = ""
)

# ============================================================
# 20. ADD KEGG PATHWAY NAMES - ROBUST VERSION
# ============================================================

library(KEGGREST)
library(dplyr)


# ============================================================
# 20.1 DOWNLOAD KEGG PATHWAY NAMES
# ============================================================

kegg_pathway_list <- KEGGREST::keggList(
  "pathway",
  "ko"
)


# ------------------------------------------------------------
# First inspect exactly what KEGG returned
# ------------------------------------------------------------

cat("\nFirst KEGG pathway identifiers returned:\n")

print(
  head(
    names(kegg_pathway_list),
    20
  )
)


cat("\nFirst KEGG pathway names returned:\n")

print(
  head(
    as.character(kegg_pathway_list),
    20
  )
)


# ============================================================
# 20.2 BUILD ROBUST PATHWAY LOOKUP TABLE
# ============================================================

pathway_lookup <- data.frame(
  
  KEGG_full_ID = names(kegg_pathway_list),
  
  Pathway_Name = as.character(kegg_pathway_list),
  
  stringsAsFactors = FALSE
)


# IMPORTANT:
# Instead of assuming the prefix is "path:ko",
# extract the final 5 digits from the KEGG identifier.
#
# Examples:
#
# path:ko00710  -> 00710
# path:map00710 -> 00710
# ko00710       -> 00710
# 00710         -> 00710

pathway_lookup <- pathway_lookup %>%
  
  mutate(
    
    ID = sub(
      ".*?([0-9]{5})$",
      "\\1",
      KEGG_full_ID
    )
    
  ) %>%
  
  filter(
    grepl(
      "^[0-9]{5}$",
      ID
    )
  ) %>%
  
  select(
    ID,
    Pathway_Name
  ) %>%
  
  distinct()


# ============================================================
# 20.3 CHECK THE LOOKUP TABLE
# ============================================================

cat(
  "\nNumber of pathways in lookup table:",
  nrow(pathway_lookup),
  "\n"
)


cat(
  "\nFirst 20 cleaned KEGG pathways:\n"
)

print(
  head(
    pathway_lookup,
    20
  )
)


# ============================================================
# 20.4 CLEAN IDs FROM YOUR ENRICHMENT RESULTS
# ============================================================

# Do the same thing to your enrichment IDs.
#
# This avoids problems with:
# 00710
# ko00710
# map00710
# path:ko00710
# etc.

combined_kegg_fixed <- combined_kegg %>%
  
  mutate(
    
    Original_ID = as.character(ID),
    
    ID = sub(
      ".*?([0-9]{5})$",
      "\\1",
      as.character(ID)
    )
    
  )


# ------------------------------------------------------------
# Check your enrichment IDs
# ------------------------------------------------------------

cat(
  "\nFirst pathway IDs in your enrichment results:\n"
)

print(
  head(
    unique(combined_kegg_fixed$ID),
    30
  )
)


# ============================================================
# 20.5 CHECK WHETHER THE TWO TABLES ACTUALLY OVERLAP
# ============================================================

shared_pathway_ids <- base::intersect(
  
  unique(
    combined_kegg_fixed$ID
  ),
  
  unique(
    pathway_lookup$ID
  )
)


cat(
  "\nNumber of shared pathway IDs:",
  length(shared_pathway_ids),
  "\n"
)


cat(
  "\nExample shared pathway IDs:\n"
)

print(
  head(
    shared_pathway_ids,
    30
  )
)


# ------------------------------------------------------------
# STOP if there is still no overlap
# ------------------------------------------------------------

if (length(shared_pathway_ids) == 0) {
  
  cat(
    "\nEnrichment IDs:\n"
  )
  
  print(
    head(
      unique(combined_kegg_fixed$ID),
      30
    )
  )
  
  
  cat(
    "\nKEGG lookup IDs:\n"
  )
  
  print(
    head(
      unique(pathway_lookup$ID),
      30
    )
  )
  
  
  stop(
    paste(
      "There is still no overlap between the enrichment",
      "pathway IDs and the KEGG pathway-name lookup.",
      "Do not continue to Sections 22 or 24."
    )
  )
}


# ============================================================
# 20.6 JOIN PATHWAY NAMES
# ============================================================

combined_kegg_named <- combined_kegg_fixed %>%
  
  left_join(
    pathway_lookup,
    by = "ID"
  )


# ============================================================
# 20.7 CHECK SUCCESS
# ============================================================

cat(
  "\n====================================\n"
)

cat(
  "PATHWAY NAME MAPPING CHECK\n"
)

cat(
  "====================================\n"
)


cat(
  "\nPathway names successfully assigned:",
  sum(
    !is.na(
      combined_kegg_named$Pathway_Name
    )
  ),
  "of",
  nrow(combined_kegg_named),
  "\n"
)


cat(
  "Missing pathway names:",
  sum(
    is.na(
      combined_kegg_named$Pathway_Name
    )
  ),
  "\n"
)


cat(
  "Percentage successfully mapped:",
  round(
    100 *
      mean(
        !is.na(
          combined_kegg_named$Pathway_Name
        )
      ),
    2
  ),
  "%\n"
)


# ============================================================
# 20.8 PRINT ACTUAL PATHWAY NAMES
# ============================================================

cat(
  "\nExample enrichment results with pathway names:\n"
)


print(
  
  combined_kegg_named %>%
    
    select(
      ID,
      Pathway_Name,
      Timepoint,
      Direction,
      p.adjust
    ) %>%
    
    distinct() %>%
    
    head(
      40
    )
)


# ============================================================
# 20.9 CHECK IMPORTANT CAM-RELATED PATHWAYS
# ============================================================

key_pathway_ids <- c(
  
  "00710",  # Carbon fixation
  "01200",  # Carbon metabolism
  "00500",  # Starch and sucrose metabolism
  "00010",  # Glycolysis / Gluconeogenesis
  "04712",  # Circadian rhythm - plant
  "04075",  # Plant hormone signal transduction
  "00190"   # Oxidative phosphorylation
)


cat(
  "\n====================================\n"
)

cat(
  "KEY CAM-RELEVANT PATHWAYS PRESENT\n"
)

cat(
  "====================================\n\n"
)


key_pathway_check <- combined_kegg_named %>%
  
  filter(
    ID %in% key_pathway_ids
  ) %>%
  
  select(
    ID,
    Pathway_Name
  ) %>%
  
  distinct()


print(
  key_pathway_check
)


# ============================================================
# 20.10 EXTRACT BIOLOGICALLY RELEVANT PATHWAYS
# ============================================================

important_pathways <- combined_kegg_named %>%
  
  filter(
    
    !is.na(Pathway_Name),
    
    grepl(
      
      paste(
        
        c(
          "Carbon fixation",
          "Carbon metabolism",
          "Starch and sucrose",
          "Glycolysis",
          "Circadian",
          "Photosynthesis",
          "Plant hormone",
          "Pyruvate",
          "Oxidative phosphorylation"
        ),
        
        collapse = "|"
      ),
      
      Pathway_Name,
      
      ignore.case = TRUE
    )
  ) %>%
  
  arrange(
    
    factor(
      Timepoint,
      levels = time_order
    ),
    
    Direction,
    
    p.adjust
  )


# ============================================================
# 20.11 SAVE RESULTS
# ============================================================

write.csv(
  
  combined_kegg_named,
  
  "KEGG_corrected_results/KEGG_all_results_with_pathway_names.csv",
  
  row.names = FALSE
)


write.csv(
  
  important_pathways,
  
  "KEGG_corrected_results/CAM_relevant_enriched_pathways.csv",
  
  row.names = FALSE
)


# ============================================================
# 20.12 PRINT IMPORTANT PATHWAYS
# ============================================================

cat(
  "\n====================================\n"
)

cat(
  "CAM-RELEVANT ENRICHED PATHWAYS\n"
)

cat(
  "====================================\n\n"
)


print(
  important_pathways %>%
    select(
      Timepoint,
      Direction,
      ID,
      Pathway_Name,
      GeneRatio,
      BgRatio,
      p.adjust,
      Count
    )
)


# ============================================================
# 20.13 FINAL SUMMARY
# ============================================================

cat(
  "\n====================================\n"
)

cat(
  "KEGG PATHWAY NAME MAPPING COMPLETE\n"
)

cat(
  "====================================\n"
)


cat(
  "\nPathway names successfully assigned to:",
  sum(
    !is.na(
      combined_kegg_named$Pathway_Name
    )
  ),
  "of",
  nrow(combined_kegg_named),
  "enrichment results.\n"
)
# ============================================================
# 21. TARGETED CAM-GENE TIME-COURSE ANALYSIS
# ============================================================

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)


# ------------------------------------------------------------
# 21.1 Obtain variance-stabilised expression
# ------------------------------------------------------------

# vsd <- vst(dds, blind = FALSE)

vst_matrix <- assay(vsd)

vst_df <- as.data.frame(vst_matrix)

vst_df$GeneID <- rownames(vst_df)


# ------------------------------------------------------------
# 21.2 Identify CAM-related genes from annotation
# ------------------------------------------------------------

cam_annotation <- anno_clean %>%
  mutate(
    
    CAM_family = case_when(
      
      grepl(
        "phosphoenolpyruvate carboxylase kinase|PPCK",
        paste(Gene_Name, Description),
        ignore.case = TRUE
      ) ~ "PPCK",
      
      grepl(
        "phosphoenolpyruvate carboxylase|PEPC",
        paste(Gene_Name, Description),
        ignore.case = TRUE
      ) ~ "PEPC",
      
      grepl(
        "pyruvate.*phosphate dikinase|PPDK",
        paste(Gene_Name, Description),
        ignore.case = TRUE
      ) ~ "PPDK",
      
      grepl(
        "malate dehydrogenase",
        paste(Gene_Name, Description),
        ignore.case = TRUE
      ) ~ "Malate dehydrogenase",
      
      grepl(
        "malic enzyme",
        paste(Gene_Name, Description),
        ignore.case = TRUE
      ) ~ "Malic enzyme",
      
      grepl(
        "carbonic anhydrase",
        paste(Gene_Name, Description),
        ignore.case = TRUE
      ) ~ "Carbonic anhydrase",
      
      grepl(
        "starch|sucrose",
        paste(Gene_Name, Description),
        ignore.case = TRUE
      ) ~ "Carbohydrate metabolism",
      
      grepl(
        "vacuolar|tonoplast|V-type.*ATPase|V-ATPase",
        paste(Gene_Name, Description),
        ignore.case = TRUE
      ) ~ "Vacuolar transport",
      
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(CAM_family)
  )


write.csv(
  cam_annotation,
  "Annotated_corrected_results/All_CAM_related_annotations.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 21.3 Rank CAM genes using corrected DESeq2 results
# ------------------------------------------------------------

all_de_long <- bind_rows(
  lapply(
    names(all_de_results),
    function(tp) {
      
      x <- all_de_results[[tp]]
      
      x$Timepoint <- tp
      
      x
    }
  )
)


cam_de <- all_de_long %>%
  inner_join(
    cam_annotation,
    by = "GeneID"
  )


# Rank genes according to their strongest statistical evidence
# across the diel cycle.

cam_rank <- cam_de %>%
  group_by(
    GeneID,
    Gene_Name,
    CAM_family
  ) %>%
  summarise(
    
    Minimum_padj = min(
      padj,
      na.rm = TRUE
    ),
    
    Maximum_abs_log2FC = max(
      abs(log2FoldChange),
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  arrange(
    Minimum_padj
  )


write.csv(
  cam_rank,
  "Annotated_corrected_results/CAM_gene_ranking.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 21.4 Select representative CAM genes
# ------------------------------------------------------------

# Select up to 2 strongest genes from each functional family.

selected_cam <- cam_rank %>%
  group_by(
    CAM_family
  ) %>%
  slice_min(
    order_by = Minimum_padj,
    n = 2,
    with_ties = FALSE
  ) %>%
  ungroup()


selected_cam_genes <- unique(
  selected_cam$GeneID
)


print(
  selected_cam
)


# ------------------------------------------------------------
# 21.5 Convert VST expression into long format
# ------------------------------------------------------------

cam_expression <- vst_df %>%
  filter(
    GeneID %in% selected_cam_genes
  ) %>%
  
  pivot_longer(
    cols = -GeneID,
    names_to = "Sample",
    values_to = "VST_expression"
  ) %>%
  
  left_join(
    data.frame(
      Sample = rownames(sample_info),
      sample_info,
      row.names = NULL
    ),
    by = "Sample"
  ) %>%
  
  left_join(
    selected_cam %>%
      select(
        GeneID,
        Gene_Name,
        CAM_family
      ),
    by = "GeneID"
  )


# ------------------------------------------------------------
# 21.6 Average biological replicates
# ------------------------------------------------------------

cam_expression_summary <- cam_expression %>%
  
  group_by(
    GeneID,
    Gene_Name,
    CAM_family,
    leaf,
    timepoint
  ) %>%
  
  summarise(
    
    Mean_expression = mean(
      VST_expression
    ),
    
    SD_expression = sd(
      VST_expression
    ),
    
    .groups = "drop"
  )


cam_expression_summary$timepoint <- factor(
  cam_expression_summary$timepoint,
  levels = time_order
)


# ------------------------------------------------------------
# 21.7 Create time-course figure
# ------------------------------------------------------------

cam_time_plot <- ggplot(
  cam_expression_summary,
  aes(
    x = timepoint,
    y = Mean_expression,
    group = leaf,
    colour = leaf
  )
) +
  
  geom_line(
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 2
  ) +
  
  geom_errorbar(
    aes(
      ymin = Mean_expression - SD_expression,
      ymax = Mean_expression + SD_expression
    ),
    width = 0.15
  ) +
  
  facet_wrap(
    ~ CAM_family + Gene_Name,
    scales = "free_y"
  ) +
  
  labs(
    title = "Diel expression of CAM-associated genes",
    subtitle = "LP1 young C3 leaves vs LP6 mature CAM leaves",
    x = "Sampling time",
    y = "Variance-stabilised expression",
    colour = "Leaf pair"
  ) +
  
  theme_minimal(
    base_size = 11
  )


print(
  cam_time_plot
)


ggsave(
  "Annotated_corrected_results/CAM_gene_timecourse.png",
  cam_time_plot,
  width = 13,
  height = 10,
  dpi = 300
)


# ============================================================
# 22. KEGG PATHWAY x TIME HEATMAP
# ============================================================
# ------------------------------------------------------------
# 22.1 Prepare enrichment strength
# ------------------------------------------------------------

heatmap_data <- combined_kegg_named %>%
  
  filter(
    !is.na(p.adjust),
    p.adjust < 0.05
  ) %>%
  
  mutate(
    
    Enrichment_strength = -log10(
      pmax(
        p.adjust,
        .Machine$double.xmin
      )
    ),
    
    Timepoint = factor(
      Timepoint,
      levels = time_order
    )
  )


# ------------------------------------------------------------
# 22.2 Select biologically relevant pathways
# ------------------------------------------------------------

cam_pathway_data <- heatmap_data %>%
  
  filter(
    
    grepl(
      
      paste(
        
        c(
          "Carbon fixation",
          "Carbon metabolism",
          "Starch and sucrose",
          "Glycolysis",
          "Pyruvate",
          "Circadian",
          "Photosynthesis",
          "Plant hormone",
          "Oxidative phosphorylation"
        ),
        
        collapse = "|"
      ),
      
      Pathway_Name,
      
      ignore.case = TRUE
    )
  )


# ------------------------------------------------------------
# 22.3 LP6 / mature CAM heatmap
# ------------------------------------------------------------

lp6_heatmap_data <- cam_pathway_data %>%
  
  filter(
    Direction == "Higher_in_LP6_CAM"
  )


lp6_heatmap <- ggplot(
  lp6_heatmap_data,
  aes(
    x = Timepoint,
    y = Pathway_Name,
    fill = Enrichment_strength
  )
) +
  
  geom_tile(
    colour = "white"
  ) +
  
  labs(
    title = "Pathway enrichment across the diel cycle",
    subtitle = "Genes expressed more highly in mature CAM LP6 leaves",
    x = "Sampling time",
    y = "KEGG pathway",
    fill = expression(-log[10](adjusted~p))
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    axis.text.y = element_text(
      size = 9
    )
  )


print(
  lp6_heatmap
)


ggsave(
  "KEGG_corrected_results/KEGG_LP6_CAM_time_heatmap.png",
  lp6_heatmap,
  width = 10,
  height = 8,
  dpi = 300
)


# ------------------------------------------------------------
# 22.4 LP1 / young C3 heatmap
# ------------------------------------------------------------

lp1_heatmap_data <- cam_pathway_data %>%
  
  filter(
    Direction == "Higher_in_LP1_C3"
  )


lp1_heatmap <- ggplot(
  lp1_heatmap_data,
  aes(
    x = Timepoint,
    y = Pathway_Name,
    fill = Enrichment_strength
  )
) +
  
  geom_tile(
    colour = "white"
  ) +
  
  labs(
    title = "Pathway enrichment across the diel cycle",
    subtitle = "Genes expressed more highly in young C3 LP1 leaves",
    x = "Sampling time",
    y = "KEGG pathway",
    fill = expression(-log[10](adjusted~p))
  ) +
  
  theme_minimal(
    base_size = 11
  ) +
  
  theme(
    axis.text.y = element_text(
      size = 9
    )
  )


print(
  lp1_heatmap
)


ggsave(
  "KEGG_corrected_results/KEGG_LP1_C3_time_heatmap.png",
  lp1_heatmap,
  width = 10,
  height = 8,
  dpi = 300
)


# ============================================================
# 23. LEAF-PAIR x TIMEPOINT INTERACTION ANALYSIS
# ============================================================

# ------------------------------------------------------------
# 23.1 Create interaction DESeq2 dataset
# ------------------------------------------------------------

dds_interaction <- DESeqDataSetFromMatrix(
  
  countData = counts_raw,
  
  colData = sample_info,
  
  design = ~ leaf + timepoint + leaf:timepoint
)


# Apply same filtering as main analysis

dds_interaction <- dds_interaction[
  
  rowSums(
    counts(dds_interaction)
  ) > 10,
  
]


# ------------------------------------------------------------
# 23.2 Likelihood-ratio test
# ------------------------------------------------------------

# Full model:
# leaf + time + leaf:time
#
# Reduced model:
# leaf + time
#
# The LRT tests whether the interaction contributes
# significantly to explaining expression.

dds_interaction <- DESeq(
  
  dds_interaction,
  
  test = "LRT",
  
  reduced = ~ leaf + timepoint
)


interaction_res <- results(
  dds_interaction,
  alpha = 0.05
)


interaction_df <- as.data.frame(
  interaction_res
)


interaction_df$GeneID <- rownames(
  interaction_df
)


# ------------------------------------------------------------
# 23.3 Significant time-dependent leaf differences
# ------------------------------------------------------------

interaction_sig <- interaction_df %>%
  
  filter(
    
    !is.na(padj),
    
    padj < 0.05
  ) %>%
  
  arrange(
    padj
  )


cat(
  "\nGenes with significant LP1/LP6 x time interaction:",
  nrow(interaction_sig),
  "\n"
)


# IMPORTANT:
#
# For the LRT, padj used to identify interaction genes.
# NOT interpreting the displayed log2FoldChange as an
# overall interaction effect size.


# ------------------------------------------------------------
# 23.4 Add biological annotation
# ------------------------------------------------------------

interaction_annotated <- interaction_sig %>%
  
  left_join(
    
    anno_clean,
    
    by = "GeneID"
  )


write.csv(
  
  interaction_annotated,
  
  "DESeq2_corrected_results/Significant_leaf_time_interaction_annotated.csv",
  
  row.names = FALSE
)


# ------------------------------------------------------------
# 23.5 Find CAM-related interaction genes
# ------------------------------------------------------------

cam_interaction <- interaction_annotated %>%
  
  filter(
    
    grepl(
      
      paste(
        
        c(
          "phosphoenolpyruvate",
          "PEPC",
          "PPCK",
          "pyruvate phosphate dikinase",
          "PPDK",
          "malate",
          "malic enzyme",
          "carbonic anhydrase",
          "circadian",
          "starch",
          "sucrose",
          "vacuolar",
          "tonoplast"
        ),
        
        collapse = "|"
      ),
      
      paste(
        Gene_Name,
        Description
      ),
      
      ignore.case = TRUE
    )
  )


write.csv(
  
  cam_interaction,
  
  "DESeq2_corrected_results/CAM_related_leaf_time_interaction_genes.csv",
  
  row.names = FALSE
)


cat(
  "CAM-related genes with significant leaf x time interaction:",
  nrow(cam_interaction),
  "\n"
)



# ============================================================
# 24. IDENTIFY GENES DRIVING IMPORTANT KEGG PATHWAYS
# ============================================================

# ------------------------------------------------------------
# 24.1 Select key CAM-relevant KEGG pathways
# ------------------------------------------------------------

key_pathway_ids <- c(
  
  "00710",   # Carbon fixation
  "01200",   # Carbon metabolism
  "00500",   # Starch and sucrose metabolism
  "00010",   # Glycolysis / Gluconeogenesis
  "04712",   # Circadian rhythm - plant
  "04075",   # Plant hormone signal transduction
  "00190"    # Oxidative phosphorylation
)


key_kegg <- combined_kegg_named %>%
  
  filter(
    ID %in% key_pathway_ids
  )


# ------------------------------------------------------------
# 24.2 Extracting gene IDs contained in enriched pathways
# ------------------------------------------------------------

# clusterProfiler stores contributing genes in geneID,
# separated by "/".

driver_genes <- key_kegg %>%
  
  select(
    
    Timepoint,
    
    Direction,
    
    ID,
    
    Pathway_Name,
    
    p.adjust,
    
    geneID
    
  ) %>%
  
  tidyr::separate_rows(
    
    geneID,
    
    sep = "/"
  ) %>%
  
  rename(
    GeneID = geneID
  )


# ------------------------------------------------------------
# 24.3 Add gene annotation
# ------------------------------------------------------------

driver_genes <- driver_genes %>%
  
  left_join(
    
    anno_clean,
    
    by = "GeneID"
  )


# ------------------------------------------------------------
# 24.4 Add corrected DESeq2 fold changes
# ------------------------------------------------------------

de_for_join <- bind_rows(
  
  lapply(
    
    names(all_de_results),
    
    function(tp) {
      
      x <- all_de_results[[tp]]
      
      x$Timepoint <- tp
      
      x
    }
  )
  
) %>%
  
  select(
    
    GeneID,
    
    Timepoint,
    
    log2FoldChange,
    
    padj
  )


driver_genes <- driver_genes %>%
  
  left_join(
    
    de_for_join,
    
    by = c(
      "GeneID",
      "Timepoint"
    )
  )


# ------------------------------------------------------------
# 24.5 Save complete driver-gene table
# ------------------------------------------------------------

write.csv(
  
  driver_genes,
  
  "KEGG_corrected_results/Genes_driving_CAM_relevant_KEGG_pathways.csv",
  
  row.names = FALSE
)


# ------------------------------------------------------------
# 24.6 Find repeatedly occurring genes
# ------------------------------------------------------------

recurrent_driver_genes <- driver_genes %>%
  
  group_by(
    
    GeneID,
    
    Gene_Name,
    
    Description,
    
    Pathway_Name
    
  ) %>%
  
  summarise(
    
    Number_of_timepoints = n_distinct(
      Timepoint
    ),
    
    Timepoints = paste(
      
      unique(Timepoint),
      
      collapse = ", "
      
    ),
    
    .groups = "drop"
  ) %>%
  
  arrange(
    desc(Number_of_timepoints)
  )


write.csv(
  
  recurrent_driver_genes,
  
  "KEGG_corrected_results/Recurrent_KEGG_driver_genes.csv",
  
  row.names = FALSE
)


head(
  recurrent_driver_genes,
  30
)

recurrent_driver_genes <- driver_genes %>%
  group_by(
    GeneID,
    Gene_Name,
    Description,
    ID,
    Pathway_Name,
    Direction
  ) %>%
  summarise(
    Number_of_timepoints = n_distinct(Timepoint),
    Timepoints = paste(
      sort(unique(Timepoint)),
      collapse = ", "
    ),
    Mean_log2FC = mean(
      log2FoldChange,
      na.rm = TRUE
    ),
    Maximum_abs_log2FC = max(
      abs(log2FoldChange),
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(
    desc(Number_of_timepoints),
    desc(Maximum_abs_log2FC)
  )

write.csv(
  recurrent_driver_genes,
  "KEGG_corrected_results/Recurrent_KEGG_driver_genes_final.csv",
  row.names = FALSE
)

print(
  recurrent_driver_genes %>%
    select(
      GeneID,
      Gene_Name,
      Description,
      Pathway_Name,
      Direction,
      Number_of_timepoints,
      Timepoints,
      Mean_log2FC,
      Maximum_abs_log2FC
    ),
  n = 50,
  width = Inf
)


# ============================================================
# 25. IDENTIFY KEY CAM CANDIDATE GENES
# ============================================================

cam_candidates <- recurrent_driver_genes %>%
  
  filter(
    Direction != "All_significant"
  ) %>%
  
  filter(
    grepl(
      paste(
        c(
          "phosphoenolpyruvate carboxylase",
          "PEP carboxylase",
          "PEPC",
          "PPCK",
          "phosphoenolpyruvate carboxylase kinase",
          "malate dehydrogenase",
          "malic enzyme",
          "pyruvate phosphate dikinase",
          "PPDK",
          "carbonic anhydrase",
          "starch",
          "sucrose",
          "phosphorylase",
          "V-type",
          "V-ATPase",
          "proton",
          "tonoplast"
        ),
        collapse = "|"
      ),
      paste(
        Gene_Name,
        Description
      ),
      ignore.case = TRUE
    )
  ) %>%
  
  arrange(
    desc(Number_of_timepoints),
    desc(Maximum_abs_log2FC)
  )


print(
  cam_candidates %>%
    select(
      GeneID,
      Gene_Name,
      Description,
      Pathway_Name,
      Direction,
      Number_of_timepoints,
      Timepoints,
      Mean_log2FC,
      Maximum_abs_log2FC
    ),
  n = 100,
  width = Inf
)


write.csv(
  cam_candidates,
  "KEGG_corrected_results/Key_CAM_candidate_genes.csv",
  row.names = FALSE
)


# ============================================================
# 26. SEARCH SPECIFICALLY FOR CORE CAM GENES
# ============================================================

core_cam_search <- anno_clean %>%
  filter(
    grepl(
      paste(
        c(
          "phosphoenolpyruvate carboxylase",
          "phosphoenolpyruvate carboxylase kinase",
          "PEPC",
          "PPCK",
          "pyruvate phosphate dikinase",
          "malate dehydrogenase",
          "malic enzyme",
          "carbonic anhydrase"
        ),
        collapse = "|"
      ),
      paste(
        Gene_Name,
        Description
      ),
      ignore.case = TRUE
    )
  ) %>%
  distinct(
    GeneID,
    Gene_Name,
    Description,
    .keep_all = TRUE
  )

core_cam_search %>%
  select(
    GeneID,
    Gene_Name,
    Description
  ) %>%
  head(100)

write.csv(
  core_cam_search,
  "KEGG_corrected_results/Core_CAM_gene_search.csv",
  row.names = FALSE
)


# ============================================================
# 27. CORE CAM GENE EXPRESSION ACROSS THE DIEL CYCLE
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)


# ------------------------------------------------------------
# 27.1 Combine LP6 vs LP1 DESeq2 results from all timepoints
# ------------------------------------------------------------

cam_timecourse <- bind_rows(
  lapply(
    names(all_de_results),
    function(tp) {
      
      x <- all_de_results[[tp]]
      
      x$Timepoint <- tp
      
      x
    }
  )
)


# ------------------------------------------------------------
# 27.2 Keep core CAM-related genes
# ------------------------------------------------------------

cam_timecourse <- cam_timecourse %>%
  
  filter(
    GeneID %in% core_cam_search$GeneID
  ) %>%
  
  left_join(
    core_cam_search %>%
      select(
        GeneID,
        Gene_Name,
        Description
      ) %>%
      distinct(GeneID, .keep_all = TRUE),
    
    by = "GeneID"
  )


# ------------------------------------------------------------
# 27.3 Set correct diel order
# ------------------------------------------------------------

cam_timecourse$Timepoint <- factor(
  cam_timecourse$Timepoint,
  levels = c(
    "10",
    "14",
    "18",
    "22",
    "02",
    "06"
  )
)


# ------------------------------------------------------------
# 27.4 Classify major CAM functions
# ------------------------------------------------------------

cam_timecourse <- cam_timecourse %>%
  
  mutate(
    
    CAM_Function = case_when(
      
      grepl(
        "carboxylase kinase",
        Description,
        ignore.case = TRUE
      ) ~ "PPCK",
      
      grepl(
        "phosphoenolpyruvate carboxylase",
        Description,
        ignore.case = TRUE
      ) ~ "PEPC",
      
      grepl(
        "pyruvate phosphate dikinase",
        Description,
        ignore.case = TRUE
      ) ~ "PPDK",
      
      grepl(
        "malate dehydrogenase",
        Description,
        ignore.case = TRUE
      ) ~ "MDH",
      
      grepl(
        "malic enzyme",
        Description,
        ignore.case = TRUE
      ) ~ "Malic enzyme",
      
      grepl(
        "carbonic anhydrase",
        Description,
        ignore.case = TRUE
      ) ~ "Carbonic anhydrase",
      
      TRUE ~ "Other"
    )
  )

table(
  cam_timecourse$CAM_Function
)

cam_timecourse %>%
  group_by(CAM_Function) %>%
  summarise(
    Unique_genes = n_distinct(GeneID),
    Rows = n(),
    .groups = "drop"
  )

cat(
  "Unique core CAM genes in time-course:",
  n_distinct(cam_timecourse$GeneID),
  "\n"
)


# ============================================================
# 27.5 SUMMARISE CORE CAM CANDIDATES
# ============================================================

cam_candidate_summary <- cam_timecourse %>%
  
  group_by(
    GeneID,
    Gene_Name,
    Description,
    CAM_Function
  ) %>%
  
  summarise(
    
    Significant_timepoints = sum(
      !is.na(padj) &
        padj < 0.05 &
        abs(log2FoldChange) > 1
    ),
    
    LP6_higher_timepoints = sum(
      !is.na(padj) &
        padj < 0.05 &
        log2FoldChange > 1
    ),
    
    LP1_higher_timepoints = sum(
      !is.na(padj) &
        padj < 0.05 &
        log2FoldChange < -1
    ),
    
    Mean_log2FC = mean(
      log2FoldChange,
      na.rm = TRUE
    ),
    
    Maximum_abs_log2FC = max(
      abs(log2FoldChange),
      na.rm = TRUE
    ),
    
    .groups = "drop"
  ) %>%
  
  arrange(
    CAM_Function,
    desc(Significant_timepoints),
    desc(Maximum_abs_log2FC)
  )


# Save full table
write.csv(
  cam_candidate_summary,
  "KEGG_corrected_results/Core_CAM_candidate_summary.csv",
  row.names = FALSE
)

cam_candidate_summary %>%
  filter(
    CAM_Function %in% c("PEPC", "PPCK", "PPDK")
  ) %>%
  select(
    GeneID,
    Gene_Name,
    CAM_Function,
    Significant_timepoints,
    LP6_higher_timepoints,
    LP1_higher_timepoints,
    Mean_log2FC,
    Maximum_abs_log2FC
  ) %>%
  as.data.frame()


# ============================================================
# CREATE CORE PEPC / PPCK / PPDK TIMECOURSE TABLE
# ============================================================

core_cam_11 <- cam_timecourse %>%
  filter(
    CAM_Function %in% c("PEPC", "PPCK", "PPDK")
  ) %>%
  select(
    GeneID,
    Gene_Name,
    CAM_Function,
    Timepoint,
    log2FoldChange,
    padj
  ) %>%
  arrange(
    CAM_Function,
    GeneID,
    factor(Timepoint, levels = time_order)
  )

head(core_cam_11, 30)

core_cam_11_wide <- core_cam_11 %>%
  select(
    GeneID,
    Gene_Name,
    CAM_Function,
    Timepoint,
    log2FoldChange
  ) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = log2FoldChange,
    names_prefix = "log2FC_"
  )

as.data.frame(core_cam_11_wide)

core_cam_11_padj_wide <- core_cam_11 %>%
  select(
    GeneID,
    Gene_Name,
    CAM_Function,
    Timepoint,
    padj
  ) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = padj,
    names_prefix = "padj_"
  )

as.data.frame(core_cam_11_padj_wide)


# ============================================================
# 27.7 PEPC / PPCK / PPDK DIEL HEATMAP
# ============================================================


library(ggplot2)
library(dplyr)
library(scales)

core_cam_plot <- core_cam_11 %>%
  mutate(
    Timepoint = factor(
      Timepoint,
      levels = c("10", "14", "18", "22", "02", "06")
    ),
    
    Gene_Label = case_when(
      !is.na(Gene_Name) ~ paste0(
        Gene_Name, " (", sub("-RA$", "", GeneID), ")"
      ),
      
      TRUE ~ paste0(
        CAM_Function, " (", sub("-RA$", "", GeneID), ")"
      )
    )
  )


p_core_cam_final <- ggplot(
  core_cam_plot,
  aes(
    x = Timepoint,
    y = Gene_Label,
    fill = log2FoldChange
  )
) +
  
  # Heatmap cells
  geom_tile(
    colour = "white",
    linewidth = 0.7
  ) +
  
  # Significance stars
  geom_text(
    aes(
      label = ifelse(
        !is.na(padj) &
          padj < 0.05 &
          abs(log2FoldChange) > 1,
        "*",
        ""
      )
    ),
    size = 5,
    fontface = "bold"
  ) +
  
  # Separate CAM gene families
  facet_grid(
    CAM_Function ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  
  # Symmetrical colour scale
  scale_fill_gradient2(
    low = "#B2182B",
    mid = "#F7F7F7",
    high = "#2166AC",
    midpoint = 0,
    limits = c(-8, 8),
    oob = scales::squish,
    breaks = c(-8, -4, 0, 4, 8),
    name = expression(
      "LP6 vs LP1 " ~ log[2] ~ " fold change"
    )
  ) +
  
  labs(
    title = "Diel variation in core CAM-associated gene expression",
    x = "Sampling time (h)",
    y = NULL
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    # Main title
    plot.title = element_text(
      face = "bold",
      size = 17,
      hjust = 0
    ),
    
    # Axis labels
    axis.title.x = element_text(
      face = "bold",
      size = 13,
      margin = margin(t = 8)
    ),
    
    axis.text.x = element_text(
      size = 11,
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    # Remove axis ticks/lines around heatmap
    axis.ticks.y = element_blank(),
    
    # Facet labels
    strip.background = element_blank(),
    
    strip.text.y = element_text(
      angle = 0,
      face = "bold",
      size = 11
    ),
    
    # Clean panel
    panel.grid = element_blank(),
    
    # Space between PEPC / PPCK / PPDK
    panel.spacing.y = unit(0.25, "cm"),
    
    # Legend
    legend.title = element_text(
      size = 11
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    legend.key.height = unit(1.2, "cm"),
    
    # Overall margins
    plot.margin = margin(
      t = 15,
      r = 20,
      b = 10,
      l = 10
    )
  )


print(p_core_cam_final)

ggsave(
  filename = "KEGG_corrected_results/Figure_Core_CAM_Diel_Heatmap_FINAL.png",
  plot = p_core_cam_final,
  width = 10.5,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ------------------------------------------------------------
# 2. Vector PDF
# ------------------------------------------------------------

ggsave(
  filename = "KEGG_corrected_results/Figure_Core_CAM_Diel_Heatmap_FINAL.pdf",
  plot = p_core_cam_final,
  width = 10.5,
  height = 7.5,
  units = "in",
  device = "pdf",
  bg = "white"
)

cat(
  "\nFinal CAM heatmap saved successfully to:",
  "\nKEGG_corrected_results/Figure_Core_CAM_Diel_Heatmap_FINAL.png",
  "\nKEGG_corrected_results/Figure_Core_CAM_Diel_Heatmap_FINAL.pdf\n"
)


# ============================================================
# 28. LEAF × TIME INTERACTION ANALYSIS
# ============================================================

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)

# ------------------------------------------------------------
# 28.1 CHECK AND SET METADATA FACTOR LEVELS
# ------------------------------------------------------------

colData(dds)

# LP1 = reference level
# LP6 = comparison level

dds$leaf <- factor(
  dds$leaf,
  levels = c("LP1", "LP6")
)

# Chronological diel order
dds$timepoint <- factor(
  dds$timepoint,
  levels = c("10", "14", "18", "22", "02", "06")
)

# Check
levels(dds$leaf)
levels(dds$timepoint)

table(
  dds$leaf,
  dds$timepoint
)

# ============================================================
# 28.2 LEAF × TIME INTERACTION MODEL
# ============================================================

# Full model:
# leaf + time + leaf × time interaction

design(dds) <- ~ leaf + timepoint + leaf:timepoint

# Check design
design(dds)

~leaf + timepoint + leaf:timepoint

# ============================================================
# 28.3 LIKELIHOOD RATIO TEST
# ============================================================

dds_interaction <- DESeq(
  dds,
  test = "LRT",
  reduced = ~ leaf + timepoint
)

interaction_res <- results(
  dds_interaction,
  alpha = 0.05
)

interaction_df <- as.data.frame(
  interaction_res
)

interaction_df$GeneID <- rownames(
  interaction_df
)

interaction_df <- interaction_df %>%
  relocate(GeneID)

head(interaction_df)

interaction_df <- interaction_df %>%
  mutate(
    Interaction_significant = ifelse(
      !is.na(padj) & padj < 0.05,
      "Significant",
      "Not significant"
    )
  )

table(
  interaction_df$Interaction_significant
)

n_interaction_sig <- sum(
  !is.na(interaction_df$padj) &
    interaction_df$padj < 0.05
)

n_interaction_tested <- sum(
  !is.na(interaction_df$padj)
)

cat(
  "\n========================================\n",
  "LEAF × TIME INTERACTION RESULTS\n",
  "========================================\n",
  "Genes with valid adjusted p-values:",
  n_interaction_tested,
  "\nSignificant interaction genes:",
  n_interaction_sig,
  "\nPercentage significant:",
  round(
    100 * n_interaction_sig /
      n_interaction_tested,
    2
  ),
  "%\n",
  "========================================\n"
)

# ============================================================
# 28.6 EXTRACT SIGNIFICANT LEAF × TIME INTERACTION GENES
# ============================================================

significant_interaction_genes <- interaction_df %>%
  filter(
    !is.na(padj),
    padj < 0.05
  ) %>%
  arrange(padj)

cat(
  "\nNumber of significant interaction genes:",
  nrow(significant_interaction_genes),
  "\n"
)

head(
  significant_interaction_genes,
  20
)

write.csv(
  significant_interaction_genes,
  "KEGG_corrected_results/Leaf_Time_Interaction_SignificantGenes.csv",
  row.names = FALSE
)


# ============================================================
# 28.7 ANNOTATE SIGNIFICANT INTERACTION GENES
# ============================================================

significant_interaction_annotated <-
  significant_interaction_genes %>%
  left_join(
    anno_clean,
    by = "GeneID"
  )

head(
  significant_interaction_annotated,
  20
)

write.csv(
  significant_interaction_annotated,
  "KEGG_corrected_results/Leaf_Time_Interaction_SignificantGenes_Annotated.csv",
  row.names = FALSE
)


# ============================================================
# 29. CORE CAM GENES × INTERACTION TEST
# ============================================================

core_cam_unique <- cam_timecourse %>%
  select(
    GeneID,
    Gene_Name,
    Description,
    CAM_Function
  ) %>%
  distinct()


cam_interaction_results <- core_cam_unique %>%
  left_join(
    interaction_df %>%
      select(
        GeneID,
        stat,
        pvalue,
        padj,
        Interaction_significant
      ),
    by = "GeneID"
  ) %>%
  arrange(
    CAM_Function,
    padj
  )


# Check
head(
  cam_interaction_results,
  30
)
cat(
  "\nCore CAM genes:",
  nrow(core_cam_unique),
  "\nCore CAM genes matched to interaction test:",
  sum(!is.na(cam_interaction_results$pvalue)),
  "\n"
)


cam_interaction_summary <- cam_interaction_results %>%
  group_by(CAM_Function) %>%
  summarise(
    Total_genes = n(),
    
    Significant_interaction = sum(
      !is.na(padj) &
        padj < 0.05
    ),
    
    Percent_significant = round(
      100 *
        Significant_interaction /
        Total_genes,
      1
    ),
    
    .groups = "drop"
  )

as.data.frame(cam_interaction_summary)

core_cam_interaction_11 <- cam_interaction_results %>%
  filter(
    CAM_Function %in% c(
      "PEPC",
      "PPCK",
      "PPDK"
    )
  ) %>%
  select(
    GeneID,
    Gene_Name,
    CAM_Function,
    stat,
    pvalue,
    padj,
    Interaction_significant
  ) %>%
  arrange(
    CAM_Function,
    padj
  )

as.data.frame(core_cam_interaction_11)


# ============================================================
# 30. KEGG ENRICHMENT OF LEAF × TIME INTERACTION GENES
# ============================================================

library(clusterProfiler)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 30.1 Genes with significant leaf × time interaction
# ------------------------------------------------------------

interaction_gene_ids <- interaction_df %>%
  filter(
    !is.na(padj),
    padj < 0.05
  ) %>%
  pull(GeneID) %>%
  unique()


# Keep only genes represented in our KEGG universe
interaction_gene_ids_kegg <- base::intersect(
  interaction_gene_ids,
  enrichment_universe
)


cat(
  "\nSignificant interaction genes:",
  length(interaction_gene_ids),
  "\n"
)

cat(
  "Significant interaction genes with KEGG pathway mapping:",
  length(interaction_gene_ids_kegg),
  "\n"
)


# ------------------------------------------------------------
# 30.2 KEGG over-representation analysis
# ------------------------------------------------------------

interaction_kegg <- clusterProfiler::enricher(
  gene = interaction_gene_ids_kegg,
  universe = enrichment_universe,
  TERM2GENE = term2gene,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.20,
  minGSSize = 10,
  maxGSSize = 1000
)


interaction_kegg_df <- as.data.frame(
  interaction_kegg
)


cat(
  "\nSignificantly enriched pathways:",
  nrow(interaction_kegg_df),
  "\n"
)


# ------------------------------------------------------------
# 30.3 ADD KEGG PATHWAY NAMES
# ------------------------------------------------------------

interaction_kegg_named <- interaction_kegg_df %>%
  mutate(
    ID = as.character(ID)
  ) %>%
  left_join(
    pathway_lookup,
    by = "ID"
  ) %>%
  arrange(
    p.adjust
  )


# Inspect top results
interaction_kegg_named %>%
  select(
    ID,
    Pathway_Name,
    GeneRatio,
    BgRatio,
    pvalue,
    p.adjust,
    Count
  ) %>%
  head(30) %>%
  as.data.frame()


# ------------------------------------------------------------
# 30.4 SAVE
# ------------------------------------------------------------

write.csv(
  interaction_kegg_named,
  "KEGG_corrected_results/Leaf_Time_Interaction_KEGG_Enrichment.csv",
  row.names = FALSE
)

# ============================================================
# 30.5 CAM-RELEVANT INTERACTION PATHWAYS
# ============================================================

interaction_cam_pathways <- interaction_kegg_named %>%
  filter(
    grepl(
      paste(
        c(
          "Carbon fixation",
          "Carbon metabolism",
          "Starch and sucrose",
          "Glycolysis",
          "Pyruvate",
          "Circadian",
          "Photosynthesis",
          "Oxidative phosphorylation",
          "Plant hormone"
        ),
        collapse = "|"
      ),
      Pathway_Name,
      ignore.case = TRUE
    )
  ) %>%
  arrange(
    p.adjust
  )


interaction_cam_pathways %>%
  select(
    ID,
    Pathway_Name,
    GeneRatio,
    BgRatio,
    p.adjust,
    Count
  ) %>%
  as.data.frame()


write.csv(
  interaction_cam_pathways,
  "KEGG_corrected_results/Leaf_Time_Interaction_CAM_Pathways.csv",
  row.names = FALSE
)


# ============================================================
# 31. GENES DRIVING CAM-RELEVANT INTERACTION KEGG PATHWAYS
# ============================================================

library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 31.1 Extract genes contributing to the CAM-relevant pathways
# ------------------------------------------------------------

interaction_cam_drivers <- interaction_cam_pathways %>%
  select(
    ID,
    Pathway_Name,
    geneID,
    p.adjust
  ) %>%
  separate_rows(
    geneID,
    sep = "/"
  ) %>%
  rename(
    GeneID = geneID
  ) %>%
  distinct()


cat(
  "\nUnique genes driving CAM-relevant interaction pathways:",
  n_distinct(interaction_cam_drivers$GeneID),
  "\n"
)


# ------------------------------------------------------------
# 31.2 Add EggNOG annotation
# ------------------------------------------------------------

interaction_cam_drivers_annotated <- interaction_cam_drivers %>%
  left_join(
    anno_clean %>%
      select(
        GeneID,
        Gene_Name,
        Description,
        KO
      ) %>%
      distinct(GeneID, .keep_all = TRUE),
    by = "GeneID"
  ) %>%
  arrange(
    Pathway_Name,
    GeneID
  )


# Inspect
interaction_cam_drivers_annotated %>%
  select(
    Pathway_Name,
    GeneID,
    Gene_Name,
    Description,
    KO,
    p.adjust
  ) %>%
  head(100) %>%
  as.data.frame()


# Save
write.csv(
  interaction_cam_drivers_annotated,
  "KEGG_corrected_results/Interaction_CAM_Pathway_Driver_Genes.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 31.3 Find overlap with your 66 CAM-associated genes
# ------------------------------------------------------------

interaction_cam_driver_overlap <- interaction_cam_drivers_annotated %>%
  inner_join(
    core_cam_unique %>%
      select(
        GeneID,
        CAM_Function
      ) %>%
      distinct(),
    by = "GeneID"
  ) %>%
  arrange(
    CAM_Function,
    Pathway_Name
  )


cat(
  "\nCore CAM genes occurring in interaction-enriched pathways:",
  n_distinct(interaction_cam_driver_overlap$GeneID),
  "\n"
)


interaction_cam_driver_overlap %>%
  select(
    GeneID,
    Gene_Name,
    CAM_Function,
    Description,
    Pathway_Name,
    KO,
    p.adjust
  ) %>%
  as.data.frame()


write.csv(
  interaction_cam_driver_overlap,
  "KEGG_corrected_results/Core_CAM_Genes_in_Interaction_Pathways.csv",
  row.names = FALSE
)

# ============================================================
# 31.4 CORE CAM ENZYME OVERLAP
# ============================================================

core_enzyme_pathway_overlap <- interaction_cam_driver_overlap %>%
  filter(
    CAM_Function %in% c(
      "PEPC",
      "PPCK",
      "PPDK"
    )
  ) %>%
  arrange(
    CAM_Function,
    GeneID
  )

core_enzyme_pathway_overlap %>%
  select(
    GeneID,
    Gene_Name,
    CAM_Function,
    Description,
    Pathway_Name,
    KO
  ) %>%
  as.data.frame()


write.csv(
  core_enzyme_pathway_overlap,
  "KEGG_corrected_results/PEPC_PPCK_PPDK_Interaction_Pathway_Overlap.csv",
  row.names = FALSE
)

# ============================================================
# 32. INTEGRATED CORE CAM CANDIDATE TABLE
# ============================================================

library(dplyr)
library(tidyr)

# ------------------------------------------------------------
# 32.1 Collapse KEGG pathways to one row per gene
# ------------------------------------------------------------

cam_pathway_summary <- interaction_cam_driver_overlap %>%
  group_by(GeneID) %>%
  summarise(
    KEGG_Pathways = paste(
      sort(unique(Pathway_Name)),
      collapse = "; "
    ),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 32.2 Combine time-course + interaction + KEGG evidence
# ------------------------------------------------------------

final_cam_candidates <- cam_candidate_summary %>%
  
  left_join(
    cam_interaction_results %>%
      select(
        GeneID,
        Interaction_stat = stat,
        Interaction_pvalue = pvalue,
        Interaction_padj = padj,
        Interaction_significant
      ),
    by = "GeneID"
  ) %>%
  
  left_join(
    cam_pathway_summary,
    by = "GeneID"
  ) %>%
  
  mutate(
    
    In_CAM_enriched_KEGG = ifelse(
      !is.na(KEGG_Pathways),
      "Yes",
      "No"
    ),
    
    Interaction_significant = ifelse(
      !is.na(Interaction_padj) &
        Interaction_padj < 0.05,
      "Yes",
      "No"
    )
    
  )


# ------------------------------------------------------------
# 32.3 Create an evidence score
# ------------------------------------------------------------

final_cam_candidates <- final_cam_candidates %>%
  mutate(
    
    Evidence_score =
      
      # significant LP6 vs LP1 differences
      Significant_timepoints +
      
      # significant leaf × time interaction
      ifelse(
        Interaction_significant == "Yes",
        2,
        0
      ) +
      
      # membership of CAM-relevant enriched KEGG pathway
      ifelse(
        In_CAM_enriched_KEGG == "Yes",
        2,
        0
      )
  )


# ------------------------------------------------------------
# 32.4 Rank candidates
# ------------------------------------------------------------

final_cam_candidates <- final_cam_candidates %>%
  arrange(
    desc(Evidence_score),
    desc(Maximum_abs_log2FC)
  )


# ------------------------------------------------------------
# 32.5 Inspect top candidates
# ------------------------------------------------------------

final_cam_candidates %>%
  select(
    GeneID,
    Gene_Name,
    CAM_Function,
    Significant_timepoints,
    LP6_higher_timepoints,
    LP1_higher_timepoints,
    Mean_log2FC,
    Maximum_abs_log2FC,
    Interaction_padj,
    Interaction_significant,
    In_CAM_enriched_KEGG,
    KEGG_Pathways,
    Evidence_score
  ) %>%
  head(30) %>%
  as.data.frame()


# ------------------------------------------------------------
# 32.6 Save full candidate table
# ------------------------------------------------------------

write.csv(
  final_cam_candidates,
  "KEGG_corrected_results/Final_Core_CAM_Candidate_Ranking.csv",
  row.names = FALSE
)


cat(
  "\n========================================\n",
  "INTEGRATED CAM CANDIDATE ANALYSIS\n",
  "========================================\n",
  "Total CAM candidates:",
  nrow(final_cam_candidates),
  "\nInteraction significant:",
  sum(
    final_cam_candidates$Interaction_significant == "Yes"
  ),
  "\nIn CAM-relevant enriched KEGG pathways:",
  sum(
    final_cam_candidates$In_CAM_enriched_KEGG == "Yes"
  ),
  "\n========================================\n"
)

# ============================================================
# 32.7 PRIORITY PEPC / PPCK / PPDK CANDIDATES
# ============================================================

priority_cam_candidates <- final_cam_candidates %>%
  filter(
    CAM_Function %in% c(
      "PEPC",
      "PPCK",
      "PPDK"
    )
  ) %>%
  arrange(
    desc(Evidence_score),
    desc(Maximum_abs_log2FC)
  )


priority_cam_candidates %>%
  select(
    GeneID,
    Gene_Name,
    CAM_Function,
    Significant_timepoints,
    LP6_higher_timepoints,
    LP1_higher_timepoints,
    Mean_log2FC,
    Maximum_abs_log2FC,
    Interaction_padj,
    Interaction_significant,
    KEGG_Pathways,
    Evidence_score
  ) %>%
  as.data.frame()


write.csv(
  priority_cam_candidates,
  "KEGG_corrected_results/Priority_PEPC_PPCK_PPDK_Candidates.csv",
  row.names = FALSE
)


# ============================================================
# 33. DIFFERENTIALLY EXPRESSED GENES ACROSS THE DIEL CYCLE
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# Correct chronological order
time_order <- c("10", "14", "18", "22", "02", "06")


# ------------------------------------------------------------
# 33.1 COUNT SIGNIFICANT GENES AT EACH TIMEPOINT
# ------------------------------------------------------------

deg_summary <- data.frame()

for (tp in time_order) {
  
  res_df <- all_de_results[[tp]]
  
  # Significant DEGs using the same threshold as the rest
  # of the analysis
  sig <- res_df %>%
    filter(
      !is.na(padj),
      padj < 0.05,
      abs(log2FoldChange) > 1
    )
  
  # Positive log2FC = higher in LP6
  lp6_higher <- sum(
    sig$log2FoldChange > 0
  )
  
  # Negative log2FC = higher in LP1
  lp1_higher <- sum(
    sig$log2FoldChange < 0
  )
  
  deg_summary <- bind_rows(
    deg_summary,
    data.frame(
      Timepoint = tp,
      LP6_higher = lp6_higher,
      LP1_higher = lp1_higher,
      Total_DEGs = nrow(sig)
    )
  )
}


# ------------------------------------------------------------
# 33.2 DISPLAY RESULTS
# ------------------------------------------------------------

deg_summary$Timepoint <- factor(
  deg_summary$Timepoint,
  levels = time_order
)

as.data.frame(
  deg_summary
)

write.csv(
  deg_summary,
  "KEGG_corrected_results/DEG_Summary_Across_Diel_Cycle.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 33.5 PREPARE DATA FOR PLOTTING
# ------------------------------------------------------------

deg_plot_data <- deg_summary %>%
  select(
    Timepoint,
    LP6_higher,
    LP1_higher
  ) %>%
  pivot_longer(
    cols = c(
      LP6_higher,
      LP1_higher
    ),
    names_to = "Direction",
    values_to = "Genes"
  ) %>%
  mutate(
    
    Plot_value = ifelse(
      Direction == "LP1_higher",
      -Genes,
      Genes
    ),
    
    Direction = recode(
      Direction,
      LP6_higher = "Higher in LP6 (CAM)",
      LP1_higher = "Higher in LP1 (C3)"
    )
  )


# ------------------------------------------------------------
# 33.6 FINAL DISSERTATION DEG FIGURE
# ------------------------------------------------------------

p_deg_diel <- ggplot(
  deg_plot_data,
  aes(
    x = Timepoint,
    y = Plot_value,
    fill = Direction
  )
) +
  
  geom_col(
    width = 0.72
  ) +
  
  geom_hline(
    yintercept = 0,
    linewidth = 0.5
  ) +
  
  geom_text(
    aes(
      label = Genes
    ),
    position = position_stack(vjust = 0.5),
    size = 3.5
  ) +
  
  scale_y_continuous(
    labels = abs,
    expand = expansion(
      mult = c(0.08, 0.08)
    )
  ) +
  
  labs(
    title = "Differential gene expression between LP1 and LP6 across the diel cycle",
    
    subtitle = paste0(
      "DEGs defined as adjusted p < 0.05 and |log2FC| > 1"
    ),
    
    x = "Sampling time (h)",
    y = "Number of differentially expressed genes",
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    plot.subtitle = element_text(
      size = 10
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    legend.position = "top"
  )

print(
  p_deg_diel
)

ggsave(
  "KEGG_corrected_results/Figure_DEGs_Across_Diel_Cycle.png",
  p_deg_diel,
  width = 9,
  height = 6,
  units = "in",
  dpi = 600
)

ggsave(
  "KEGG_corrected_results/Figure_DEGs_Across_Diel_Cycle.pdf",
  p_deg_diel,
  width = 9,
  height = 6,
  units = "in",
  device = "pdf"
)


# ============================================================
# 34. KEGG PATHWAY × TIMEPOINT HEATMAP
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# 34.1 CHECK THE AVAILABLE DIRECTIONS
# ------------------------------------------------------------

unique(combined_kegg_named$Direction)

# ------------------------------------------------------------
# 34.2 KEEP DIRECTIONAL ENRICHMENT RESULTS
# ------------------------------------------------------------

kegg_directional <- combined_kegg_named %>%
  filter(
    Direction %in% c(
      "Higher_in_LP6_CAM",
      "Higher_in_LP1_C3"
    ),
    !is.na(Pathway_Name),
    !is.na(p.adjust),
    p.adjust < 0.05
  )

# Check how many significant pathways are present
table(
  kegg_directional$Timepoint,
  kegg_directional$Direction
)

# ------------------------------------------------------------
# 34.3 SELECT CAM-RELEVANT BIOLOGICAL PATHWAYS
# ------------------------------------------------------------

pathway_patterns <- c(
  "Carbon fixation",
  "Carbon metabolism",
  "Starch and sucrose",
  "Glycolysis",
  "Pyruvate",
  "Photosynthesis",
  "Circadian",
  "Oxidative phosphorylation",
  "Plant hormone",
  "Pentose phosphate"
)

kegg_cam_heatmap <- kegg_directional %>%
  filter(
    grepl(
      paste(
        pathway_patterns,
        collapse = "|"
      ),
      Pathway_Name,
      ignore.case = TRUE
    )
  )

sort(
  unique(
    kegg_cam_heatmap$Pathway_Name
  )
)

# ============================================================
# 34.4 DEFINE FINAL DISSERTATION PATHWAYS
# ============================================================

final_kegg_pathways <- c(
  "Carbon fixation by Calvin cycle",
  "Carbon metabolism",
  "Glycolysis / Gluconeogenesis",
  "Pentose phosphate pathway",
  "Starch and sucrose metabolism",
  "Photosynthesis",
  "Photosynthesis - antenna proteins",
  "Oxidative phosphorylation",
  "Circadian rhythm - plant",
  "Plant hormone signal transduction"
)

kegg_cam_final <- kegg_directional %>%
  filter(
    Pathway_Name %in% final_kegg_pathways
  )

# ------------------------------------------------------------
# 34.5 CHECK SIGNIFICANCE ACROSS TIME AND DIRECTION
# ------------------------------------------------------------

kegg_pathway_presence <- kegg_cam_final %>%
  count(
    Pathway_Name,
    Direction,
    Timepoint
  ) %>%
  pivot_wider(
    names_from = Timepoint,
    values_from = n,
    values_fill = 0
  ) %>%
  arrange(
    Direction,
    Pathway_Name
  )

as.data.frame(
  kegg_pathway_presence
)


# ============================================================
# 34.6 FINAL KEGG PATHWAY × DIEL TIME HEATMAP
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)

time_order <- c(
  "10", "14", "18", "22", "02", "06"
)

# ------------------------------------------------------------
# 34.6.1 CREATE SIGNED ENRICHMENT SCORE
# ------------------------------------------------------------

kegg_signed <- kegg_cam_final %>%
  mutate(
    
    Timepoint = factor(
      Timepoint,
      levels = time_order
    ),
    
    # -log10 adjusted p-value
    Enrichment_strength = -log10(p.adjust),
    
    # Positive = LP6/CAM
    # Negative = LP1/C3
    Signed_score = case_when(
      
      Direction == "Higher_in_LP6_CAM" ~
        Enrichment_strength,
      
      Direction == "Higher_in_LP1_C3" ~
        -Enrichment_strength,
      
      TRUE ~ NA_real_
    )
  )

# ------------------------------------------------------------
# 34.6.2 COLLAPSE TO ONE VALUE PER PATHWAY × TIME
# ------------------------------------------------------------

kegg_signed_summary <- kegg_signed %>%
  group_by(
    Pathway_Name,
    Timepoint
  ) %>%
  summarise(
    
    Signed_score = Signed_score[
      which.max(abs(Signed_score))
    ],
    
    .groups = "drop"
  )

# ------------------------------------------------------------
# 34.6.3 COMPLETE ALL PATHWAY × TIME COMBINATIONS
# ------------------------------------------------------------

kegg_signed_complete <- kegg_signed_summary %>%
  complete(
    Pathway_Name = final_kegg_pathways,
    Timepoint = factor(
      time_order,
      levels = time_order
    )
  )


# ------------------------------------------------------------
# 34.7 BIOLOGICAL PATHWAY ORDER
# ------------------------------------------------------------

pathway_order <- c(
  "Circadian rhythm - plant",
  "Photosynthesis - antenna proteins",
  "Photosynthesis",
  "Carbon fixation by Calvin cycle",
  "Carbon metabolism",
  "Glycolysis / Gluconeogenesis",
  "Pentose phosphate pathway",
  "Starch and sucrose metabolism",
  "Oxidative phosphorylation",
  "Plant hormone signal transduction"
)

kegg_signed_complete <- kegg_signed_complete %>%
  mutate(
    
    Pathway_Name = factor(
      Pathway_Name,
      levels = rev(pathway_order)
    ),
    
    Timepoint = factor(
      Timepoint,
      levels = time_order
    )
  )

# ------------------------------------------------------------
# 34.8 FINAL DISSERTATION KEGG HEATMAP
# ------------------------------------------------------------

p_kegg_final <- ggplot(
  kegg_signed_complete,
  aes(
    x = Timepoint,
    y = Pathway_Name,
    fill = Signed_score
  )
) +
  
  geom_tile(
    colour = "grey85",
    linewidth = 0.5
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    na.value = "grey95",
    name = expression(
      "Directional enrichment\n" ~ -log[10] * "(adjusted p)"
    )
  ) +
  
  labs(
    title = "Diel enrichment of pathways associated with the C3-to-CAM transition",
    
    subtitle = paste0(
      "Red: enriched among LP6-higher genes; ",
      "blue: enriched among LP1-higher genes"
    ),
    
    x = "Sampling time (h)",
    y = NULL,
    
    caption = paste0(
      "Significant enrichment: adjusted p < 0.05. ",
      "Grey cells indicate no significant directional enrichment."
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    plot.subtitle = element_text(
      size = 10,
      margin = margin(
        b = 10
      )
    ),
    
    plot.caption = element_text(
      size = 8,
      hjust = 0
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 11
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 10
    ),
    
    axis.text.y = element_text(
      size = 10
    ),
    
    legend.title = element_text(
      size = 9
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    panel.grid = element_blank(),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )

print(p_kegg_final)

# ------------------------------------------------------------
# 34.9 SAVED FIGURE
# ------------------------------------------------------------

ggsave(
  "KEGG_corrected_results/Figure_KEGG_Diel_Pathway_Heatmap.png",
  p_kegg_final,
  width = 11,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  "KEGG_corrected_results/Figure_KEGG_Diel_Pathway_Heatmap.pdf",
  p_kegg_final,
  width = 11,
  height = 7,
  units = "in",
  device = "pdf"
)

write.csv(
  kegg_signed_complete,
  "KEGG_corrected_results/Figure_KEGG_Diel_Pathway_Heatmap_Data.csv",
  row.names = FALSE
)

as.data.frame(
  kegg_signed_complete
)

# ============================================================
# 34.10 FINAL DISSO KEGG HEATMAP
# ============================================================

library(ggplot2)
library(scales)

max_score <- ceiling(
  max(
    abs(kegg_signed_complete$Signed_score),
    na.rm = TRUE
  )
)

p_kegg_final <- ggplot(
  kegg_signed_complete,
  aes(
    x = Timepoint,
    y = Pathway_Name,
    fill = Signed_score
  )
) +
  
  geom_tile(
    colour = "white",
    linewidth = 0.6
  ) +
  
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-max_score, max_score),
    oob = scales::squish,
    na.value = "grey92",
    name = expression(
      "Directional enrichment " ~ -log[10] * "(adjusted p)"
    )
  ) +
  
  labs(
    title = "Diel pathway enrichment associated with the C3-to-CAM transition",
    
    subtitle = paste0(
      "Red = enriched among LP6-higher genes; ",
      "blue = enriched among LP1-higher genes"
    ),
    
    x = "Sampling time (h)",
    y = NULL,
    
    caption = paste0(
      "Only significantly enriched pathways are shown ",
      "(adjusted p < 0.05). Grey cells indicate no significant enrichment."
    )
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    
    plot.subtitle = element_text(
      size = 10,
      margin = margin(b = 8)
    ),
    
    plot.caption = element_text(
      size = 8,
      hjust = 0
    ),
    
    axis.title.x = element_text(
      face = "bold",
      size = 11
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 10
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    legend.title = element_text(
      size = 9
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    panel.grid = element_blank()
  )

print(p_kegg_final)

ggsave(
  "KEGG_corrected_results/Figure_KEGG_Diel_Pathway_Heatmap_FINAL.png",
  p_kegg_final,
  width = 11,
  height = 7,
  units = "in",
  dpi = 600,
  bg = "white"
)

ggsave(
  "KEGG_corrected_results/Figure_KEGG_Diel_Pathway_Heatmap_FINAL.pdf",
  p_kegg_final,
  width = 11,
  height = 7,
  units = "in",
  device = "pdf",
  bg = "white"
)


# ============================================================
# KEGG MAPPER METABOLIC MAP INPUT FILES
# ============================================================
# Purpose:
# Prepare KO-level differential-expression files for KEGG Mapper
# to visualise global metabolic shifts between LP6 (mature CAM)
# and LP1 (young C3) across the diel cycle.
#
# Positive log2FC = higher in LP6 / mature CAM
# Negative log2FC = higher in LP1 / young C3

dir.create(
  "KEGG_Mapper_inputs",
  showWarnings = FALSE
)

for (tp in time_order) {
  
  # Get DESeq2 results for this timepoint
  res_df <- all_de_results[[tp]]
  
  # Keep significant differentially expressed genes
  sig <- res_df %>%
    filter(
      !is.na(padj),
      padj < 0.05,
      abs(log2FoldChange) > 1
    )
  
  # Add KEGG Orthology (KO) identifiers
  mapper_df <- sig %>%
    inner_join(
      gene_to_ko,
      by = "GeneID"
    ) %>%
    
    # Keep the information useful for interpretation
    select(
      GeneID,
      KO,
      log2FoldChange,
      padj
    ) %>%
    
    # If several genes map to the same KO,
    # retain the gene with the strongest expression difference
    arrange(
      desc(abs(log2FoldChange))
    ) %>%
    
    distinct(
      KO,
      .keep_all = TRUE
    ) %>%
    
    # Assign direction for KEGG Mapper
    mutate(
      Direction = ifelse(
        log2FoldChange > 0,
        "LP6_higher",
        "LP1_higher"
      ),
      
      Colour = ifelse(
        log2FoldChange > 0,
        "red",
        "blue"
      )
    )
  
  # ----------------------------------------------------------
  # Save a detailed table for reproducibility
  # ----------------------------------------------------------
  
  write.csv(
    mapper_df,
    file.path(
      "KEGG_Mapper_inputs",
      paste0(
        "KEGG_Mapper_",
        tp,
        "h_detailed.csv"
      )
    ),
    row.names = FALSE
  )
  
  # ----------------------------------------------------------
  # Create simple KO + colour file for KEGG Mapper
  # ----------------------------------------------------------
  
  mapper_input <- mapper_df %>%
    select(
      KO,
      Colour
    )
  
  write.table(
    mapper_input,
    file.path(
      "KEGG_Mapper_inputs",
      paste0(
        "KEGG_Mapper_",
        tp,
        "h.txt"
      )
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
  
  cat(
    tp,
    ":00 -",
    nrow(mapper_input),
    "KO identifiers written for KEGG Mapper\n"
  )
}


