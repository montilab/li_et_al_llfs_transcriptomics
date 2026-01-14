#### preprocessing and quality filtering ####

library(tidyverse)
library(DESeq2)
library(edgeR)

# read in data ----------------------------------------------------------------
dds <- readRDS("...../batch5_202311/dds_gene_20231023.rds")
design(dds) <- ~1

## excluding pool controls
dds <- dds[, dds$purpose == "experiment"]

# filter ----------------------------------------------------------------------

## filter applied at 2025/1/23:
## cpm > 3 in 3% of samples (125 samples) of the whole experiment
num_samples <- floor(ncol(dds) * .03) # 125
mid_expression_filter <- rowSums(cpm(counts(dds)) > 3) >= num_samples

sample_library_ids <- colData(dds) %>%
  as.data.frame() %>%
  filter(!suspicious_sex, relabel == mislabel, percent_intergenic < 0.08) %>%
  ## in the experiment set, a sample should either be mislabeled and relabeled, i.e. both==1
  ## or not mislabeled and not relabeled, i.e. both==0

  ## RN7SL1 and RN7SL2
  mutate(total_minus_rn7sl = protein_coding_total - rn7sl_total) %>%
  mutate(total_minus_rn7sl_quantile = ecdf(.$total_minus_rn7sl)(total_minus_rn7sl)) %>%
  group_by(subject, visitcode) %>%
  arrange(percent_intergenic, .by_group = TRUE) %>%
  mutate(
    total_minus_rn7sl_max_diff =
      total_minus_rn7sl[percent_intergenic == dplyr::first(percent_intergenic)][1] -
        total_minus_rn7sl[total_minus_rn7sl == max(total_minus_rn7sl)][1]
  ) %>%
  mutate(
    select_flag = if_else(
      condition = total_minus_rn7sl_max_diff > -1e6 | dplyr::first(total_minus_rn7sl_quantile) >= 0.2,
      true = percent_intergenic == min(percent_intergenic),
      false = total_minus_rn7sl == max(total_minus_rn7sl)
    )
  ) %>%
  filter(select_flag) %>%
  pull(library_id)

## Apply filters ---------------------------------------------------------------
# filter the deseq data object on the rows (genes) and samples (column)
dds_filtered <- dds[mid_expression_filter, dds$library_id %in% sample_library_ids]

# add size factors to the dds object
dds_filtered <- estimateSizeFactors(dds_filtered)

# extract visit 1 samples ------------------------------------------------------
sample_v1 <- dds_filtered$visitcode %in% c(1, 4, 7)
dds_filtered_v1 <- dds_filtered[, sample_v1]

# extract data & log2 transform & make it to gene x sample matrix --------------
norm_fltr_counts_v1 <- counts(dds_filtered_v1, normalized = TRUE)
rownames(norm_fltr_counts_v1) <- rownames(dds_filtered_v1)
colnames(norm_fltr_counts_v1) <- dds_filtered_v1$subject
norm_fltr_counts_v1 <- log2(norm_fltr_counts_v1 + 1)

# write ------------------------------------------------------------------------
saveRDS(as.data.frame(colData(dds_filtered_v1)), "RNAseq_sample_meta_v1.rds")
saveRDS(norm_fltr_counts_v1, "expdat_fltrnorm_v1.rds")
