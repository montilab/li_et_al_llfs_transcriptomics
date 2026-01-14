library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(hypeR)

res <- readRDS("survl_RNA_v1_250514.rds")

## "up" - "increase mortality", markers HR>1, Est>0, increased value is worse survival.\
## "dn" - "decrease mortality", markers HR<1, Est<0, increased value is better survival.\

thresh_fdr <- 0.01

sig <- list(
  "up" = res %>%
    filter(symbol != "", Mortality_qval <= thresh_fdr, Mortality_HR > 1) %>%
    arrange(Mortality_qval) %>% distinct(symbol) %>%
    pull(symbol),
  "dn" = res %>%
    filter(symbol != "", Mortality_qval <= thresh_fdr, Mortality_HR < 1) %>%
    arrange(Mortality_qval) %>% distinct(symbol) %>%
    pull(symbol)
)

#### Fetch Gene Sets and perform Hypergeometric Enrichment ---------------------
## Background: 19116 (total number of protein-coding genes)\
geneset_H <- hypeR::msigdb_gsets(species = "Homo sapiens", collection = "H", subcollection = "", clean = TRUE)
geneset_K <- hypeR::msigdb_gsets("Homo sapiens", "C2", "CP:KEGG_LEGACY", clean = TRUE)
geneset_R <- hypeR::msigdb_gsets("Homo sapiens", "C2", "CP:REACTOME", clean = TRUE)

sig_hypergeo_H <- hypeR(signature = sig, genesets = geneset_H, test = "hypergeometric", background = 19116)
sig_hypergeo_K <- hypeR(signature = sig, genesets = geneset_K, test = "hypergeometric", background = 19116)
sig_hypergeo_R <- hypeR(signature = sig, genesets = geneset_R, test = "hypergeometric", background = 19116)

### Plots
hyp_dots(sig_hypergeo_H, merge = TRUE, top = 30, fdr = 0.05, title = "HALLMARK") + theme(text = element_text(size = 14))
hyp_dots(sig_hypergeo_K, merge = TRUE, top = 30, fdr = 0.05, title = "KEGG") + theme(text = element_text(size = 14))
hyp_dots(sig_hypergeo_R, merge = TRUE, top = 30, fdr = 0.05, title = "REACTOME") + theme(text = element_text(size = 14))
