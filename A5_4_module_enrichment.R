library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(hypeR)

RNAexpr <- readRDS("../../data/expdat_fltrnorm_v1.rds")

# subset to the subjects analyzed
samples_used <- readRDS("../../data/frz6_allcvrt_resid/residual_2503.rds") %>% pull(subject)
RNAexpr <- RNAexpr[, samples_used]
RNAsymbol <- read.csv("../../data/RNA_symbol.csv")


## gene sets
hallmark <- msigdb_gsets("Homo sapiens", "H", "", clean = TRUE)
kegg <- msigdb_gsets("Homo sapiens", "C2", "KEGG_LEGACY", clean = TRUE)
reactome <- msigdb_gsets("Homo sapiens", "C2", "REACTOME", clean = TRUE)

module_assignment <- readRDS("module.rds")

## manually assign module sequence
module_name_list <- c(
  "turquoise", "green", "black", "magenta", # up
  "blue", "brown", "yellow", "red", "pink" # dn
)
module_assignment <- module_assignment[module_name_list]


## Change ensembl id to symbol.
mod_symbol <- list()
for (i in c(seq(length(module_assignment)))) {
  mod_symbol[[i]] <- RNAsymbol[which(RNAsymbol$RNA %in% module_assignment[[i]]), "symbol"]
}
names(mod_symbol) <- names(module_assignment)

#### Enrichment ----------------------------------------------------------------
h <- hypeR(mod_symbol, hallmark, background = 19116, test = "hypergeometric", fdr = 0.05, plotting = FALSE, quiet = TRUE)
k <- hypeR(mod_symbol, kegg, background = 19116, test = "hypergeometric", fdr = 0.05, plotting = FALSE, quiet = TRUE)
r <- hypeR(mod_symbol, reactome, background = 19116, test = "hypergeometric", fdr = 0.05, plotting = FALSE, quiet = TRUE)

#### Plotting ------------------------------------------------------------------
hyp_dots(h, merge = TRUE, top = 20, fdr = 0.05, title = "HALLMARK") +
  theme(text = element_text(size = 15))
hyp_dots(k, merge = TRUE, top = 20, fdr = 0.05, title = "KEGG") +
  theme(text = element_text(size = 15))
hyp_dots(r, merge = TRUE, top = 20, fdr = 0.05, title = "REACTOME") +
  theme(text = element_text(size = 15))
