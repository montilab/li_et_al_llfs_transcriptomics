suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(hypeR))

## Read data -------------------------------------------------------------------
## manually edit the input and output file at the end ##
module_assignment <- readRDS("module.rds")

RNAexpr <- readRDS("../../data/expdat_fltrnorm_v1.rds")

# subset to the subjects analyzed
samples_used <- readRDS("../../data/frz6_allcvrt_resid/residual_2503.rds") %>% pull(subject)
RNAexpr <- RNAexpr[, samples_used]

RNAsymbol <- read.csv("../../data/RNA_symbol.csv")


## Hub Genes -------------------------------------------------------------------
### Extract Hub genes
hub_input_RNAexpr <- t(RNAexpr[unlist(module_assignment), ])
hub_input_module_length <- unlist(lapply(module_assignment, length))
hub_input_module_assignment <- rep(names(hub_input_module_length), hub_input_module_length)

module_hub <- WGCNA::chooseTopHubInEachModule(datExpr = hub_input_RNAexpr, hub_input_module_assignment)
module_hub_df <- data.frame("mod" = names(module_hub), "RNA" = module_hub)
module_hub_df <- merge(x = module_hub_df, y = RNAsymbol, by = "RNA") %>%
  select(mod, RNA, symbol) %>%
  rename(hub_RNA = RNA, hub_symbol = symbol)


## Top Geneset Annotation ------------------------------------------------------

## change ensembl id to symbol

mod_symbol <- list()
for (i in c(seq(length(module_assignment)))) {
  mod_symbol[[i]] <- RNAsymbol[which(RNAsymbol$RNA %in% module_assignment[[i]]), "symbol"]
}
names(mod_symbol) <- names(module_assignment)

## gene sets
hallmark <- msigdb_gsets("Homo sapiens", "H", "", clean = TRUE)
kegg <- msigdb_gsets("Homo sapiens", "C2", "CP:KEGG", clean = TRUE)
reactome <- msigdb_gsets("Homo sapiens", "C2", "CP:REACTOME", clean = TRUE)
GOMF <- msigdb_gsets("Homo sapiens", "C5", "GO:MF", clean = TRUE)
GOBP <- msigdb_gsets("Homo sapiens", "C5", "GO:BP", clean = TRUE)
GOCC <- msigdb_gsets("Homo sapiens", "C5", "GO:CC", clean = TRUE)

## annotate module with top enriched gene set
annot_df <- data.frame()

for (i in c(seq(length(mod_symbol)))) {
  cat("-- Start module", i, "\n")
  h0 <- hypeR(mod_symbol[[i]], hallmark, background = 19116, test = "hypergeometric", fdr = 0.05, plotting = FALSE, quiet = TRUE)
  k0 <- hypeR(mod_symbol[[i]], kegg, background = 19116, test = "hypergeometric", fdr = 0.05, plotting = FALSE, quiet = TRUE)
  r0 <- hypeR(mod_symbol[[i]], reactome, background = 19116, test = "hypergeometric", fdr = 0.05, plotting = FALSE, quiet = TRUE)
  GOBP0 <- hypeR(mod_symbol[[i]], GOBP, background = 19116, test = "hypergeometric", fdr = 0.05, plotting = FALSE, quiet = TRUE)
  GOMF0 <- hypeR(mod_symbol[[i]], GOMF, background = 19116, test = "hypergeometric", fdr = 0.05, plotting = FALSE, quiet = TRUE)
  GOCC0 <- hypeR(mod_symbol[[i]], GOCC, background = 19116, test = "hypergeometric", fdr = 0.05, plotting = FALSE, quiet = TRUE)

  # fill no-result with NA, so each df is 1 row instead of 0 row, so they can be cbind together
  for (df0 in c(h0, k0, r0, GOBP0, GOMF0, GOCC0)) {
    if (nrow(df0$data) == 0) df0$data[1, ] <- NA
  }

  annot_df0 <- cbind(
    "mod" = names(mod_symbol)[i],
    "size" = length(mod_symbol[[i]]),
    h0$data %>% select(label, pval, fdr, hits) %>% arrange(fdr) %>% slice_head(n = 1) %>% rename_with(~ paste0("annot_HALLMARK_", .)),
    k0$data %>% select(label, pval, fdr, hits) %>% arrange(fdr) %>% slice_head(n = 1) %>% rename_with(~ paste0("annot_KEGG_", .)),
    r0$data %>% select(label, pval, fdr, hits) %>% arrange(fdr) %>% slice_head(n = 1) %>% rename_with(~ paste0("annot_REACTOME_", .)),
    GOBP0$data %>% select(label, pval, fdr, hits) %>% arrange(fdr) %>% slice_head(n = 1) %>% rename_with(~ paste0("annot_GOBP_", .)),
    GOMF0$data %>% select(label, pval, fdr, hits) %>% arrange(fdr) %>% slice_head(n = 1) %>% rename_with(~ paste0("annot_GOMF_", .)),
    GOCC0$data %>% select(label, pval, fdr, hits) %>% arrange(fdr) %>% slice_head(n = 1) %>% rename_with(~ paste0("annot_GOCC_", .))
  )

  annot_df <- rbind(annot_df, annot_df0)
}

rownames(annot_df) <- NULL

## Save ------------------------------------------------------
res <- merge(x = module_hub_df, y = annot_df, by = "mod") %>% relocate(size, .after = "mod")

saveRDS(res, "module_annot.rds")
