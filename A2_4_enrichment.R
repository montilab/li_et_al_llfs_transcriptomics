library(dplyr)
library(ggplot2)
library(ggpubr)
library(fgsea)
library(hypeR)

res <- readRDS("age_coeff_250407.rds")
allGene <- read.csv("../../data/RNA_annot_all.csv") %>%
  pull(hgnc_symbol) %>%
  unique()

#### Hypergeom enrichment ------------------------------------------------------
thresh_fdr <- 0.01
thresh_Z <- 6

sig <- list(
  "up" = res %>%
    filter(symbol != "", Age_qval < thresh_fdr, Age_Z > thresh_Z) %>%
    arrange(Age_qval) %>% distinct(symbol) %>%
    pull(symbol),
  "dn" = res %>%
    filter(symbol != "", Age_qval < thresh_fdr, Age_Z < (-1) * thresh_Z) %>%
    arrange(Age_qval) %>% distinct(symbol) %>%
    pull(symbol)
)

#### Fetch gene sets -----------------------------------------------------------
## background set to 19116 for the total number of genes
hallmark <- hypeR::msigdb_gsets("Homo sapiens", "H", "", clean = TRUE)
kegg <- hypeR::msigdb_gsets("Homo sapiens", "C2", "KEGG_LEGACY", clean = TRUE)
reactome <- hypeR::msigdb_gsets("Homo sapiens", "C2", "REACTOME", clean = TRUE)

GOBP <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:BP", clean = TRUE)
GOCC <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:CC", clean = TRUE)
GOMF <- hypeR::msigdb_gsets("Homo sapiens", "C5", "GO:MF", clean = TRUE)

#### Perform enrichment --------------------------------------------------------
sig_hypergeo_H <- hypeR(
  signature = sig, genesets = hallmark, test = "hypergeometric",
  fdr = 0.05, background = 19116
)
sig_hypergeo_K <- hypeR(
  signature = sig, genesets = kegg, test = "hypergeometric",
  fdr = 0.05, background = 19116
)
sig_hypergeo_R <- hypeR(
  signature = sig, genesets = reactome, test = "hypergeometric",
  fdr = 0.05, background = 19116
)

sig_hypergeo_GOBP <- hypeR(
  signature = sig, genesets = GOBP, test = "hypergeometric",
  fdr = 0.05, background = 19116
)
sig_hypergeo_GOCC <- hypeR(
  signature = sig, genesets = GOCC, test = "hypergeometric",
  fdr = 0.05, background = 19116
)
sig_hypergeo_GOMF <- hypeR(
  signature = sig, genesets = GOMF, test = "hypergeometric",
  fdr = 0.05, background = 19116
)


### Plots ----------------------------------------------------------------------
## Only Gene Sets at the top and reached FDR \< 0.05 are included in the plot.
hyp_dots(sig_hypergeo_H, merge = TRUE, top = 30, title = "HALLMARK", size_by = "none", fdr = 0.01) + theme(text = element_text(size = 13))

hyp_dots(sig_hypergeo_K, merge = TRUE, top = 30, title = "KEGG", size_by = "none", fdr = 0.01) + theme(text = element_text(size = 13))

hyp_dots(sig_hypergeo_R, merge = TRUE, top = 30, title = "REACTOME", size_by = "none", fdr = 0.01) + theme(text = element_text(size = 13))

hyp_dots(sig_hypergeo_GOBP, merge = TRUE, top = 30, abrv = 55, title = "GOBP", size_by = "none", fdr = 0.01) + theme(text = element_text(size = 13))

hyp_dots(sig_hypergeo_GOCC, merge = TRUE, top = 30, abrv = 55, title = "GOCC", size_by = "none", fdr = 0.01) + theme(text = element_text(size = 13))

hyp_dots(sig_hypergeo_GOMF, merge = TRUE, top = 30, abrv = 55, title = "GOMF", size_by = "none", fdr = 0.01) + theme(text = element_text(size = 13))


#### SenMayo KS test -----------------------------------------------------------
# Paper: https://www.nature.com/articles/s41467-022-32552-1
# Gene set downloaded from: https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/SAUL_SEN_MAYO

SenMayo_gs <- list("SenMayo" = unique(c("ACVR1B", "ANG", "ANGPT1", "ANGPTL4", "AREG", "AXL", "BEX3", "BMP2", "BMP6", "C3", "CCL1", "CCL13", "CCL16", "CCL2", "CCL20", "CCL24", "CCL26", "CCL3", "CCL3L1", "CCL4", "CCL5", "CCL7", "CCL8", "CD55", "CD9", "CSF1", "CSF2", "CSF2RB", "CST4", "CTNNB1", "CTSB", "CXCL1", "CXCL10", "CXCL12", "CXCL16", "CXCL2", "CXCL3", "CXCL8", "CXCR2", "DKK1", "EDN1", "EGF", "EGFR", "EREG", "ESM1", "ETS2", "FAS", "FGF1", "FGF2", "FGF7", "GDF15", "GEM", "GMFG", "HGF", "HMGB1", "ICAM1", "ICAM3", "IGF1", "IGFBP1", "IGFBP2", "IGFBP3", "IGFBP4", "IGFBP5", "IGFBP6", "IGFBP7", "IL10", "IL13", "IL15", "IL18", "IL1A", "IL1B", "IL2", "IL32", "IL6", "IL6ST", "IL7", "INHA", "IQGAP2", "ITGA2", "ITPKA", "JUN", "KITLG", "LCP1", "MIF", "MMP1", "MMP10", "MMP12", "MMP13", "MMP14", "MMP2", "MMP3", "MMP9", "NAP1L4", "NRG1", "PAPPA", "PECAM1", "PGF", "PIGF", "PLAT", "PLAU", "PLAUR", "PTBP1", "PTGER2", "PTGES", "RPS6KA5", "SCAMP4", "SELPLG", "SEMA3F", "SERPINB4", "SERPINE1", "SERPINE2", "SPP1", "SPX", "TIMP2", "TNF", "TNFRSF10C", "TNFRSF11B", "TNFRSF1A", "TNFRSF1B", "TUBGCP2", "VEGFA", "VEGFC", "VGF", "WNT16", "WNT2")))

ranked_signature_dat <- res %>%
  filter(symbol != "") %>%
  arrange(desc(Age_Z))
ranked_signature_Z <- ranked_signature_dat %>% pull(Age_Z)
names(ranked_signature_Z) <- ranked_signature_dat %>% pull(symbol)

fgseaRes <- fgsea(
  pathways = SenMayo_gs, # gene set list, can get this by `hypeR::msigdb_gsets("Homo sapiens", "H", "")$genesets`
  stats = ranked_signature_Z, # numeric vector of score, the `names` are gene symbols
  minSize = 5,
  maxSize = 500
)

#### KS plotting function ------------------------------------------------------
## original reference from fgsea function:
## https://rdrr.io/bioc/fgsea/src/R/plot.R

plotEnrichment <- function(pathway, stats,
                           gseaParam = 1,
                           ticksSize = 0.2, line_color = "blue4") {
  rnk <- rank(-stats)
  ord <- order(rnk)

  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))

  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)

  gseaRes <- calcGseaStat(statsAdj,
    selectedStats = pathway,
    returnAllExtremes = TRUE
  )

  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops

  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))

  diff <- (max(tops) - min(bottoms)) / 8

  # Getting rid of NOTEs
  x <- y <- NULL
  g <- ggplot(toPlot, aes(x = x, y = y)) +
    geom_point(color = line_color, size = 0.1) +
    geom_hline(yintercept = max(tops), colour = "firebrick2", linetype = "dashed") +
    geom_hline(yintercept = min(bottoms), colour = "firebrick2", linetype = "dashed") +
    geom_hline(yintercept = 0, colour = "black") +
    geom_line(color = line_color) +
    geom_segment(
      data = data.frame(x = pathway),
      mapping = aes(
        x = x, y = -diff / 2,
        xend = x, yend = diff / 2
      ),
      lwd = ticksSize
    ) +
    theme(
      panel.border = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(x = "rank", y = "enrichment score")
  g
}

## KS plot
plotEnrichment(SenMayo_gs$SenMayo, ranked_signature_Z, ticksSize = 0.3) +
  labs(title = "SenMayo") +
  theme_pubclean() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, color = "blue4"),
    plot.caption = element_text(size = 11, color = "blue4")
  ) + geom_line(color = "blue4")
