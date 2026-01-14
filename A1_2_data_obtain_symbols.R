library("biomaRt")
library(dplyr)

#### read the original dds object before QC ------------------------------------
dat <- readRDS(".../data/dds_analysis_gene.rds")

ensemblid <- rownames(dat)

# get rid of the version number after the dot
ensemblid_s <- unlist(strsplit(ensemblid, "\\."))
ensemblid_s <- ensemblid_s[grep("ENSG", ensemblid_s)]

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

annot <- getBM(
  mart = mart,
  attributes = c(
    "ensembl_gene_id",
    "hgnc_symbol",
    "gene_biotype",
    "external_gene_name",
    "uniprot_gn_symbol",
    "uniprot_gn_id"
  ),
  values = ensemblid_s,
  uniqueRows = TRUE
)

dat <- merge(
  x = data.frame("ensembl_full" = ensemblid, "ensembl" = ensemblid_s), by.x = "ensembl",
  y = annot, by.y = "ensembl_gene_id",
  all.x = T, all.y = F
)

dat$hgnc_symbol[which(is.na(dat$hgnc_symbol))] <- ""
dat <- dat %>% distinct()

sum(unique(dat$hgnc_symbol) != "")
# 40189 unique symbols before QC filtering

write.csv(dat, "RNA_annot_all.csv", row.names = F)

#### obtain gene symbols in the analyzed data ----------------------------------
RNAs <- rownames(readRDS("expdat_fltrnorm_v1.rds"))
dat_hgnc_analyzed <- dat %>%
  filter(ensembl_full %in% RNAs) %>%
  select(ensembl_full, hgnc_symbol) %>%
  rename(RNA = ensembl_full, symbol = hgnc_symbol) %>%
  distinct()

write.csv(dat_hgnc_analyzed, "RNA_symbol.csv", row.names = F, quote = F)
