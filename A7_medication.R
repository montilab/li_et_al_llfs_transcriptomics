#### Analysis model
## The analysis method and model is exact same as the aging signature, the only difference is 4 medications are added to the model.

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(Biobase))

#### parameters & logs ####
root <- ".../batch5_16k/"
output.dir <- paste0(root, "Aging_Signature/v1_frz6_allsubject_wmed/result_each_task/")
RNA_file <- paste0(root, "data/RNA_symbol.csv")
input.RNAseq <- paste0(root, "data/expdat_fltrnorm_v1.rds")
input.pheno <- paste0(root, "data/pheno_v1_wmed_1845samp_2503.rds")

input.gds <- "...LLFS.WGS.freeze6.chr21ck04.BI-SNP.gds"
grm <- as.matrix(readRDS(paste0(root, "data/grm_v1_2101samp_2503.rds")))

task_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
RNA_per_task <- 500

sink(paste0(output.dir, "task", task_id, ".log"), append = FALSE, split = TRUE)
date()

RNAall <- read.csv(RNA_file) %>% pull(RNA)
RNAall <- RNAall[((task_id - 1) * RNA_per_task + 1):(task_id * RNA_per_task)]
RNAall <- RNAall[complete.cases(RNAall)]

#### Read Data ####
pheno.dat <- readRDS(input.pheno)
RNAseq.dat <- as.data.frame(t(readRDS(input.RNAseq)[RNAall, ]))
RNAseq.dat$subject <- rownames(RNAseq.dat)
analysis.dat <- merge(x = RNAseq.dat, y = pheno.dat, by = "subject", all = F)

analysis.dat <- analysis.dat %>% mutate(
  batch = as.factor(batch),
  fc = relevel(as.factor(fc), ref = "BU"),
  Sex = relevel(as.factor(Sex), ref = "Male")
)

gds <- seqOpen(input.gds)
id.gds <- data.frame(sample.id = seqGetData(gds, "sample.id"))

#### Run Association ####
annot <- left_join(id.gds, analysis.dat, by = c("sample.id" = "subject"))
seqData <- SeqVarData(gds, sampleData = AnnotatedDataFrame(annot))

res <- c()

for (idx in c(1:length(RNAall))) {
  RNA0 <- RNAall[idx]
  cat("-- Start", idx, RNA0, "\n")

  nullmod <- try(fitNullModel(seqData,
    outcome = RNA0,
    sample.id = colnames(grm), ## need to subset seqData to subjects available in grm
    covars = c(
      "Age.enrollment", "fc", "Sex", "Education",
      "percent_intergenic", "batch",
      "htn", "lipid", "nitro", "diab",
      "PC1", "PC2", "PC3", "PC4",
      "PC5", "PC6", "PC7", "PC8",
      "PC9", "PC10"
    ),
    cov.mat = grm,
    family = "gaussian", verbose = T
  ))

  if (class(nullmod) != "try-error") {
    cat("   Sample size:", nrow(nullmod$fit), "\n")
    coeff <- nullmod$fixef

    var_names <- rownames(coeff)
    stat_names <- colnames(coeff)
    new_col_names <- paste(rep(var_names, each = length(stat_names)), stat_names, sep = "_")
    res0 <- data.frame(matrix(unlist(t(coeff)), nrow = 1))
    colnames(res0) <- new_col_names
    res0 <- res0 %>% mutate(RNA = RNA0, .before = everything())
    res <- rbind(res, res0)

    cat("-- Finish", idx, RNA0, "\n\n")
  } else {
    cat("-- Error:", idx, RNA0, "\n\n")
  }
}

#### output ####
saveRDS(res, paste0(output.dir, "coeff_task", task_id, ".rds"))
saveRDS(nullmod, paste0(output.dir, "nullmod_example_task", task_id, ".rds"))

closefn.gds(gds)

#### Plots ---------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ggvenn)

pheno <- readRDS("../../data/pheno_v1_wmed_1845samp_2503.rds") %>%
  select(
    subject, Age.enrollment, fc, Sex, Education,
    percent_intergenic, batch,
    htn, lipid, nitro, diab,
    PC1, PC2, PC3, PC4,
    PC5, PC6, PC7, PC8,
    PC9, PC10
  ) %>%
  filter(complete.cases(.)) %>%
  mutate(
    htn = as.factor(as.character(htn)),
    lipid = as.factor(as.character(lipid)),
    nitro = as.factor(as.character(nitro)),
    diab = as.factor(as.character(diab))
  )

venn_fill_color <- RColorBrewer::brewer.pal(4, "Pastel1")
venn_font_color <- c("coral3", "skyblue4", "darkseagreen4", "mediumpurple3")

res_o <- readRDS("../v1_frz6_allsubject/age_coeff_250407.rds") %>%
  select(RNA, symbol, starts_with("Age_"))
res_m <- readRDS("age_coeff_wmed_250613.rds") %>%
  select(RNA, symbol, starts_with("Age_"))

med_coeff <- readRDS("full_coeff_wmed_250613.rds") %>%
  select(
    RNA, symbol, starts_with("Age_"),
    starts_with("htn_"), starts_with("lipid_"), starts_with("nitro_"), starts_with("diab_")
  )

compare <- merge(
  x = res_o,
  y = res_m %>%
    rename_with(~ paste0(., "_med"), .cols = starts_with("Age_")) %>% mutate(symbol = NULL),
  by = "RNA", all = F
)

## summary table
top_table <- med_coeff %>%
  filter(htn_qval <= 0.01 | lipid_qval <= 0.01 | nitro_qval <= 0.01 | diab_qval <= 0.01) %>%
  select(RNA, symbol, Age_Z, Age_qval, htn_Z, htn_qval, lipid_Z, lipid_qval, nitro_Z, nitro_qval, diab_Z, diab_qval) %>%
  rename_with(~ sub("_qval", "_q", .), everything())

top_table <- merge(
  x = res_o %>% select(RNA, symbol, Age_Z, Age_qval) %>%
    rename(Age_Z_ori = Age_Z, Age_q_ori = Age_qval),
  y = top_table, by = c("RNA", "symbol"), all = F
) %>%
  arrange(lipid_q) %>%
  mutate_if(is.numeric, formatC, 4)

DT::datatable(top_table)

## compare age effect on each gene with vs without med
p1 <- ggplot(compare, aes(x = Age_Est, y = Age_Est_med)) +
  geom_point(alpha = 0.8, size = 1) +
  geom_abline(slope = 1, lwd = 0.7, intercept = 0, lty = 2, col = "deepskyblue3") +
  labs(title = "Age Effect", caption = "blue line: x=y", x = "Age Effect", y = "Age Effect with Medication") +
  theme_pubr() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14)
  )

p2 <- ggplot(compare, aes(x = -log10(Age_pval), y = -log10(Age_pval_med))) +
  geom_point(alpha = 0.8, size = 1) +
  geom_abline(slope = 1, lwd = 0.7, intercept = 0, lty = 2, col = "deepskyblue3") +
  labs(title = "Age -log10 p-value", caption = "blue line: x=y", x = "-log10 Age P", y = "-log10 Age P with Medication") +
  theme_pubr() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14)
  )

plot <- ggarrange(p1, p2, ncol = 2, nrow = 1)
annotate_figure(plot, top = text_grob("Compare Age Effect on Each RNA With vs Without Medication",
  face = "bold", size = 14
))

## effect of medications on gene expression

p1 <- ggplot(
  med_coeff %>% mutate(sig_symbol = ifelse(htn_qval < 0.01, yes = symbol, no = "")),
  aes(x = htn_Est, y = -log10(htn_pval), label = sig_symbol)
) +
  geom_point(alpha = 0.8, size = 1) +
  # geom_hline(yintercept = -log10(0.05 / 11229), lty = 2, col = "deepskyblue3") +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "seagreen3") +
  geom_vline(xintercept = 0, lty = 3, col = "grey") +
  geom_label_repel(box.padding = 0.5, point.padding = 0.2, segment.color = "deepskyblue3", min.segment.length = 0, max.overlaps = Inf) +
  ylim(0, 20) +
  labs(title = "Hypertension Treatment", x = "Effect on gene expression", y = "-log10 p-value") +
  theme_pubclean() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14)
  )

p2 <- ggplot(
  med_coeff %>% mutate(sig_symbol = ifelse(lipid_qval < 0.01, yes = symbol, no = "")),
  aes(x = lipid_Est, y = -log10(lipid_pval), label = sig_symbol)
) +
  geom_point(alpha = 0.8, size = 1) +
  # geom_hline(yintercept = -log10(0.05 / 11229), lty = 2, col = "deepskyblue3") +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "seagreen3") +
  geom_vline(xintercept = 0, lty = 3, col = "grey") +
  geom_label_repel(box.padding = 0.5, point.padding = 0.2, segment.color = "deepskyblue3", min.segment.length = 0, max.overlaps = Inf) +
  ylim(0, 20) +
  labs(title = "High Cholesterol Treatment", x = "Effect on gene expression", y = "-log10 p-value") +
  theme_pubclean() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14)
  )

p3 <- ggplot(
  med_coeff %>% mutate(sig_symbol = ifelse(nitro_qval < 0.01, yes = symbol, no = "")),
  aes(x = nitro_Est, y = -log10(nitro_pval), label = sig_symbol)
) +
  geom_point(alpha = 0.8, size = 1) +
  # geom_hline(yintercept = -log10(0.05 / 11229), lty = 2, col = "deepskyblue3") +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "seagreen3") +
  geom_vline(xintercept = 0, lty = 3, col = "grey") +
  geom_label_repel(box.padding = 0.5, point.padding = 0.2, segment.color = "deepskyblue3", min.segment.length = 0, max.overlaps = Inf) +
  ylim(0, 20) +
  labs(title = "Heart Disease Treatment (nitroglycerin)", x = "Effect on gene expression", y = "-log10 p-value") +
  theme_pubclean() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14)
  )

p4 <- ggplot(
  med_coeff %>% mutate(sig_symbol = ifelse(diab_qval < 0.01, yes = symbol, no = "")),
  aes(x = diab_Est, y = -log10(diab_pval), label = sig_symbol)
) +
  geom_point(alpha = 0.8, size = 1) +
  # geom_hline(yintercept = -log10(0.05 / 11229), lty = 2, col = "deepskyblue3") +
  geom_hline(yintercept = -log10(0.05), lty = 2, col = "seagreen3") +
  geom_vline(xintercept = 0, lty = 3, col = "grey") +
  geom_label_repel(box.padding = 0.5, point.padding = 0.2, segment.color = "deepskyblue3", min.segment.length = 0, max.overlaps = Inf) +
  ylim(0, 20) +
  labs(title = "Diabetes Treatment", x = "Effect on gene expression", y = "-log10 p-value") +
  theme_pubclean() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14)
  )

plot <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)

annotate_figure(plot,
  top = text_grob("Effect of Medications on Gene Expressions",
    face = "bold", size = 14
  ),
  bottom = text_grob("green line: p = 0.05",
    size = 12, hjust = 1, x = 1
  )
)
