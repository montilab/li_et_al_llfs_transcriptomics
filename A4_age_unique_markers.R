library(dplyr)
library(tidyverse)
library(ggplot2)
library(hypeR)

#### Read analysis results -----------------------------------------------------

## LLFS & ILO
res <- read.csv("../../Compare_ILO/ILO_LLFS_Age_CombinedResults.csv")

## Peters
## for duplicated symbols, only keep the most significant one
Peters <- read.csv("../../../Gene_Signature_Compare/transcriptional_landscape_of_age/supplementary_data_1_table_useful_columns.csv") %>%
  mutate(Peters_qval = p.adjust(Peters_pval, method = "BH"), .after = "Peters_pval") %>%
  arrange(Peters_pval) %>%
  distinct(Peters_symbol, .keep_all = TRUE)

## Mortality risk result
survl <- readRDS("../../Survival/v1_batch_resid_frz6/survl_RNA_v1_250514.rds") %>%
  select(RNA, Mortality_HR, Mortality_pval) %>%
  rename(mortality_llfs_HR = Mortality_HR, mortality_llfs_pval = Mortality_pval)

compare <- merge(x = res, y = Peters, by.x = "symbol_llfs", by.y = "Peters_symbol", all = T)

## Unique aging markers: both LLFS+ILO age q<=0.01 & Peters p > 0.25
compare %>%
  filter(fdr_llfs <= 0.01, fdr_ilo <= 0.01, sign(est_llfs) == sign(est_ilo), !is.na(Peters_pval)) %>%
  nrow()

age_uniq_df <- compare %>%
  filter(fdr_llfs <= 0.01, fdr_ilo <= 0.01, sign(est_llfs) == sign(est_ilo), Peters_pval > 0.25)

sum_df <- data.frame(matrix(numeric(0),
  ncol = 4, nrow = 1,
  dimnames = list(
    c("number"), c("RNA_up", "RNA_dn", "Gene_Symbol_up", "Gene_Symbol_dn")
  )
))
sum_df <- data.frame()
sum_df["number", "RNA_up"] <- age_uniq_df %>%
  filter(est_llfs > 0) %>%
  pull(RNA) %>%
  unique() %>%
  length()
sum_df["number", "RNA_dn"] <- age_uniq_df %>%
  filter(est_llfs < 0) %>%
  pull(RNA) %>%
  unique() %>%
  length()

sum_df["number", "Gene_Symbol_up"] <- age_uniq_df %>%
  filter(symbol_llfs != "", est_llfs > 0) %>%
  pull(symbol_llfs) %>%
  unique() %>%
  length()
sum_df["number", "Gene_Symbol_dn"] <- age_uniq_df %>%
  filter(symbol_llfs != "", est_llfs < 0) %>%
  pull(symbol_llfs) %>%
  unique() %>%
  length()

sum_df <- sum_df %>%
  mutate(RNA_all = RNA_up + RNA_dn, .after = RNA_dn) %>%
  mutate(Gene_Symbol_all = Gene_Symbol_up + Gene_Symbol_dn, .after = Gene_Symbol_dn)

DT::datatable(sum_df)


#### Check mortality risk of those markers -------------------------------------
age_uniq_df <- merge(x = age_uniq_df, y = survl, by = "RNA", all.x = T, all.y = F) %>%
  mutate(mortality_llfs_qval = p.adjust(mortality_llfs_pval, method = "BH"), .after = "mortality_llfs_pval")

sig_all <- list(
  "LLFS_up" = age_uniq_df %>% filter(est_llfs > 0) %>% pull(symbol_llfs) %>% unique(),
  "LLFS_dn" = age_uniq_df %>% filter(est_llfs < 0) %>% pull(symbol_llfs) %>% unique(),
  "HR>1" = age_uniq_df %>% filter(mortality_llfs_HR > 1, mortality_llfs_qval < 0.05) %>% pull(symbol_llfs) %>% unique(),
  "HR<1" = age_uniq_df %>% filter(mortality_llfs_HR < 1, mortality_llfs_qval < 0.05) %>% pull(symbol_llfs) %>% unique()
)

#### Plot ----------------------------------------------------------------------
color <- RColorBrewer::brewer.pal(4, "Pastel1")
font_color <- c("coral3", "skyblue4", "darkseagreen4", "mediumpurple3")

ggvenn::ggvenn(
  sig_all,
  stroke_size = 0.5,
  set_name_size = 5, set_name_color = font_color,
  text_size = 5,
  show_percentage = F,
  fill_color = color
) + coord_cartesian(ylim = c(-2, 1.5), xlim = c(-2.5, 2))
