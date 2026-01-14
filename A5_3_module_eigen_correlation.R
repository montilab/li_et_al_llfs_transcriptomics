library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)

mod_assign_a <- readRDS("module.rds")
genes_cor_a <- readRDS("average_genes_cor_with_eigen.rds")

ags <- readRDS("../../Aging_Signature/v1_frz6_allsubject/age_coeff_250407.rds")

mod_assign_a_df <- map2_dfr(.x = names(mod_assign_a), .y = mod_assign_a, ~ data.frame(mod = .x, RNA = .y))
ags_a <- merge(x = ags, y = mod_assign_a_df, all.x = T, all.y = T, by = "RNA") %>%
  mutate(mod = ifelse(is.na(mod), "grey", mod))

#### Select the markers correlated with eigen.
cor_cutoff <- 0.75
mod_assign_a_df <- map2_dfr(.x = names(mod_assign_a), .y = mod_assign_a, ~ data.frame(mod = .x, RNA = .y))
ags_a <- merge(x = ags, y = mod_assign_a_df, all.x = F, all.y = T, by = "RNA")
ags_a <- merge(x = ags_a, y = genes_cor_a %>% select(!mod), by = "RNA")

ags_a_label <- ags_a %>%
  filter(cor >= cor_cutoff) %>%
  group_by(mod) %>%
  mutate(label = ifelse(cor == max(cor), yes = paste0(mod, " - ", symbol), no = "")) %>%
  ungroup() %>%
  as.data.frame()


### Volcano plot ---------------------------------------------------------------
ggplot(ags_a_label, aes(x = Age_Est, y = -log10(Age_pval))) +
  geom_point(shape = 19, alpha = 0.6, size = 2.5, aes(color = mod)) +
  scale_color_identity() +
  geom_vline(xintercept = 0, lty = 3, col = "grey70") +
  geom_label_repel(aes(label = label),
    family = "Helvetica", size = 4, label.padding = 0.3, point.padding = 0.3,
    force = 20, segment.color = "grey60", max.overlaps = Inf,
  ) +
  labs(y = "-log10(p)", x = "Age Effect") +
  theme_pubclean() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold")
  ) +
  # coord_cartesian(xlim = c(-0.036, 0.036), ylim = c(-5, 85))
  coord_cartesian(xlim = c(-0.06, 0.06), ylim = c(-5, 85))
