library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)

vi <- readRDS("A05_ENet_Variable_Importance.rds")
ags <- readRDS("../../Aging_Signature/v1_frz6_allsubject/age_coeff_250407.rds")

## add gene symbol to top important markers
vi_imp_label_cutoff <- vi %>%
  arrange(desc(Importance)) %>%
  head(10) %>%
  pull(Importance)
vi_dat <- merge(x = vi, y = ags, by = "RNA", all = F) %>%
  mutate(vi_imp_label = ifelse(Importance >= vi_imp_label_cutoff, yes = symbol, no = "")) %>%
  mutate(
    ags_sig_color = "grey50",
    ags_sig_color = ifelse(Age_qval <= 0.01 & Age_Z > 0, yes = "lightcoral", no = ags_sig_color),
    ags_sig_color = ifelse(Age_qval <= 0.01 & Age_Z < 0, yes = "cadetblue3", no = ags_sig_color),
    ags_sig_color = factor(ags_sig_color)
  ) %>%
  arrange(desc(Importance))

sum_df <- data.frame(matrix(numeric(0),
  nrow = 5, ncol = 2,
  dimnames = list(c("vi>0.1", "vi>0.05", "vi>0.01", "vi>0", "vi_0"), c("total", "aging_markers"))
))

sum_df["vi>0.1", "total"] <- sum(vi_dat$Importance >= 0.1)
sum_df["vi>0.05", "total"] <- sum(vi_dat$Importance >= 0.05)
sum_df["vi>0.01", "total"] <- sum(vi_dat$Importance >= 0.01)
sum_df["vi>0", "total"] <- sum(vi_dat$Importance > 0)
sum_df["vi_0", "total"] <- sum(vi_dat$Importance == 0)

sum_df["vi>0.1", "aging_markers"] <- sum(vi_dat$Importance >= 0.1 & vi_dat$Age_qval <= 0.01)
sum_df["vi>0.05", "aging_markers"] <- sum(vi_dat$Importance >= 0.05 & vi_dat$Age_qval <= 0.01)
sum_df["vi>0.01", "aging_markers"] <- sum(vi_dat$Importance >= 0.01 & vi_dat$Age_qval <= 0.01)
sum_df["vi>0", "aging_markers"] <- sum(vi_dat$Importance > 0 & vi_dat$Age_qval <= 0.01)
sum_df["vi_0", "aging_markers"] <- sum(vi_dat$Importance == 0 & vi_dat$Age_qval <= 0.01)

sum_df <- sum_df %>% mutate(aging_marker_percentage = paste0(round(aging_markers / total * 100, 1), "%"), non_aging_markers = total - aging_markers)

sum_df

#### Plot top important markers. Colored by aging marker. ----------------------
label_cutoff_rows <- 15
vi_dat_plot <- vi_dat %>%
  arrange(desc(Importance)) %>%
  head(label_cutoff_rows)

ggplot(
  vi_dat_plot,
  aes(x = reorder(symbol, Importance), y = Importance)
) +
  geom_segment(aes(x = symbol, xend = symbol, y = 0, yend = Importance), color = "slategray", alpha = 0.5) +
  geom_point(color = vi_dat_plot$ags_sig_color, size = 4, alpha = 0.9) +
  geom_text(aes(label = sprintf("%.2f", Importance)), hjust = -0.3, vjust = 0.5, color = "black", size = 4) +
  ylim(c(0, max(vi_dat$Importance) + 0.1)) +
  labs(x = "Gene Symbol", y = "Importance") +
  coord_flip() +
  theme_linedraw() +
  theme(
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.border = element_blank()
  )
