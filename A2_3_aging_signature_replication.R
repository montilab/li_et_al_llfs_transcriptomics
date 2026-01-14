library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)

#### Read data -----------------------------------------------------------------
## For duplicated symbols, take the more significant one (smaller pval)
LLFS <- readRDS("age_coeff_250407.rds") %>%
  select(RNA, symbol, starts_with("Age")) %>%
  filter(symbol != "") %>%
  group_by(symbol) %>%
  arrange(Age_pval) %>%
  slice(1) %>%
  ungroup() %>%
  rename_with(~ sub("^Age_", "LLFS_", .), .cols = !c(RNA, symbol)) %>%
  as.data.frame()

NatC <- read.csv("...transcriptional_landscape_of_age/supplementary_data_1.csv") %>%
  mutate(NatC_q = p.adjust(NatC_pval, method = "BH"), .after = "NatC_pval") %>%
  group_by(NatC_symbol) %>%
  filter(NatC_pval == min(NatC_pval)) %>%
  ungroup() %>%
  as.data.frame()

compare_shared <- merge(x = LLFS, y = NatC, by.x = "symbol", by.y = "NatC_symbol", all = F)

#### Z-score comparison plot ---------------------------------------------------

p1 <- ggplot(compare_shared, aes(x = LLFS_Z, y = NatC_Z)) +
  geom_point(alpha = 0.7, size = 0.8) +
  geom_smooth(method = "lm", formula = y ~ x) +
  stat_cor() +
  stat_regline_equation(label.y = 20) +
  theme_pubr() +
  theme(
    text = element_text(family = "Helvetica"),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 10), legend.title = element_text(size = 12)
  ) +
  labs(x = "LLFS Z-scores", y = "Peters et.al. Z-scores", caption = paste0("Matched symbols: ", length(unique(compare_shared$symbol))))

#### Venn Diagram --------------------------------------------------------------
color <- RColorBrewer::brewer.pal(4, "Pastel1")
font_color <- c("coral3", "skyblue4", "darkseagreen4", "mediumpurple3")

### Filter abs(Z)>3 & q<=0.01 from each

set_up <- list(
  "LLFS up" = LLFS %>% filter(LLFS_Z > 3, LLFS_qval <= 0.01) %>% pull(symbol) %>% unique(),
  "Peters up" = NatC %>% filter(NatC_Z > 3, NatC_q <= 0.01) %>% pull(NatC_symbol) %>% unique()
)
p2 <- ggvenn::ggvenn(
  set_up,
  stroke_size = 0.5,
  set_name_size = 5, set_name_color = font_color[1:2],
  text_size = 5,
  show_percentage = F,
  fill_color = color[1:2],
) + coord_cartesian(ylim = c(-1, 2)) +
  labs(caption = paste0("Fisher's P < 2.2e-16\nBackground = ", length(genes_total_union))) +
  theme(plot.caption = element_text(size = 10))

set_dn <- list(
  "LLFS down" = LLFS %>% filter(LLFS_Z < -3, LLFS_qval <= 0.01, symbol != "") %>% pull(symbol) %>% unique(),
  "Peters down" = NatC %>% filter(NatC_Z < -3, NatC_q <= 0.01) %>% pull(NatC_symbol) %>% unique()
)
p3 <- ggvenn::ggvenn(
  set_dn,
  stroke_size = 0.5,
  set_name_size = 5, set_name_color = font_color[3:4],
  text_size = 5,
  show_percentage = F,
  fill_color = color[3:4]
) + coord_cartesian(ylim = c(-1, 2)) +
  labs(caption = paste0("Fisher's P < 2.2e-16\nBackground = ", length(genes_total_union))) +
  theme(plot.caption = element_text(size = 10))

ggarrange(
  p1,
  ggarrange(p2, p3, ncol = 1, nrow = 2),
  nrow = 1, ncol = 2,
  widths = c(4, 3)
)
