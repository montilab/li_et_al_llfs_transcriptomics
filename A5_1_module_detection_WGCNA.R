library(dplyr)
library(tidyverse)
library(ggplot2)

source("../utils/wgcna_wrapper.R")
source("../utils/wgcna_create_modules.R")

ags <- readRDS("../../Aging_Signature/v1_frz6_allsubject/age_coeff_250407.rds") %>%
  filter(abs(Age_Z) > 6, Age_qval <= 0.01) %>%
  pull(RNA)

# need sample x gene matrix
dat <- readRDS("../../data/frz6_allcvrt_resid/residual_2503.rds") %>%
  column_to_rownames(var = "subject")
dat <- dat[, ags]

wgcna_average <- wgcna.wrapper(
  dat = dat, # sample x gene matrix
  min.size = 20,
  cores = 4,
  cor.fn = c("bicor"),
  merging = TRUE,
  merging.cut = 0.1,
  cor.type = "signed",
  hclust.method = "average",
  do.plot = TRUE
)

mod_average <- createModules(dat = dat, wgcnaRes = wgcna_average)

saveRDS(wgcna_average, "wgcna_result.rds")
saveRDS(mod_average$mods, "module.rds")
