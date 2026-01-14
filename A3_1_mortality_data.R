## create residual of batches before running survival analysis

library(tidyverse)
library(limma)

exp <- readRDS("../expdat_fltrnorm_v1.rds")
pheno <- readRDS("../pheno_v1_1856samp_2501.rds")
RNAsymbol <- read.csv("../RNA_symbol.csv")
RNA <- rownames(exp)

## Run a basic lm to obtain residual of the batch
## RNA ~ batch

analysis_dat <- merge(
  x = pheno %>% select(subject, batch),
  y = exp %>% t() %>% as.data.frame() %>% rownames_to_column(var = "subject"),
  by = "subject", all = F
)

coef <- c()
resid <- data.frame("subject" = analysis_dat$subject)

for (idx in seq(RNA)) {
  RNA0 <- RNA[idx]
  formula0 <- paste0(RNA0, "~ batch")
  m0 <- lm(formula0, analysis_dat)
  s0 <- summary(m0)
  coef0 <- s0$coefficients %>%
    as.data.frame() %>%
    rownames_to_column(var = "batch") %>%
    mutate(RNA = RNA0, .before = everything())
  resid0 <- s0$residuals
  resid[, RNA0] <- resid0
  coef <- rbind(coef, coef0)
}

saveRDS(coef, "lm_batch_coeff_2501.rds")
saveRDS(resid, "lm_batch_resid_2501.rds")
