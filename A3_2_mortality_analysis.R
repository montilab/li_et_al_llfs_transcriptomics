## use residual after regress out batch covariates

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(survival))

#### Read data -----------------------------------------------------------------
covariates <- c(
  "Age.enrollment", "Sex", "Education", "fc", "percent_intergenic",
  "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"
)

pheno_dat <- readRDS("../../data/pheno_v1_2167samp_2503.rds") %>%
  select(subject, all_of(covariates), pedid, Alive, Age.last.contact) %>%
  mutate(Sex = relevel(Sex, ref = "Male")) %>%
  mutate(Age_diff = Age.last.contact - Age.enrollment, .after = Age.last.contact) %>%
  mutate(status = ifelse(Alive == "Yes", yes = 0, no = 1), .after = Alive)

pheno_dat$fam_id <- 0
unique_pedid_v1 <- unique(pheno_dat$pedid)
for (i in seq(length(unique_pedid_v1))) {
  pheno_dat$fam_id[which(pheno_dat$pedid == unique_pedid_v1[i])] <- i
}

RNAall <- read.csv("../../data/RNA_symbol.csv") %>% pull(RNA)
RNA_symbol <- read.csv("../../data/RNA_symbol.csv")

## read residual after regress on batch
expr <- readRDS("../../data/frz6_batch_resid/lm_batch_resid_2503.rds")

## normalize expression matrix so the input value is standardized
## should be: rows are genes and columns are samples
rownames(expr) <- expr$subject
expr$subject <- NULL
expr <- t(expr)
standardized_expr <- apply(expr, 1, function(x) (x - mean(x)) / sd(x))
standardized_expr <- as.data.frame(standardized_expr)
standardized_expr$subject <- rownames(standardized_expr)

#### Prepare data --------------------------------------------------------------
analysis_dat <- merge(
  x = pheno_dat,
  y = standardized_expr,
  by = "subject"
)
analysis_dat_no_missing <- na.omit(analysis_dat)

## GEE analysis of each RNA ----------------------------------------------------
res <- c()

for (RNA0 in RNAall) {
  fix.eff <- paste("Surv(Age_diff, status) ~", RNA0)
  if (!is.null(covariates)) {
    for (covi in covariates) fix.eff <- paste(fix.eff, "+", covi)
  }
  fix.eff <- paste(fix.eff, "+cluster(fam_id)")
  fix.eff <- formula(fix.eff)
  fit <- try(coxph(fix.eff, data = analysis_dat_no_missing))

  coeff <- coefficients(summary(fit)) %>%
    as.data.frame() %>%
    rename(
      eff = coef, HR = `exp(coef)`, se = `se(coef)`,
      robust_se = `robust se`, Z = z, pval = `Pr(>|z|)`
    )
  rownames(coeff)[1] <- "RNA"
  rownames(coeff)[which(rownames(coeff) == "Age.enrollment")] <- "AgeEnrollment"

  var_names <- rownames(coeff)
  stat_names <- colnames(coeff)
  new_col_names <- paste(rep(var_names, each = length(stat_names)), stat_names, sep = "_")
  res0 <- data.frame(matrix(unlist(t(coeff)), nrow = 1))
  colnames(res0) <- new_col_names
  res0 <- res0 %>% mutate(RNA = RNA0, .before = everything())

  res <- rbind(res, res0)
}

#### Add gene symbol and FDR ---------------------------------------------------
res <- res %>% rename_with(~ sub("^RNA_", "Mortality_", .), starts_with("RNA_"))
res <- res %>%
  add_column(Mortality_qval = p.adjust(res$Mortality_pval, method = "BH"), .after = "Mortality_pval")

res <- merge(
  x = RNA_symbol, by.y = "RNA", all.x = F,
  y = res, by.x = "RNA", all.y = T
) %>% arrange(Mortality_pval)

#### Save ----------------------------------------------------------------------

## remove the genes with PAR_Y
## remove the genes that failed to converge in aging signature model
age_RNA <- readRDS("../../Aging_Signature/v1_frz6_allsubject/age_coeff_250407.rds") %>% pull(RNA)
res <- res %>% filter(RNA %in% age_RNA)

saveRDS(res, paste0("survl_full_v1_", format(Sys.Date(), format = "%y%m%d"), ".rds"))
saveRDS(res %>% select(
  RNA, symbol, Mortality_eff, Mortality_HR, Mortality_robust_se, Mortality_Z, Mortality_pval, Mortality_qval
), paste0("survl_RNA_v1_", format(Sys.Date(), format = "%y%m%d"), ".rds"))
