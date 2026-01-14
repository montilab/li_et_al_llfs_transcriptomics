library(tidyverse)
library(ggplot2)
library(ggpubr)
library(adjustedCurves)
library(survival)

CPA <- readRDS("A02_ENet_Cole_CPA_v1.rds") %>%
  rename_with(~ sub("_v1", "", .), everything()) %>%
  mutate(Age_Accel = CPA - CA)

covariates <- c("AgeEnrollment", "Sex", "Education", "fc", "percent_intergenic", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

#### Age_Accel as continuous ---------------------------------------------------

## Read phenotype data, keep only rows with covariates information
pheno_dat <- readRDS("../../data/pheno_v1_2167samp_2503.rds") %>%
  select(subject, pedid, Sex, fc, Education, Alive, Age.enrollment, Age.last.contact, PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10, percent_intergenic) %>%
  filter(complete.cases(.)) %>%
  mutate(Sex = relevel(Sex, ref = "Male")) %>%
  mutate(Age_diff = Age.last.contact - Age.enrollment, .after = Age.last.contact) %>%
  mutate(status = ifelse(Alive == "Yes", yes = 0, no = 1), .after = Alive)

pheno_dat$fam_id <- 0
unique_pedid_v1 <- unique(pheno_dat$pedid)
for (i in seq(length(unique_pedid_v1))) {
  pheno_dat$fam_id[which(pheno_dat$pedid == unique_pedid_v1[i])] <- i
}

dat <- merge(
  x = pheno_dat,
  y = CPA,
  by = "subject", all = F
) %>%
  rename(AgeEnrollment = Age.enrollment)

rownames(dat) <- dat$subject

#### Run CoxPH -----------------------------------------------------------------
formula1 <- paste("Surv(Age_diff, status) ~ Age_Accel +", paste(covariates, collapse = " + "), "+ cluster(fam_id)")
formula1

formula1 <- formula(formula1)

model1 <- coxph(formula1, data = dat)

coeff1 <- coefficients(summary(model1)) %>%
  as.data.frame() %>%
  rename(
    eff = coef, HR = `exp(coef)`, se = `se(coef)`,
    robust_se = `robust se`, Z = z, pval = `Pr(>|z|)`
  )

## coefficients
DT::datatable(coeff1 %>% mutate_if(is.numeric, formatC, 4),
  options = list(pageLength = 15)
)

#### Age_Accel as Group --------------------------------------------------------
cutoff <- 3
dat_group <- dat %>%
  mutate(
    Age_Accel_Group = ifelse(Age_Accel > cutoff, yes = "Accel",
      no = ifelse(Age_Accel < (0 - cutoff), yes = "Decel", no = "Normal")
    ),
    Age_Accel_Group = factor(Age_Accel_Group, levels = c("Normal", "Decel", "Accel"))
  )

table(dat_group$Age_Accel_Group)
## Normal  Decel  Accel
##   1284    386    387

#### Run CoxPH -----------------------------------------------------------------
formula1_g <- paste("Surv(Age_diff, status) ~ Age_Accel_Group +", paste(covariates, collapse = " + "), "+ cluster(fam_id)")
formula1_g

formula1_g <- formula(formula1_g)

model1_g <- coxph(formula1_g, data = dat_group, x = TRUE)

coeff1_g <- coefficients(summary(model1_g)) %>%
  as.data.frame() %>%
  rename(
    eff = coef, HR = `exp(coef)`, se = `se(coef)`,
    robust_se = `robust se`, Z = z, pval = `Pr(>|z|)`
  )

## coefficient
DT::datatable(coeff1_g %>% mutate_if(is.numeric, formatC, 4),
  options = list(pageLength = 15)
)

#### Adjusted KM curve ---------------------------------------------------------
model1_g_km <- adjustedsurv(
  data = dat_group,
  variable = "Age_Accel_Group",
  ev_time = "Age_diff",
  event = "status",
  method = "direct",
  outcome_model = model1_g,
  conf_int = TRUE,
  bootstrap = FALSE
)

plot(model1_g_km, conf_int = TRUE) +
  theme_pubr() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, color = "blue4"),
    plot.caption = element_text(size = 11, color = "blue4"),
    legend.position = "right"
  )
