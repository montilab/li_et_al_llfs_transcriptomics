suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(vip))
suppressPackageStartupMessages(library(glmnet))

#### Read data -----------------------------------------------------------------
model <- readRDS("A01_ENet_model_v1.rds")
ags <- readRDS("../../Aging_Signature/v1_frz6_allsubject/age_coeff_250407.rds")
pheno <- readRDS("../../data/pheno_v1_2167samp_2503.rds") %>%
  select(subject, pedid, Age.enrollment, Sex, Education, fc, PC1, PC2, PC3, PC4, batch, percent_intergenic, Alive, Age.last.contact)
pheno <- pheno[complete.cases(pheno), ]
expr <- readRDS("../../data/expdat_fltrnorm_v1.rds")[ags$RNA, pheno$subject]

standardized_expr <- apply(expr, 1, function(x) (x - mean(x)) / sd(x))
y <- pheno$Age.enrollment

#### Feature importance using permutation --------------------------------------
cat("-- Start VI. \n")
Sys.time()

vip_allsubject <- vi(model,
  method = "permute", target = pheno$Age.enrollment, train = standardized_expr,
  metric = "MAE", nsim = 10,
  pred_wrapper = function(model, newdata) as.vector(predict(model, newdata))
)

cat("-- Finished VI. \n")
Sys.time()

## Formatting the VI output
vi <- vip_allsubject %>%
  as.data.frame() %>%
  rename(RNA = Variable)

saveRDS(vi, "A04_ENet_Variable_Importance.rds")
