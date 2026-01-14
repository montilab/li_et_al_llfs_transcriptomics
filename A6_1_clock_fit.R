library(tidyverse)
library(ggplot2)
library(ggpubr)
library(glmnet)

## Read Data -------------------------------------------------------------------
pheno <- readRDS("../../data/pheno_v1_2167samp_2503.rds") %>%
  select(subject, Age.enrollment) %>%
  filter(complete.cases(.))

## read ags result just to get the 16284 genes
genes <- readRDS("../../Aging_Signature/v1_frz6_allsubject/age_coeff_250407.rds") %>% pull(RNA)

## read expr and re-arrange the column as the same order as pheno
expr <- readRDS("../../data/expdat_fltrnorm_v1.rds")[genes, pheno$subject]

#### Elastic Net Regression ----------------------------------------------------
## Prepare input
standardized_expr <- apply(expr, 1, function(x) (x - mean(x)) / sd(x))
y <- pheno$Age.enrollment

## check subject match
## x: standardized_expr
## y: pheno$Age.enrollment
all.equal(pheno$subject, rownames(standardized_expr))

## Try different alpha
all_model_alpha <- seq(0.05, 0.95, 0.05)
all_model_rmse <- c()
all_model_mae <- c()
all_model_lambda <- c()
all_model_list <- list()

for (idx in seq(all_model_alpha)) {
  alpha <- all_model_alpha[idx]

  set.seed(104)
  cv_enet <- cv.glmnet(standardized_expr, y, alpha = alpha)
  best_lambda <- cv_enet$lambda.min
  final_enet_model <- glmnet(standardized_expr, y, alpha = alpha, lambda = best_lambda)
  all_model_list[[idx]] <- final_enet_model
  all_model_lambda[idx] <- best_lambda

  predict <- as.numeric(predict(final_enet_model, standardized_expr))
  all_model_rmse[idx] <- caret::RMSE(predict, y)
  all_model_mae[idx] <- caret::MAE(predict, y)
}

all_model_statistics <- data.frame("idx" = seq(all_model_alpha), "alpha" = all_model_alpha, "lambda" = all_model_lambda, "RMSE" = all_model_rmse, "MAE" = all_model_mae)
DT::datatable(all_model_statistics)

## Manually select the model with lowest RMSE.
best_idx <- all_model_statistics %>%
  arrange(RMSE) %>%
  pull(idx) %>%
  head(1)
all_model_statistics %>% filter(idx == best_idx)

final_enet_model <- all_model_list[[best_idx]]
coefficients <- coef(final_enet_model)

saveRDS(final_enet_model, "A01_ENet_model_v1.rds")
saveRDS(coefficients, "A01_ENet_coeff_v1.rds")


#### Predicted age (uncorrected) -----------------------------------------------
dat <- data.frame(
  "subject" = pheno$subject,
  "age" = pheno$Age.enrollment,
  "predict_age" = as.numeric(predict(final_enet_model, standardized_expr))
)

ggplot(dat, aes(x = age, y = predict_age)) +
  geom_point(size = 1.3, alpha = 0.7) +
  geom_smooth(method = "lm", lwd = 1, formula = y ~ x) +
  stat_regline_equation(color = "blue2") +
  geom_abline(slope = 1, intercept = 0, lwd = 0.9, col = "firebrick2") +
  labs(caption = "red: x=y; blue: regression") +
  theme_linedraw() +
  theme(text = element_text(size = 14))

saveRDS(dat, "A01_ENet_predicted_v1.rds")
