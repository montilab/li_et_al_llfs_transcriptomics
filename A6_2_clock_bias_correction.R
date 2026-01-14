library(tidyverse)
library(limma)
library(ggplot2)
library(ggpubr)

pred <- readRDS("A01_ENet_predicted_v1.rds") %>%
  rename(CA = age, PA = predict_age) %>%
  select(subject, CA, PA)

all_coeff <- as.matrix(readRDS("A01_ENet_coeff_v1.rds"))
coeff <- all_coeff[2:nrow(all_coeff), ]

## Histogram of non-zero coefficients
ggplot(data.frame("coefficients" = coeff) %>% filter(coefficients != 0), aes(x = coefficients)) +
  geom_histogram(aes(y = after_stat(count)), color = "skyblue2", fill = "skyblue2", bins = 500) +
  geom_density(aes(y = after_stat(density * n * (max(x) - min(x)) / 500)), lwd = 0.8, lty = 1, color = "royalblue", alpha = 0) +
  labs(title = "Histogram of non-zero coefficients") +
  theme_pubclean() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, color = "blue4"),
    plot.caption = element_text(size = 11, color = "blue4")
  )

#### Model ---------------------------------------------------------------------
# CA <- Chronological Age
# RA <- Regression Age
# PA <- Predicted Age
# CPA <- Corrected Predicted Age

# CA is the data. PA comes from gene exp.
# Regress PA ~ CA, gives us a regression line, which is RA.
# Ideally, PA = CA with individual noise. But now we see a systematic bias.

## Before Correction
p0 <- ggplot(pred, aes(x = CA, y = PA)) +
  geom_point(alpha = 0.7, size = 0.85) +
  geom_smooth(method = "lm", lwd = 0.8, formula = y ~ x) +
  stat_cor(label.y = 105) +
  stat_regline_equation(label.y = 95, color = "mediumblue") +
  geom_abline(slope = 1, intercept = 0, lwd = 0.8, lty = 2, col = "firebrick2") +
  labs(
    title = "Predicted Age Before Correction",
    caption = "red dash line: x=y \nblue line: regression line",
    y = "Predicted Age", x = "Chronological Age"
  ) +
  theme_pubr() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, color = "blue4"),
    plot.caption = element_text(size = 11, color = "blue4")
  )
p0


## Method 1 and 2 require the same parameters and have the same level of variance.
## Method 3 exclude CA to obtain CPA after obtain alpha and beta using M2, which introduces more variance but less biased by CA. \

# Offset_Beheshti <- PA - CA
# Offset_Beheshti = $\alpha$ × CA + $\beta$
# CPA = PA - ($\alpha$ × CA + $\beta$)

# Detailed steps:
# 1. Offset_Beheshti = PA - RA
# 2. fit linear regression line: Offset_Beheshti ~ CA, obtain $\alpha$ and $\beta$
# 3. CPA = PA - ($\alpha$ × CA + $\beta$)

# **Method 2: Lange**
# PA = $\alpha$ × CA + $\beta$
# CPA = PA + [CA - ($\alpha$ × CA + $\beta$)]

# Detailed steps:
# 1. fit linear regression line: PA ~ CA, obtain $\alpha$ and $\beta$
# 2. CPA = PA + [CA - ($\alpha$ × CA + $\beta$)]
#
# **Method 3: Cole**
# PA = $\alpha$ × CA + $\beta$ (same as method 2)
# CPA = (PA - $\beta$) / $\alpha$
# This way CPA is calculated without CA, which introduces more variance.

#### M1 Beheshti ---------------------------------------------------------------

# Offset_Beheshti <- PA - CA
# Offset_Beheshti = alpha × CA + beta
# CPA = PA - (alpha × CA + beta)

# 1. Offset_Beheshti = PA - CA
dat_m1 <- pred %>% mutate(offset = PA - CA)

# 2. fit linear regression line: Offset_Beheshti ~ CA, obtain alpha and beta \
m1 <- lm(offset ~ CA, dat_m1)
beta_m1 <- m1$coefficients[1]
alpha_m1 <- m1$coefficients[2]
cat("alpha:", alpha_m1, "\n", "beta:", beta_m1)

# 3. CPA = PA - (alpha × CA + beta) \
dat_m1 <- dat_m1 %>% mutate(CPA = PA - (alpha_m1 * CA + beta_m1))


# 4. visualize: CPA ~ CA <- ideally, "randomly distributed"
p1 <- ggplot(dat_m1, aes(x = CA, y = CPA)) +
  geom_point(alpha = 0.7, size = 0.85) +
  geom_smooth(method = "lm", lwd = 0.8, formula = y ~ x) +
  stat_cor(label.y = 105) +
  stat_regline_equation(label.y = 95, color = "mediumblue") +
  geom_abline(slope = 1, intercept = 0, lwd = 0.8, lty = 2, col = "firebrick2") +
  labs(
    title = "Beheshti Corrected Predicted Age",
    caption = "red dash line: x=y \nblue line: regression line",
    y = "Corrected Predicted Age", x = "Chronological Age"
  ) +
  theme_pubr() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, color = "blue4"),
    plot.caption = element_text(size = 11, color = "blue4")
  )
p1

dat_m1 <- dat_m1 %>%
  select(subject, CA, PA, CPA) %>%
  rename_with(.fn = ~ paste0(., "_v1"), .cols = !subject)
saveRDS(dat_m1, "A02_ENet_Beheshti_CPA_v1.rds")

#### M2 Lange ------------------------------------------------------------------
# PA = alpha × CA + beta
# CPA = PA + [CA - (alpha × CA + beta)]

# 1. fit linear regression line: PA ~ CA, obtain alpha and beta
dat_m2 <- pred
m2 <- lm(PA ~ CA, dat_m2)
beta_m2 <- m2$coefficients[1]
alpha_m2 <- m2$coefficients[2]
cat("alpha:", alpha_m2, "\n", "beta:", beta_m2)

# 2. CPA = PA + [CA - (alpha × CA + beta)]
dat_m2 <- dat_m2 %>% mutate(CPA = PA + (CA - (alpha_m2 * CA + beta_m2)))

# 2. visualize: CPA ~ CA <- ideally, "randomly distributed"
p2 <- ggplot(dat_m2, aes(x = CA, y = CPA)) +
  geom_point(alpha = 0.7, size = 0.85) +
  geom_smooth(method = "lm", lwd = 0.8, formula = y ~ x) +
  stat_cor(label.y = 105) +
  stat_regline_equation(label.y = 95, color = "mediumblue") +
  geom_abline(slope = 1, intercept = 0, lwd = 0.8, lty = 2, col = "firebrick2") +
  labs(
    title = "Lange Corrected Predicted Age",
    caption = "red dash line: x=y \nblue line: regression line",
    y = "Corrected Predicted Age", x = "Chronological Age"
  ) +
  theme_pubr() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, color = "blue4"),
    plot.caption = element_text(size = 11, color = "blue4")
  )
p2

dat_m2 <- dat_m2 %>%
  select(subject, CA, PA, CPA) %>%
  rename_with(.fn = ~ paste0(., "_v1"), .cols = !subject)
saveRDS(dat_m2, "A02_ENet_Lange_CPA_v1.rds")

#### M3 Cole -------------------------------------------------------------------
## Use the alpha and beta from M2 Lange.
# CPA = (PA - beta) / alpha
# alpha and beta are the same as method 2

alpha_m3 <- alpha_m2
beta_m3 <- beta_m2
cat("alpha:", alpha_m3, "\n", "beta:", beta_m3)
dat_m3 <- pred %>% mutate(CPA = (PA - beta_m3) / alpha_m3)

# 2. visualize: CPA ~ CA <- ideally, "randomly distributed"
p3 <- ggplot(dat_m3, aes(x = CA, y = CPA)) +
  geom_point(alpha = 0.7, size = 0.85) +
  geom_smooth(method = "lm", lwd = 0.8, formula = y ~ x) +
  stat_cor(label.y = 105) +
  stat_regline_equation(label.y = 95, color = "mediumblue") +
  geom_abline(slope = 1, intercept = 0, lwd = 0.8, lty = 2, col = "firebrick2") +
  labs(
    title = "Cole Corrected Predicted Age",
    caption = "red dash line: x=y \nblue line: regression line",
    y = "Corrected Predicted Age", x = "Chronological Age"
  ) +
  theme_pubr() +
  theme(
    text = element_text(family = "Helvetica", size = 14),
    axis.text = element_text(size = 13),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 14, color = "blue4"),
    plot.caption = element_text(size = 11, color = "blue4")
  )
p3

dat_m3 <- dat_m3 %>%
  select(subject, CA, PA, CPA) %>%
  rename_with(.fn = ~ paste0(., "_v1"), .cols = !subject)
saveRDS(dat_m3, "A02_ENet_Cole_CPA_v1.rds")
