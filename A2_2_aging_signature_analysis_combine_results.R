suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))

output.dir <- "result_each_task/"
task_id <- 1
res <- readRDS(paste0(output.dir, "coeff_task", task_id, ".rds"))
for (task_id in c(2:33)) {
  res0 <- readRDS(paste0(output.dir, "coeff_task", task_id, ".rds"))
  res <- rbind(res, res0)
}

RNA_symbol <- read.csv("../../data/RNA_symbol.csv")
res <- merge(x = RNA_symbol, y = res, by.x = "RNA", by.y = "RNA", all = F)

# remove PAR_Y genes
res <- res %>% filter(!str_detect(RNA, "PAR_Y$"))

res <- res %>%
  rename(Age_Est = Age.enrollment_Est, Age_SE = Age.enrollment_SE, Age_Stat = Age.enrollment_Stat, Age_pval = Age.enrollment_pval) %>%
  rename(Female_Est = SexFemale_Est, Female_SE = SexFemale_SE, Female_Stat = SexFemale_Stat, Female_pval = SexFemale_pval) %>%
  mutate(Age_qval = p.adjust(Age_pval, method = "BH"), .after = "Age_pval") %>%
  mutate(Age_Z = Age_Est / Age_SE, .after = "Age_SE")

saveRDS(res, paste0("full_coeff_", format(Sys.Date(), format = "%y%m%d"), ".rds"))
saveRDS(
  res %>% select(RNA, symbol, starts_with("Age_")),
  paste0("age_coeff_", format(Sys.Date(), format = "%y%m%d"), ".rds")
)
