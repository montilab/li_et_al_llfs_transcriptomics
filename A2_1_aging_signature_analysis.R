## script used to run assoc in bash ##
## use all subjects, including EL ##

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GENESIS))
suppressPackageStartupMessages(library(SeqArray))
suppressPackageStartupMessages(library(SeqVarTools))
suppressPackageStartupMessages(library(Biobase))

#### parameters & logs #### ----------------------------------------------------
root <- ".../batch5_16k/"
output.dir <- paste0(root, "Aging_Signature/v1_frz6_allsubject/result_each_task/")
RNA_file <- paste0(root, "data/RNA_symbol.csv")
input.RNAseq <- paste0(root, "data/expdat_fltrnorm_v1.rds")
input.pheno <- paste0(root, "data/pheno_v1_2167samp_2503.rds")

input.gds <- "..../LLFS.WGS.freeze6.chr21ck04.BI-SNP.gds"
grm <- as.matrix(readRDS(paste0(root, "data/grm_v1_2101samp_2503.rds")))

## split 16k genes into 500 genes per each task
task_id <- as.numeric(Sys.getenv("SGE_TASK_ID"))
RNA_per_task <- 500

sink(paste0(output.dir, "task", task_id, ".log"), append = FALSE, split = TRUE)
date()

RNAall <- read.csv(RNA_file) %>% pull(RNA)
RNAall <- RNAall[((task_id - 1) * RNA_per_task + 1):(task_id * RNA_per_task)]
RNAall <- RNAall[complete.cases(RNAall)]

#### Read Data #### ------------------------------------------------------------
pheno.dat <- readRDS(input.pheno)
RNAseq.dat <- as.data.frame(t(readRDS(input.RNAseq)[RNAall, ]))
RNAseq.dat$subject <- rownames(RNAseq.dat)
analysis.dat <- merge(x = RNAseq.dat, y = pheno.dat, by = "subject", all = F)

analysis.dat <- analysis.dat %>% mutate(
  batch = as.factor(batch),
  fc = relevel(as.factor(fc), ref = "BU"),
  Sex = relevel(as.factor(Sex), ref = "Male")
)

gds <- seqOpen(input.gds)
id.gds <- data.frame(sample.id = seqGetData(gds, "sample.id"))

#### Run Association #### ------------------------------------------------------
annot <- left_join(id.gds, analysis.dat, by = c("sample.id" = "subject"))
seqData <- SeqVarData(gds, sampleData = AnnotatedDataFrame(annot))

res <- c()

## loop through all RNAs
for (idx in c(1:length(RNAall))) {
  RNA0 <- RNAall[idx]
  cat("-- Start", idx, RNA0, "\n")

  nullmod <- try(fitNullModel(seqData,
    outcome = RNA0,
    sample.id = colnames(grm), ## need to subset seqData to subjects available in grm
    covars = c(
      "Age.enrollment", "fc", "Sex", "Education",
      "percent_intergenic", "batch",
      "PC1", "PC2", "PC3", "PC4",
      "PC5", "PC6", "PC7", "PC8",
      "PC9", "PC10"
    ),
    cov.mat = grm,
    family = "gaussian", verbose = T
  ))

  if (class(nullmod) != "try-error") {
    cat("   Sample size:", nrow(nullmod$fit), "\n")
    coeff <- nullmod$fixef

    var_names <- rownames(coeff)
    stat_names <- colnames(coeff)
    new_col_names <- paste(rep(var_names, each = length(stat_names)), stat_names, sep = "_")
    res0 <- data.frame(matrix(unlist(t(coeff)), nrow = 1))
    colnames(res0) <- new_col_names
    res0 <- res0 %>% mutate(RNA = RNA0, .before = everything())
    res <- rbind(res, res0)

    cat("-- Finish", idx, RNA0, "\n\n")
  } else {
    cat("-- Error:", idx, RNA0, "\n\n")
  }
}

#### output #### ---------------------------------------------------------------
saveRDS(res, paste0(output.dir, "coeff_task", task_id, ".rds"))
saveRDS(nullmod, paste0(output.dir, "nullmod_example_task", task_id, ".rds"))

closefn.gds(gds)

date()
sink()
