library(ProCESS)
library(ggplot2)
library(dplyr)
library(optparse)
library(tidyr)

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN04'),
                    make_option(c("--purity"), type = "character", default = '0.9'),
                    make_option(c("--coverage"), type = "character", default = '50')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

spn <- opt$spn_id #"SPN01"
cov <- as.numeric(opt$coverage)
pur <- as.numeric(opt$purity)
combination <- paste0(cov,"x_",pur,"p")
# source("compute_FGA.R")

outidr <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/mutations_counts/"
bases <- c("A","T","G","C")


message(paste0("Processing ",spn))
message(paste0("Processing combination ",combination))
process_muts <- readRDS(
  get_mutations(spn = spn, coverage = cov, purity = pur, type = "tumour")
) %>%
  dplyr::filter(classes != "germinal") %>%
  mutate(
    type = case_when(
      (ref %in% bases & alt %in% bases) ~ "SNV",
      TRUE ~ "INDEL"
    )
  ) %>%
  select(type, contains("occurrences")) %>%
  group_by(type) %>%
  summarise(across(
    contains("occurrences"),
    ~ sum(.x > 0),
    .names = "{.col}"
  ), .groups = "drop") %>%
  rename_with(~ gsub("\\.occurrences$", "", .), ends_with(".occurrences")) %>%
  pivot_longer(
    cols = starts_with("SPN"),   # adjust if your sample columns differ
    names_to = "sample",
    values_to = "mutation_count"
  ) %>% 
  dplyr::mutate(coverage=cov,
                purity=pur)
dir.create(path = paste0(outidr,spn),recursive = T)
saveRDS(process_muts,file=paste0(outidr,spn,"/mutations_counts_",combination,".rds"))
