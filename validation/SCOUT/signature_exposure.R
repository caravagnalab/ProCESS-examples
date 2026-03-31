library(ProCESS)
library(ggplot2)
library(dplyr)
library(optparse)
library(tidyr)

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/nf-validation/bin/somatic/utils/plot_utils.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/nf-validation/bin/signatures/utils_getters.R")
option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN04'),
                                make_option(c("--purity"), type = "character", default = '0.9'),
                                make_option(c("--coverage"), type = "character", default = '50'))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

spn <- opt$spn_id #"SPN01"
coverage <- as.numeric(opt$coverage)
purity <- as.numeric(opt$purity)
combination <- paste0(coverage,"x_",purity,"p")

df_all_SPN_exposure_sbs <- list()
df_all_SPN_exposure_id <- list()

message(paste0("Processing ",spn))
validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION"

validation_dir_somatic <- file.path(validation_dir,spn,"signature")

# coverage <- 100
# purity <- 0.9
process_exposures <- get_process_exposures(spn = spn,coverage = coverage,purity = purity)
process_exposures_sbs_m <- process_exposures$SBS %>% 
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  select(!c("Sample_ID")) %>% 
  as.matrix() %>% t()
m_sbs = cbind(process_exposures_sbs_m,process_exposures_sbs_m)
process_exposures_id_m <- process_exposures$ID %>% 
  mutate(across(everything(), ~replace_na(., 0))) %>% 
  select(!c("Sample_ID")) %>% 
  as.matrix() %>% t()
m_id = cbind(process_exposures_id_m,process_exposures_id_m)
outdir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/signature_exposures/"
dir.create(path = paste0(outdir,spn),recursive = T)
saveRDS(object = process_exposures_sbs_m,file = paste0(outdir,spn,"/signature_exposure_",combination,"_SBS.rds"))
saveRDS(object = process_exposures_id_m,file = paste0(outdir,spn,"/signature_exposure_",combination,"_ID.rds"))

