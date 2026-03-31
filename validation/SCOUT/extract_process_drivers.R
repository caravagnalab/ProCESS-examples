library(ProCESS)
library(ggplot2)
library(dplyr)
library(optparse)
library(tidyr)

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN04'))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

spn <- opt$spn_id #"SPN01"

outidr <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/drivers/"


phylo_forest <- load_phylogenetic_forest(get_phylo_forest(spn = spn))
drivers_snv <- phylo_forest$get_driver_mutations() %>% 
  dplyr::mutate(SPN=spn)

dir.create(path = paste0(outidr,spn),recursive = T)
saveRDS(drivers_snv,file=paste0(outidr,spn,"/process_drivers.rds"))

############### JUST AFTER RUNNING THE UPPER PART ######################

#spns <- c("SPN01","SPN03","SPN04","SPN06","SPN07")
#driver_table_all <- lapply(spns, function(x){
#  dr_t <- readRDS(file = file.path(outidr,x,"process_drivers.rds"))
#}) %>% bind_rows() %>% 
#  mutate(code=case_when(SPN=="SPN07" & type=="SID" & is.na(code) ~ "MSH6 p.R361H",
#                        TRUE ~ code))
#saveRDS(object = driver_table_all,file = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/drivers/all_drivers.rds")
