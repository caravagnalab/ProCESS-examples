library(ProCESS)
library(dplyr)

set.seed(9)

p_info <- ps::ps_handle()
start_time <- Sys.time()
initial_cpu <- ps::ps_cpu_times(p_info)
initial_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3

#outdir <- "/scratch/salvatore.milite/SNP05_building_scripts/SPN05/process/"
outdir <- '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/SCOUT/time/'
forest <- load_sample_forest(paste0(outdir,"sample_forest.sff"))


#setwd("/scratch/salvatore.milite/SNP05_building_scripts/")
setwd("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/")
m_engine <- MutationEngine(setup_code = "GRCh38",tumour_type = "BRCA",
                           context_sampling = 20)

mu_SNV = 1e-8
mu_SNV_hyper = 1e-7
mu_CNA = 1e-10
mu_INDELs = 1e-9
mu_INDELs_hyper = 1e-8



CNA_Clone3 = CNA(type = "D", "13",
                 chr_pos = 32315086 - 1e6, len = 1e7,allele = 0)

CNA_Clone4 = CNA(type = "D", "17",
                 chr_pos = 7661779 - 1e6, len = 1e7,allele = 0)

CNA_Clone5 = CNA(type = "D", "9",
                 chr_pos = 21967752 - 1e6, len = 1e7,allele = 0)

## Drivers for the tumors
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA, indel = mu_INDELs),
                    drivers = list(list("TP53 R248W", allele = 1)))

m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(SNV = mu_SNV, CNA = mu_CNA, indel = mu_INDELs),
                    drivers = list(list("BRCA2 E2028*", allele = 1)))

m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV = mu_SNV_hyper, CNA = mu_CNA, indel = mu_INDELs_hyper),
                    drivers = list(CNA_Clone3))

m_engine$add_mutant(mutant_name = "Clone 4",
                    passenger_rates = c(SNV = mu_SNV_hyper, CNA = mu_CNA, indel = mu_INDELs_hyper),
                    drivers = list(CNA_Clone4))


m_engine$add_mutant(mutant_name = "Clone 5",
                    passenger_rates = c(SNV = mu_SNV_hyper, CNA = mu_CNA, indel = mu_INDELs_hyper),
                    drivers = list(CNA_Clone5))

# Mutational signatures
m_engine$add_exposure(
  time = 0 %>% as.numeric(),
  coefficients = c(
    "SBS1" = 0.2,
    "SBS5" = 0.8, "ID1" = 0.6, "ID2" = 0.4)
)

clone3_time <- readRDS(file.path(outdir, "clone3_clock.rds"))

m_engine$add_exposure(
  time = clone3_time %>% as.numeric(),
  coefficients = c(
    SBS1 = 0.1,
    SBS5 = 0.3,
    SBS3 = 0.6, ID7 = 0.7, ID1 = 0.2, ID2= 0.1)
)



m_engine$set_germline_subject("NA20514")

print("Mutation engine created")
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save(paste0(outdir,"phylo_forest.sff"))
# dir.create(paste0(outdir,"cna_data"))
# sample_names <- phylo_forest$get_samples_info()[["name"]]
# lapply(sample_names,function(s){
#   cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
#   saveRDS(file=paste0(outdir,"cna_data/",s,"_cna.rds"),object=cna)
# })
# 


end_time <- Sys.time()
final_cpu <- ps::ps_cpu_times(p_info)
final_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3
elapsed_time <- end_time - start_time
elapsed_time <- as.numeric(elapsed_time, units = "mins")
cpu_used <- (final_cpu["user"] + final_cpu["system"]) - (initial_cpu["user"]+ initial_cpu["system"])
mem_used <- final_mem - initial_mem
resource_usage <- data.frame(
  elapsed_time_mins =  elapsed_time,
  cpu_time_secs = cpu_used,
  memory_used_MB = mem_used
)
saveRDS(object = resource_usage, file = file.path(outdir, 'SPN05_tissue.rds'))

