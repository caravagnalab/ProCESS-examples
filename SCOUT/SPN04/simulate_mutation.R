rm(list = ls())
library(ProCESS)
library(dplyr)
seed <- 12345
set.seed(seed)

outdir <- "/orfeo/scratch/cdslab/shared/SCOUT/SPN04/process/"

forest <- load_sample_forest(paste0(outdir,"sample_forest.sff"))
treatment_info <- readRDS(paste0(outdir,"treatment_info.rds"))

setwd("/orfeo/scratch/cdslab/shared/SCOUT/")
#setwd(outdir)
# Set tumour type and gender

m_engine <- MutationEngine(setup_code = "GRCh38", tumour_type = "AML", context_sampling = 20)
m_engine$set_germline_subject("NA18940")

mu_SNV = 1e-9
mu_ID = 1e-10
mu_CNA = 1e-12

driver_Clone1 = list("IDH1 R132C")
driver_Clone2 = CNA(type='A', chr='8', chr_pos=46000001, len=99138636)
driver_Clone3 = list("NRAS Q61K")

# Add drivers ####

m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA, indel=mu_ID),
                    drivers = list(driver_Clone1))

m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA, indel=mu_ID),
                    drivers = list(driver_Clone2))

m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV=mu_SNV, CNA=mu_CNA, indel=mu_ID),
                    drivers = list(driver_Clone3))


# Mutational Signatures ####
m_engine$add_exposure(time = 0, coefficients = c(SBS5 = 0.5, SBS1 = 0.2, SBS18 = 0.3, ID1 = 0.5, ID2 = 0.5))
m_engine$add_exposure(time = treatment_info$treatment_start, coefficients = c(SBS25 = 1, ID1 = 0.5, ID2 = 0.5))
m_engine$add_exposure(time = treatment_info$treatment_end, coefficients = c(SBS5 = 0.5, SBS1 = 0.2, SBS18 = 0.3, ID1 = 0.5, ID2 = 0.5))

# Phylo forest ######
print("Mutation engine created")
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs=800, num_of_preneoplatic_indels=200)
phylo_forest$save(paste0(outdir,"phylo_forest.sff"))

print("Saving cna")
dir.create(paste0(outdir,"cna_data"))
sample_names <- phylo_forest$get_samples_info()[["name"]]
print(sample_names)

for (s in sample_names) {
  cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
  print(paste0(outdir,"cna_data/",s,"_cna.rds"))
  saveRDS(file=paste0(outdir,"cna_data/",s,"_cna.rds"),object=cna)
}

