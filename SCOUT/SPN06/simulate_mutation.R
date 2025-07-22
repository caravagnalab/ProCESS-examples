

rm(list = ls())
library(ProCESS)
library(dplyr)
seed <- 19999
set.seed(seed)

outdir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN06/process/"

forest <- load_samples_forest(paste0(outdir, "sample_forest.sff"))

timing <- readRDS(paste0(outdir, "chemo_timing.rds"))

#-------------------------------------------------------------------------------
#-------------------------- set up Mutation Engine -----------------------------
#-------------------------------------------------------------------------------
setwd("/orfeo/cephfs/scratch/cdslab/shared/SCOUT")
m_engine <- MutationEngine(setup_code = "GRCh38", tumour_type = "LUAD", context_sampling = 20)


#-------------------------------------------------------------------------------
#---------------------------- passenger mutations ------------------------------
#-------------------------------------------------------------------------------
passengers <- c(SNV = 1e-8, CNA = 1e-12, indel = 1e-9)

#-------------------------------------------------------------------------------
#------------------------------ driver mutations -------------------------------
#-------------------------------------------------------------------------------

# TP53 R248Q/W/L
#-------------------------------------------------------------------------------
SNV_C1 = "TP53 R248W"


# STK11 LOH
#-------------------------------------------------------------------------------
CNA_C2 <- CNA(
  type = "D", 
  chr = "19", 
  chr_pos = 702994, 
  len = 2e7, 
  allele = 0
)


# EGFR amp
#-------------------------------------------------------------------------------
CNA_C3 <- CNA(
  type = "A", 
  chr = "7", 
  chr_pos = 54615322,  
  len = 2e7
)


# KEAP1 R413C/H/L
#-------------------------------------------------------------------------------
#SNV_C4 <- "KEAP1 R460M"
SNV_C4 <- "KEAP1 R413H"


# KRAS G12D
#-------------------------------------------------------------------------------
SNV_C5 <- "KRAS G12D"


# KRAS amp (mutant)
#-------------------------------------------------------------------------------
CNA_C6 <- CNA(
  type = "A", 
  chr = "12", 
  chr_pos = 24728091, 
  len = 2e7
)

# whole genome doubling
#-------------------------------------------------------------------------------
#CNA_C7 <- WGD


#-------------------------------------------------------------------------------
#-------------------------- add mutations to mutants ---------------------------
#-------------------------------------------------------------------------------

# Clone 1 : TP53 R248Q/W/L
#-------------------------------------------------------------------------------
m_engine$add_mutant(
  mutant_name = "Clone 1", 
  passenger_rates = passengers, 
  drivers = list(SNV_C1)
)

# Clone 2 : STK11 LOH
#-------------------------------------------------------------------------------
m_engine$add_mutant(
  mutant_name = "Clone 2", 
  passenger_rates = passengers, 
  drivers = list(CNA_C2)
)

# Clone 3 : EGFR amp
#-------------------------------------------------------------------------------
m_engine$add_mutant(
  mutant_name = "Clone 3", 
  passenger_rates = passengers, 
  drivers = list(CNA_C3)
)

# Clone 4 : KEAP1 R413C/H/L
#-------------------------------------------------------------------------------
m_engine$add_mutant(
  mutant_name = "Clone 4", 
  passenger_rates = passengers, 
  drivers = list(SNV_C4)
)

# Clone 5 : KRAS G12D
#-------------------------------------------------------------------------------
m_engine$add_mutant(
  mutant_name = "Clone 5", 
  passenger_rates = passengers, 
  drivers = list(SNV_C5)
)

# Clone 6 : KRAS amp (mutant)
#-------------------------------------------------------------------------------
m_engine$add_mutant(
  mutant_name = "Clone 6", 
  passenger_rates = passengers, 
  drivers = list(CNA_C6)
)

# Clone 7 : whole genome doubling
#-------------------------------------------------------------------------------
m_engine$add_mutant(
  mutant_name = "Clone 7", 
  passenger_rates = passengers, 
  drivers = list(WGD)
)


#-------------------------------------------------------------------------------
#----------------------------- mutational signatures ---------------------------
#-------------------------------------------------------------------------------
# SBS1, SBS4, SBS5, ID1 and ID2 are always active
m_engine$add_exposure(time = 0, coefficients = c(SBS1 = 0.2, SBS4 = 0.5, SBS5 = 0.3, ID1 = 0.6, ID2 = 0.4))

# SBS11 active after timolozomide
m_engine$add_exposure(time = timing$chemo1_start, coefficients = c(SBS1 = 0.15, SBS4 = 0.45, SBS5 = 0.25, SBS11 = 0.15, ID1 = 0.7, ID2 = 0.3))

# SBS25 active during chemotherapy
m_engine$add_exposure(time = timing$chemo2_start, coefficients = c(SBS1 = 0.15, SBS4 = 0.45, SBS5 = 0.25, SBS11 = 0.1, SBS25 = 0.05, ID1 = 0.7, ID2 = 0.3))

# SBS25 de-active after chemotherapy
m_engine$add_exposure(time = timing$chemo2_end, coefficients = c(SBS1 = 0.15, SBS4 = 0.45, SBS5 = 0.25, SBS11 = 0.15, ID1 = 0.8, ID2 = 0.2))


#-------------------------------------------------------------------------------
#------------------------------ place mutations --------------------------------
#-------------------------------------------------------------------------------
phylo_forest <- m_engine$place_mutations(forest, num_of_preneoplatic_SNVs = 800, num_of_preneoplatic_indels = 200)

phylo_forest$save(paste0(outdir, "phylo_forest.sff"))

dir.create(paste0(outdir,"cna_data"), recursive = TRUE)

sample_names <- phylo_forest$get_samples_info()[["name"]]
lapply(sample_names, function(s) {
  cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
  saveRDS(file = paste0(outdir, "cna_data/", s, "_cna.rds"), object=cna)
})

message("ALL DONE!")



