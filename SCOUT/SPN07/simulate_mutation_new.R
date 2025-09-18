rm(list = ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
set.seed(06117)

setwd("/orfeo/scratch/cdslab/shared/SCOUT/")

m_engine <- MutationEngine(setup_code = "GRCh38",
                           tumour_type = "GBM",
                           #tumour_study = "US",
                           context_sampling = 20)
m_engine$set_germline_subject("NA18941")

## 1. Add Drivers
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = 1e-8, CNA = 1e-9, indel = 1e-9),
                    drivers = list(list("PTEN R130G", allele=0))
) # list(list('PTEN R130*', allele=0)

m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(SNV = 1e-8, CNA = 1e-12, indel=1e-9),
                    drivers = list(
                      CNA(type="D", chr="10", chr_pos=41545820+1e6, len=91251602 ,allele = 1)
                    )
)

m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV = 1e-8, CNA = 1e-12, indel=1e-9),
                    drivers = list(
                      list('NF1 R192*', allele=0)
                    )
)

# Convert mutation in protein coordinates to genome coordinates : https://bibliome.ai/GRCh38/gene/ATRX
# m_engine$add_mutant(mutant_name = "Clone 4",
#                     passenger_rates = c(SNV = 1e-8, CNA = 1e-12, indel=1e-9),
#                     drivers = list(
#                       SNV('X', 77682537, 'A') # ARTX R907*
#                     )
# )

m_engine$add_mutant(mutant_name = "Clone 4",
                    passenger_rates = c(SNV = 1e-8, CNA = 1e-12, indel=1e-9),
                    drivers = list(
                      list('ATRX R907*', allele=0)
                    )
)

m_engine$add_mutant(mutant_name = "Clone 5",
                    passenger_rates = c(SNV = 1e-7, CNA = 1e-12, indel=1e-8),
                    drivers = list(
                      SNV('2', 47799065,'A') # MSH6 c.1082G>A   p.R361H 2-47799065-G-A
                    )
)
m_engine$add_mutant(mutant_name = "Clone 6",
                    passenger_rates = c(SNV = 1e-7, CNA = 1e-12, indel=1e-8),
                    drivers = list(
                      list('TP53 R248W')
                    )
)

## 2. Add exposures
m_engine$add_exposure(c(SBS1 = .2, SBS5 = .45, SBS3 = .35))
m_engine$add_exposure(c(ID1 = .4, ID2=.15, ID4=.25, ID8=0, ID9=.2))
m_engine$add_exposure(time = 158.128336812493, c(ID1 = .1, ID2=.1, ID8=.5, ID9=.3))
m_engine$add_exposure(time = 158.128336812493, c(SBS11= .9, SBS5=.1))
m_engine$add_exposure(time = 173, c(SBS1 = .1, SBS5 = .3, SBS3 = .35, SBS26=.12, SBS25=.13))


samples_forest <- load_sample_forest("/orfeo/scratch/cdslab/shared/SCOUT/SPN07_last/process/sample_forest.sff")
phylo_forest <- m_engine$place_mutations(samples_forest, 
                                         num_of_preneoplatic_SNVs = 800,
                                         num_of_preneoplatic_indels = 200)
phylo_forest$save("/orfeo/scratch/cdslab/shared/SCOUT/SPN07_last/process/phylo_forest.sff")

