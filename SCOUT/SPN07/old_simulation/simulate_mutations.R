rm(list = ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
set.seed(06117)

setwd("/orfeo/cephfs/scratch/cdslab/shared/ProCESS/GRCh38")

m_engine <- MutationEngine(setup_code = "GRCh38",
                           tumour_type = "GBM",
                           tumour_study = "US")
m_engine$set_germline_subject("NA18941")

## 1. Add Drivers
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = 1e-8, CNA = 1e-9, indel = 1e-9),
                    drivers = list(list('PTEN R130*', allele=0))
)
m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(SNV = 1e-8, CNA = 1e-12, indel=1e-9),
                    drivers = list(
                      CNA("D", "10", 87862638 , 109292 ,allele = 1)
                    )
)
m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV = 1e-8, CNA = 1e-12, indel=1e-9),
                    drivers = list(
                      list('NF1 Q1174*', allele=0)
                    )
)
# Convert mutation in protein coordinates to genome coordinates : https://bibliome.ai/GRCh38/gene/ATRX
m_engine$add_mutant(mutant_name = "Clone 4",
                    passenger_rates = c(SNV = 1e-8, CNA = 1e-12, indel=1e-9),
                    drivers = list(
                      SNV('X', 77682537, 'A') # ARTX R907*
                    )
)
m_engine$add_mutant(mutant_name = "Clone 5",
                    passenger_rates = c(SNV = 1e-7, CNA = 1e-12, indel=1e-9),
                    drivers = list(
                      SNV('2', 47799065,'A') # MSH6 c.1082G>A	p.R361H
                    )
)
m_engine$add_mutant(mutant_name = "Clone 6",
                    passenger_rates = c(SNV = 1e-7, CNA = 1e-12, indel=1e-9),
                    drivers = list(
                      list('TP53 R248W')
                    )
)

## 2. Add exposures
m_engine$add_exposure(c(SBS1 = .8, SBS5 = .2))
m_engine$add_exposure(c(ID1 = 1))
m_engine$add_exposure(time = 108.4,
                      c(SBS11= 1)) # Chemotherapy active from 7.85 to 10.85
m_engine$add_exposure(time = 110.4, c(ID1 = 1))
# Hypermutant signatures : SBS6, SBS14, SBS15, SBS20, SBS21, and SBS44.
m_engine$add_exposure(time = 110.4, c(SBS1 = .1, SBS5 = .1,SBS6=.4,
                                      SBS21=.2, SBS44=.2))


samples_forest <- load_samples_forest("/orfeo/scratch/cdslab/shared/SCOUT/SPN07/process/sample_forest.sff")
phylo_forest <- m_engine$place_mutations(samples_forest, 500, 200)
phylo_forest$save("/orfeo/scratch/cdslab/shared/SCOUT/SPN07/process/phylo_forest.sff")


sample_names <- phylo_forest$get_samples_info()[["name"]]
lapply(sample_names,function(s){
  cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
  saveRDS(file=paste0("/orfeo/scratch/cdslab/shared/SCOUT/SPN07/process/cna_data/",s,"_cna.rds"),object=cna)
})

forest = samples_forest
annot_forest <- plot_forest(forest) %>%
  annotate_forest(phylo_forest,
                  samples = T,
                  MRCAs = T,
                  exposures = T,
                  drivers=T,
                  add_driver_label = T)

exp_timeline <- plot_exposure_timeline(phylo_forest)

labels <- get_relevant_branches(forest)
sticks <- plot_sticks(forest, labels)

st = 'AAA
      AAA
      AAA
      AAA
      BBB
      BBB
      CCC'
pl = patchwork::wrap_plots(
  annot_forest , sticks , exp_timeline,
  design = st
)
#pl <- annot_forest + sticks + exp_timeline + plot_layout(nrow = 3, design = 'A\nA\nB\nB\nC')
ggsave("/orfeo/scratch/cdslab/antonelloa/ProCESS-examples/SCOUT/SPN07/final_simulations/mutations.png",plot = pl, dpi = 300, height = 30, width = 15)

