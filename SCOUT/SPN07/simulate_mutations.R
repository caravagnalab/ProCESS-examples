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
                    passenger_rates = c(SNV = 1e-7, CNA = 1e-12, indel=1e-9),
                    drivers = list(
                      SNV('2', 47799065,'A') # MSH6 c.1082G>A	p.R361H 2-47799065-G-A
                    )
)
m_engine$add_mutant(mutant_name = "Clone 6",
                    passenger_rates = c(SNV = 1e-7, CNA = 1e-12, indel=1e-9),
                    drivers = list(
                      list('TP53 R248W')
                    )
)

## 2. Add exposures
m_engine$add_exposure(c(SBS1 = .35, SBS5 = .4, SBS3 = .25))
m_engine$add_exposure(c(ID1 = .4, ID2=.15, ID4=.25, ID8=0, ID9=.2))

m_engine$add_exposure(time = 96.1,
                      c(SBS11= .9, SBS5=.1)) 
m_engine$add_exposure(time = 96.1,c(ID1 = .1, ID2=.1, ID8=.5, ID9=.3))

m_engine$add_exposure(time = 100, c(SBS1 = .1, SBS5 = .3, SBS3 = .35, SBS26=.12, SBS25=.13))


samples_forest <- load_sample_forest("/orfeo/scratch/cdslab/shared/SCOUT/SPN07/process/sample_forest.sff")
phylo_forest <- m_engine$place_mutations(samples_forest, 500, 200)
phylo_forest$save("/orfeo/scratch/cdslab/shared/SCOUT/SPN07/process/phylo_forest.sff")


#sample_names <- phylo_forest$get_samples_info()[["name"]]
#chromosomes = c(paste0("chr", 1:22), "chrX")
#plot_cna_all_segments = function(cna){
#  cna = cna %>% mutate(chr = paste0("chr", chr))
#  shifted_cna = lapply(chromosomes, function(c){
#    chr_from = CNAqc:::get_reference("hg38") %>% filter(chr == c) %>% pull(from)
#    cna %>% filter(chr == c) %>% mutate(begin_abs = begin + chr_from, end_abs = end + chr_from)
#  })
#  shifted_cna = Reduce(rbind, shifted_cna)
#  CNAqc:::blank_genome() + #ggplot(shifted_cna) +
#    geom_segment(data = shifted_cna %>% mutate(major = ifelse(ratio<1, major+.1,major)), aes(x = begin_abs, xend = end_abs, y = major, yend = major), color = "indianred", size = 1)+
#    geom_segment(data = shifted_cna %>% mutate(minor = ifelse(ratio<1, minor-.1,minor)), aes(x = begin_abs, xend = end_abs, y = minor, yend = minor), color = "steelblue", size = 1)+
#    geom_rect(data = shifted_cna %>% filter(ratio < 1) %>% mutate(a = ifelse(ratio > .5, 1-ratio, ratio)), 
 #             aes(xmin=begin_abs, xmax=end_abs, ymin=0,ymax=Inf, alpha = a), fill="purple")
#}
#plot_cna_clonal = function(cna){
 # cna = cna %>% mutate(chr = paste0("chr", chr))
  #shifted_cna = lapply(chromosomes, function(c){
   # chr_from = CNAqc:::get_reference("hg38") %>% filter(chr == c) %>% pull(from)
    #cna %>% filter(chr == c) %>% mutate(begin_abs = begin + chr_from, end_abs = end + chr_from)
  #})
  #shifted_cna = Reduce(rbind, shifted_cna)
  #CNAqc:::blank_genome() + #ggplot(shifted_cna) +
   # geom_segment(data = shifted_cna %>% filter(ratio==1), aes(x = begin_abs, xend = end_abs, y = major, yend = major), color = "indianred", size = 1)+
   # geom_segment(data = shifted_cna %>% filter(ratio==1), aes(x = begin_abs, xend = end_abs, y = minor, yend = minor), color = "steelblue", size = 1)
#}

#cna_plots = lapply(sample_names,function(s){
 # cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
 # plot_cna(cna) + ggtitle(s)
  # saveRDS(file=paste0("/orfeo/scratch/cdslab/shared/SCOUT/SPN07/process/cna_data/",s,"_cna.rds"),object=cna)
#})

#ggpubr::ggarrange(plotlist = cna_plots, ncol = 2, nrow = 3)

#cna_plots_clonal = lapply(sample_names,function(s){
 # cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
  plot_cna_clonal(cna) + ggtitle(s)
  # saveRDS(file=paste0("/orfeo/scratch/cdslab/shared/SCOUT/SPN07/process/cna_data/",s,"_cna.rds"),object=cna)
#})

#ggpubr::ggarrange(plotlist = cna_plots_clonal, ncol = 2, nrow = 3)

#forest = samples_forest
#annot_forest <- plot_forest(forest) %>%
 # annotate_forest(phylo_forest,
  #                samples = T,
   #               MRCAs = T,
    #              exposures = T,
     #             drivers=T,
      #            add_driver_label = T)

#exp_timeline <- plot_exposure_timeline(phylo_forest)

#labels <- get_relevant_branches(forest)
#sticks <- plot_sticks(forest, labels)

#st = 'AAA
 #     AAA
  #    AAA
   #   AAA
    #  BBB
     # BBB
      #CCC'
#pl = patchwork::wrap_plots(
 # annot_forest , sticks , exp_timeline,
 # design = st
#)
#pl <- annot_forest + sticks + exp_timeline + plot_layout(nrow = 3, design = 'A\nA\nB\nB\nC')
#ggsave("/orfeo/scratch/cdslab/antonelloa/ProCESS-examples/SCOUT/SPN07/mutations.png",plot = pl, dpi = 300, height = 30, width = 15)

