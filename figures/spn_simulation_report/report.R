library(ggplot2)
library(patchwork)
library(ProCESS)
library(dplyr)
library(tidyverse)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/report/")
source("plotting/signature_ProCESS.R")
source("plotting/plot_genome_wide.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
source("plotting/dynamics_ProCESS.R",local =T)
source("plotting/tables.R", local = knitr::knit_global())
metadata <- read.table(file = "SCOUT_metadata.csv",header = T,sep = "\t")

spn <- "SPN04"
sample_forest <- load_sample_forest(get_sample_forest(spn = spn))
phylo_forest <- load_phylogenetic_forest(get_phylo_forest(spn = spn))

sample_names <- get_sample_names(spn = spn) 
number_of_samples <- length(sample_names)
cna_data <- lapply(sample_names,function(s){
  readRDS(get_process_cna(spn = spn,sample = s))
}) %>% bind_rows()


info_spn <- metadata %>% filter(SPN_ID==spn)
nrow = ceiling(number_of_samples/3)
height = nrow * 3

basedir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
spndir <- file.path(basedir,spn,"process")
sim <- recover_simulation(spn)
setwd(spndir)
plot <- plot_tumour_dynamics(spn, sample_forest)
plot$plot_dynamics
plot$plot_sampling

color_map_clones <-  get_clone_map(sample_forest)



#sample_forest <- load_samples_forest(forest)
info = sim$get_samples_info() ## requested from either the simulation recovery or as saved table

nodes = sample_forest$get_nodes()
clones = nodes %>% 
  dplyr::filter(!is.na(sample)) %>% 
  dplyr::group_by(sample, mutant) %>% 
  dplyr::pull(mutant) %>% 
  unique()
clones_prop = nodes %>%
  dplyr::filter(!is.na(sample)) %>% 
  dplyr::group_by(sample, mutant) %>% 
  # dplyr::mutate(mutant = gsub(" ", "_", mutant)) %>% 
  dplyr::count(mutant) %>% 
  group_by(sample) %>%
  mutate(proportion = n / sum(n))

# plot pie charts
sample_composition <- ggplot(clones_prop, aes(x = "", y = proportion, fill = mutant)) +
  geom_col(width = 1, color = "white") +
  scale_fill_manual(values=color_map_clones)+
  coord_polar(theta = "y") +
  facet_wrap(~ sample) +
  theme_void()

plt_sample <- plot_forest(sample_forest,
                          color_map = color_map_clones) +
  theme(legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing = unit(0.01, "cm"))
plt_sample_annotated <- annotate_forest(plt_sample, phylo_forest, samples = TRUE, MRCAs = TRUE, driver = F)

exposure_table <- phylo_forest$get_exposures()
df_sign <- get_exposure_ends(phylo_forest) %>% 
  mutate(start_time = round(time),
         end_time = round(end_time)-1) %>% 
  select(start_time, end_time, signature, exposure, type)

plot_signatures <- plot_exposure_evolution(sample_forest, phylo_forest, file.path(spndir,spn))
plot_signatures$sign_tree
plot_signatures$sign_muller

setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/spn_simulation_report")
ggsave(filename = paste0(spn,"_sampling.png"),plot = plot$plot_sampling,width = number_of_samples*3,height = 3,units = "in")
ggsave(filename = paste0(spn,"_dynamics.pdf"),plot = plot$plot_dynamics,width = 8,height = 8,units = "in")
ggsave(filename = paste0(spn,"_phylo.png"),plot = plt_sample_annotated,width = 4,height = 6,units = "in")

ggsave(filename = paste0(spn,"_exp_muller.pdf"),plot = plot_signatures$sign_muller,width = 4,height = 4,units = "in")
ggsave(filename = paste0(spn,"_sample_composition.pdf"),plot = sample_composition,width = number_of_samples*3,height = 3,units = "in")
#final_plt <- wrap_plots(plot$plot_sampling,plot$plot_dynamics,sample_composition,plot_signatures$sign_muller,plt_sample_annotated)
