# adding mutations
setwd('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN02/process')
rm(list=ls())

# Mutation generation ####
library(ProCESS)
library(dplyr)

# load time in which is born the third clone and the samples forest
clone3_born = readRDS("clone3_clock.rds")
sampled_phylogeny = load_sample_forest("sample_forest.sff")

curr_dir = getwd()
setwd('/orfeo/cephfs/scratch/cdslab/shared/SCOUT')

# setting the germline subject as columbian woman and tumor type to colonrectal cancer
mu_SNV = 1e-8
mu_ID = 1e-9
mu_CNA = 1e-11
 
# last clone is hypermutant (only for snvs and ids, lower CNAs)
mu_SNV_clone3 = 1e-7
mu_ID_clone3 = 7e-9
mu_CNA_clone3 = 1e-11

m_engine <- MutationEngine(setup_code = "GRCh38", 
                           context_sampling = 20,
                           tumour_type = "COADREAD",
                           germline_subject = "HG01113")

m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = mu_SNV, indel = mu_ID, CNA = mu_CNA),
                    drivers = list("BRAF V600E"))

m_engine$add_mutant(mutant_name = "Clone 2", 
                    passenger_rates = c(SNV = mu_SNV, indel = mu_ID, CNA = mu_CNA),
                    drivers = list("PIK3CA E545K"))

m_engine$add_mutant(mutant_name = "Clone 3", 
                    passenger_rates = c(SNV = mu_SNV_clone3, indel = mu_ID_clone3, CNA = mu_CNA_clone3), # modelled as an hypermutant, so with higher passenger_rates
                    # drivers = list(SNV("2", 47799065, "A")))
                    drivers = list(SNV("12", 132673703, "A"))) # POLE somatic mutation POLE p.V411L 
                                             # G>C
# indels and sbs have different coefficients summing (up to 1 for IDs and up to 1 for SBSs)
m_engine$add_exposure(
    time = 0, 
    coefficients = c(
        ID1 = 0.7,
        ID2 = 0.3,
        SBS1 = 0.2,
        SBS5 = 0.8)
)

m_engine$add_exposure(
    time = round(clone3_born),
    coefficients = c(
        SBS1 = 0.2, SBS5 = 0.2, SBS10b = 0.3, SBS10a = 0.3, # hypermutant POLE signature
        ID1 = 0.2, ID7 = 0.3, ID2 = 0.5)
)

# Sys.time()
# add mutations on the forest with 1000 pre-neoplastic mutations
phylo_forest <- m_engine$place_mutations(sampled_phylogeny, 
    num_of_preneoplatic_SNVs = 1000, 
    num_of_preneoplatic_indels = 200
)

print("Mutations placed")
# save the phylogenetic forest in the file "phylo_forest.sff"
phylo_forest$save(paste(curr_dir, "phylo_forest.sff", sep = "/"))

print("Forest saved")

Sys.time()

print('saving cna')

dir.create(paste0(curr_dir,"/cna_data"))
sample_names <- phylo_forest$get_samples_info()[["name"]]
lapply(sample_names,function(s) {

 print(Sys.time())	       
 print(paste0('saving cnas sample ', s))

 cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
 print(Sys.time())
 saveRDS(file=paste0(curr_dir,"/cna_data/",s,"_cna.rds"),object=cna)
 print(paste0('cnas saved sample', s))
})

# # plotting
# tree_plot = plot_forest(sampled_phylogeny) 
# mut_forest = annotate_forest(tree_plot,
#                 phylo_forest,
#                 samples = TRUE,
#                 MRCAs = TRUE,
#                 exposures = T,
#                 facet_signatures = F,
#                 drivers = T,
#                 add_driver_label = T) +
#             ggplot2::facet_wrap( ~ signature, scales = "free")
# 
# ggsave(filename = "/orfeo/cephfs/scratch/area/vgazziero/CDSlab/rRaces/rRACES-examples/SPN02/phylogenetic_forest_cycles.png", plot = mut_forest, width = 20, height = 15, bg = "white")
# print("done everything")
# 
# exposure_time = plot_exposure_timeline(
#   phylo_forest,
#   linewidth = 0.8,
#   emphatize_switches = TRUE
# ) 
# ggsave(filename = "/orfeo/scratch/area/vgazziero/CDSlab/rRaces/rRACES-examples/SPN02/exposure_plot.png", plot = exposure_time, bg = "white")
# 
# all = mut_forest / exposure_time
# ggsave(filename = "/orfeo/scratch/area/vgazziero/CDSlab/rRaces/rRACES-examples/SPN02/recap_phylo.png", plot = all)
