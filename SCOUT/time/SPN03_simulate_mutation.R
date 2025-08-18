library(ProCESS)
library(dplyr)

#base="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN03/process/"
base ="/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/SCOUT/test/"
setwd(base)
set.seed(1806)

p_info <- ps::ps_handle()
start_time <- Sys.time()
initial_cpu <- ps::ps_cpu_times(p_info)
initial_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3

forest <- load_sample_forest("sample_forest.sff")

setwd('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/')
m_engine <- MutationEngine(setup_code = "GRCh38", 
                           tumour_type = 'CLLSLL', 
                           context_sampling = 20)

SNV_rate = 1e-9
indel_rate = 1e-10
CNA_rate = 1e-12

# Clone 1 
# NOTCH1p2514*fs*4  
# 13q14.2  deletion
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = SNV_rate,
                                        CNA = CNA_rate,
                                        indel = indel_rate),
                    drivers = list("NOTCH1 P2514Rfs*4", 
                                   CNA(chr = "13",
                                       chr_pos = 45000000, 
                                       len = 1e7,  
                                       type = "D"))
)

# Clone 2 
# KRAS G12D
m_engine$add_mutant(mutant_name = "Clone 2",
                    passenger_rates = c(SNV = SNV_rate,
                                        CNA = CNA_rate,
                                        indel = indel_rate),
                    driver = list("TP53 R175H")
)

# Clone 3
# Unknown
# MAP2K1 P124L
# 15	66436825	66436825	C	G
# https://pmc.ncbi.nlm.nih.gov/articles/PMC4815041/
m_engine$add_mutant(mutant_name = "Clone 3",
                    passenger_rates = c(SNV = SNV_rate,
                                        CNA = CNA_rate,
                                        indel = indel_rate),
                    driver = list(SNV("15", 66436825, alt = "G", ref = 'C'))
)

# Signatures
# SBS1, SBS5
m_engine$add_exposure(c(SBS5 = 0.3, SBS1 = 0.2, SBS9 = 0.5,
                        ID5 = 0.6, ID1 = 0.2, ID2 = 0.2))

phylo_forest <- m_engine$place_mutations(forest, 
                                         num_of_preneoplatic_SNVs = 800,
                                         num_of_preneoplatic_indels = 200)

phylo_forest$save(paste0(base,"phylo_forest.sff"))


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
saveRDS(object = resource_usage, file = 'SPN03_mutations.rds')

# dir.create(paste0(base, "cna_data"),recursive = T)
# sample_names <- phylo_forest$get_samples_info()[["name"]]
# lapply(sample_names,function(s){
#   cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
#   saveRDS(file=paste0(base, "cna_data/",s,"_cna.rds"), object=cna)
# })

