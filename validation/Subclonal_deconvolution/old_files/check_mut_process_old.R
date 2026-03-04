library(ggplot2)
library(tidyverse)
library(ProCESS)
library(CNAqc)

spn = 'SPN03'
coverage = 100
purity = 0.9
vcf_caller = "mutect2"
cna_caller = "ascat"

github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"

source("/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/generate_table_main.R")
get_cell_id = function(mutation_object) {
  tryCatch(
    expr = { phylo_forest$get_first_occurrences(mutation_object)[[1]] },
    error = function(e) return(NA)
  )
}

mut_process = get_mutations(spn=spn, type="tumour", coverage=coverage, purity=purity)
seq_results = readRDS(mut_process)

# mut_process = seq_results %>%
#   mutate(mutation_id=paste0(spn, ":", chr, ":", chr_pos,  ":", alt),
#          is_driver_process=classes=="driver")

sample_name = 'SPN03_1.1'

sample_forest = load_sample_forest(get_sample_forest(spn))
phylo_forest = load_phylogenetic_forest(get_phylo_forest(spn=spn))
phylo_forest$get_sampled_cell_mutations() %>% view()

# Extract samples MRCA
samples_info <- sample_forest$get_samples_info()
sample_names <- samples_info %>% dplyr::pull(.data$name)


MRCAs_cells <- lapply(sample_names,
                      function(s) {
                        sample_forest$get_coalescent_cells(
                          sample_forest$get_nodes() %>%
                            dplyr::filter(sample %in% s) %>%
                            dplyr::pull(.data$cell_id)
                        ) %>%
                          dplyr::mutate(sample = s)
                      }) %>%
  Reduce(f = dplyr::bind_rows) %>%
  dplyr::group_by(.data$cell_id) %>%
  dplyr::mutate(
    cell_id = paste(.data$cell_id)
  ) %>%
  dplyr::summarise(
    label = paste0("    ", .data$sample, collapse = "\n")
  )

# time_mrca = 0
time_mrca = as.numeric(MRCAs_cells %>% filter(grepl(sample_name, label)) %>% pull(cell_id))

# Time of the sample
sample_time = phylo_forest$get_samples_info() %>% filter(name == sample_name) %>% pull(time)

# Cell id and birth times
sample_cell_ids = phylo_forest$get_nodes() %>% filter((birth_time < sample_time) & (birth_time >= time_mrca))

mut_process_with_clusterid = seq_results %>% 
  filter(classes != "germinal") %>%
  rowwise() %>%
  mutate(cell_id=get_cell_id(Mutation(chr, chr_pos, ref, alt))) %>%
  ungroup() %>% 
  select(cell_id, chr, chr_pos, ref, alt, causes, contains(".VAF")) %>%
  pivot_longer(
    cols=ends_with(".VAF"),
    names_to="sample_id",
    names_pattern="(.*)\\.VAF", # remove matching text "VAF" from the start of each variable name
    values_to="vaf_process" # this is the VAF!
  )

sample_mutations = mut_process_with_clusterid %>% 
  filter(sample_id == sample_name, vaf_process > 0, cell_id %in% sample_cell_ids$cell_id) %>% 
  mutate(from = chr_pos, to = chr_pos, VAF = vaf_process, NV = VAF, DP =VAF)

sample_mutations %>% filter((VAF > vaf_min) & (VAF < vaf_max)) %>% nrow()


n_sample_mutations = nrow(sample_mutations)

# Expected number of mutations ####

# compute length of 1:1 genome

segments = phylo_forest$get_bulk_allelic_fragmentation(sample_name) %>% 
  mutate(from = begin, to = end, Major = major)

x = init(
  mutations = sample_mutations, 
  cna = segments,
  purity = 0.9,
  ref = 'Ghr38'
)

diploid_mutations = x$mutations %>% filter(karyotype == '1:1')

diploid_segments = segments %>% 
  mutate(k = paste0(Major, ':', minor)) %>% 
  filter(k == '1:1') %>% 
  mutate(len = to-from)

l = sum(diploid_segments$len)

b = 0.3
d = 0.01
w = b-d

mu = 1e-9

vaf_min = 0.05
vaf_max = 0.2

1e-9 * l * ((1/vaf_min) - (1/vaf_max))





