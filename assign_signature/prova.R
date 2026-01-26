reticulate::use_python("/orfeo/scratch/area/lvaleriani/myconda/bin")
reticulate::py_config()
library(SigProfilerMatrixGeneratorR)
library(SigProfilerAssignmentR)

muts_process <- readRDS('/orfeo/scratch/cdslab/shared/SCOUT/SPN02/sequencing/tumour/purity_0.9/data/mutations/seq_results_muts_merged_coverage_150x.rds')
base = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature/prova/'

sbs <- c('SBS10a') #'SBS5',

for (s in sbs){
  print(s)
  muts_process <- readRDS('/orfeo/scratch/cdslab/shared/SCOUT/SPN02/sequencing/tumour/purity_0.9/data/mutations/seq_results_muts_merged_coverage_150x.rds')
  
  muts_process <- muts_process %>% 
    dplyr::filter(causes == s)
  
  data <- muts_process %>% 
    ungroup() %>% 
    filter(SPN02_1.1.coverage > 0) %>% 
    select(chr,  chr_pos, ref, alt, causes) %>% 
    rename(cluster = causes) %>% 
    dplyr::mutate(Project = 'spn', Genome = 'GRCh38', mut_type = 'SNP', Type = 'SOMATIC', ID = cluster, Sample = cluster) %>%
    dplyr::rename(chrom = chr, pos_start = chr_pos) %>%
    rowwise() %>%
    dplyr::mutate(pos_end = pos_start + abs(str_count(ref) - str_count(alt))) %>%
    dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type) %>%
    filter(ref != alt) %>% 
    distinct()
  
  out = paste0(base, s)
  dir.create(out, showWarnings = F, recursive = T)
  
  write.table(data, file = paste0(out, '/', s, '.txt'), quote = F, sep = '\t', row.names = F)
  SigProfilerMatrixGeneratorR(project = s, genome = "GRCh38", matrix_path = out, plot=T)
  
  cosmic_fit(samples = paste0(out, '/output/SBS/', s, '.SBS96.all'), 
             output = out, 
             input_type='matrix', 
             context_type='96',
             collapse_to_SBS96=F, 
             cosmic_version=3.3, 
             genome_build="GRCh38", 
             export_probabilities=TRUE,
             make_plots=TRUE)
}
