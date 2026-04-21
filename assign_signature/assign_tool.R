setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")
library(ProCESS)
library(tidyverse)
library(optparse)
source('../getters/tumourevo_getters.R')
source('../getters/process_getters.R')

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN04'),
                    make_option(c("--purity"), type = "double", default = 0.9),
                    make_option(c("--coverage"), type = "integer", default = 50),
                    make_option(c("--cna_caller"), type = "character", default = 'ascat'),
                    make_option(c("--vcf_caller"), type = "character", default = 'mutect2'),
                    make_option(c("--signature"), type = "character", default = 'BASCULE'),
                    make_option(c("--new"), type = "logical", default = TRUE)
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

signature_tool <- opt$signature

if (signature_tool == 'SigProfiler'){
  reticulate::use_python("/orfeo/scratch/area/lvaleriani/myconda/bin")
  reticulate::py_config()
  library(SigProfilerMatrixGeneratorR)
  library(SigProfilerAssignmentR)
} else if (signature_tool == 'BASCULE'){
  reticulate::use_condaenv("/orfeo/scratch/cdslab/ggandolfi/miniconda/envs/bascule-env")
  py = reticulate::import("pybascule")
  library(bascule)
}

base = '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/validation_subclonal_new/tables_interpreted/'
spn = opt$spn_id
cov = opt$coverage
pur = opt$purity

cna_caller = opt$cna_caller
mut_caller = opt$vcf_caller

out = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assign_signature/"
for (tool in c('pyclonevi', 'viber')){
  
  out_data_raw = paste0(out, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', signature_tool, '_raw/')
  out_data_int = paste0(out, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', signature_tool, '_int/')
  dir.create(out_data_raw, recursive = T, showWarnings = F)
  dir.create(out_data_int, recursive = T, showWarnings = F)
  
  if (opt$new & signature_tool == 'SigProfiler'){
    print('Deleting old files')
    unlink(paste0(out_data_raw, '/*'), recursive = T)
    unlink(paste0(out_data_int, '/*'), recursive = T)
  }
  
  
  if (tool == 'viber'){
    tool_table <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/", spn, "/tumourevo/", cov, "x_",pur,"p_",mut_caller,"_",cna_caller,"/subclonal_deconvolution/viber/SCOUT/",spn, "/SCOUT_", spn,"_viber_best_fit.rds"))
    mut_ids <- tool_table$data %>% 
      select(chr, from, ref, alt) %>% 
      distinct() %>% 
      separate(col = chr, sep = 'chr', into = c('tmp', 'chr')) %>% 
      mutate(mutation_id = paste(spn, chr, from, alt, sep = ':')) %>% 
      select(-tmp)
    
  } else if (tool == 'pyclonevi'){
    tool_table <- get_tumourevo_subclonal(spn = spn, 
                                          coverage = cov, 
                                          purity = pur, 
                                          tool = tool,
                                          vcf_caller = mut_caller,
                                          cna_caller = cna_caller)
    
    fit <- read.table(tool_table$best_fit_txt, header = T, sep = '\t') %>% 
      select(mutation_id, cluster_id) %>% 
      distinct()
    
    sample = get_sample_names(spn)[[1]]
    vcf <- readRDS(get_tumourevo_driver(spn = spn, 
                                        coverage = cov, 
                                        purity = pur, 
                                        vcf_caller = mut_caller,
                                        cna_caller = cna_caller,
                                        sample = sample))
    
    mut_ids <- vcf[[paste0(spn, '_', sample)]][['mutations']] %>% 
      mutate(mutation_id = paste(spn, str_replace(chr, patter = 'chr', ''), from, alt, sep = ':')) %>% 
      select(chr, from, to, ref, alt, mutation_id) %>% 
      distinct() %>% 
      separate(col = chr, sep = 'chr', into = c('tmp', 'chr'))  %>% 
      select(-tmp)
  }
    
    table <- readRDS(paste0(base, tool, '_', spn, '_', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '.rds'))
    raw <- table %>% select(patient_id, mutation_id, cluster_id_tool)
    int <- table %>% select(patient_id, mutation_id, cluster_id_tool_interpreted)

    
    if (signature_tool == 'SigProfiler'){ 
      table_raw <- raw %>% 
        left_join(mut_ids) %>% 
        select(chr,  from, ref, alt, cluster_id_tool) %>% 
        dplyr::rename(cluster = cluster_id_tool) %>% 
        dplyr::mutate(Project = spn, Genome = 'GRCh38', mut_type = 'SNP', Type = 'SOMATIC', ID = cluster, Sample = cluster) %>%
        dplyr::rename(chrom = chr, pos_start = from) %>%
        rowwise() %>%
        dplyr::mutate(pos_end = pos_start + abs(str_count(ref) - str_count(alt))) %>%
        dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref,alt, Type) %>%
        filter(ref != alt) %>% 
        distinct()
      
      table_int <- int %>% 
        left_join(mut_ids) %>% 
        select(chr,  from, ref, alt, cluster_id_tool_interpreted) %>% 
        dplyr::rename(cluster = cluster_id_tool_interpreted) %>% 
        dplyr::mutate(Project = spn, Genome = 'GRCh38', mut_type = 'SNP', Type = 'SOMATIC', ID = cluster, Sample = cluster) %>%
        dplyr::rename(chrom = chr, pos_start = from) %>%
        rowwise() %>%
        dplyr::mutate(pos_end = pos_start + abs(str_count(ref) - str_count(alt))) %>%
        dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref,alt, Type) %>%
        filter(ref != alt) %>% 
        distinct()
      
      if (!dir.exists(paste0(out, 'output'))){
        write.table(table_raw, file = paste0(out_data_raw, spn,'.txt'), quote = F, sep = '\t', row.names = F)
        write.table(table_int, file = paste0(out_data_int, spn,'.txt'), quote = F, sep = '\t', row.names = F)
        
        
        SigProfilerMatrixGeneratorR(project = spn, genome = "GRCh38", matrix_path = out_data_raw, plot=T)
        SigProfilerMatrixGeneratorR(project = spn, genome = "GRCh38", matrix_path = out_data_int, plot=T)
        
      }
    }
    
    for (context in c('SBS96', 'ID83')){
      print(context)
      out_raw = paste0(out_data_raw, '/', context)
      out_int = paste0(out_data_int, '/', context)
      dir.create(out_raw, showWarnings = F, recursive = T)
      dir.create(out_int, showWarnings = F, recursive = T)
      
      if (context == 'SBS96'){
        data_id_raw = paste0(out_data_raw, '/output/SBS/', spn, '.SBS96.all')
        data_id_int = paste0(out_data_int, '/output/SBS/', spn, '.SBS96.all')
        c_type = "96"
      }  else if (context == 'ID83'){
        data_id_raw = paste0(out_data_raw, '/output/ID/', spn, '.ID83.all')
        data_id_int = paste0(out_data_int, '/output/ID/', spn, '.ID83.all')
        c_type = "83"
      }
      
      if (signature_tool == 'SigProfiler'){
        signature_tool_te = 'sigprofiler'
      } else {
        signature_tool_te = signature_tool
      }
      
      data_signature <- get_tumourevo_signatures(spn = spn, 
                                                 coverage = cov, 
                                                 purity = pur, 
                                                 tool = signature_tool_te,
                                                 vcf_caller = mut_caller,
                                                 cna_caller = cna_caller,
                                                 context = context)
      
      if (signature_tool == 'SigProfiler'){
        fit <- read.table(data_signature$COSMIC_signatures, header = T)
        signature <- colnames(fit)[2:ncol(fit)]  
      } else if (signature_tool == 'BASCULE'){
        fit <- readRDS(data_signature$refined_fit[[1]])
        if (context == 'SBS96'){
          c = 'SBS'
        } else{
          c = 'ID'
        }
        signature <- unique(fit$nmf[[c]]$exposure$sigs)
      }
      
      if (signature_tool == 'SigProfiler'){
        if (context == 'SBS96'){
          cosmic <- read.table('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/GRCh38/SBS_signatures.txt', header = T)
        } else if (context == 'ID83'){
          cosmic <- read.table('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/GRCh38/indel_signatures.txt', header = T)
        }
        subset_cosmic <- cosmic[c('Type', signature)]
        write.table(subset_cosmic, file = paste0(out_data_int, '/subset_cosmic.txt'), quote = F, sep = '\t', row.names = F)
        write.table(subset_cosmic, file = paste0(out_data_raw, '/subset_cosmic.txt'), quote = F, sep = '\t', row.names = F)
        
        # get fit
        cosmic_fit(samples = data_id_raw, 
                   output = out_raw, 
                   input_type='matrix', 
                   context_type=c_type,
                   collapse_to_SBS96=F, 
                   cosmic_version=3.3, 
                   genome_build="GRCh38", 
                   signature_database=paste0(out_data_raw, '/subset_cosmic.txt'),
                   export_probabilities=TRUE,
                   make_plots=TRUE)
        
        cosmic_fit(samples = data_id_int, 
                   output = out_int, 
                   input_type='matrix', 
                   context_type=c_type,
                   collapse_to_SBS96=F, 
                   cosmic_version=3.3, 
                   genome_build="GRCh38", 
                   signature_database=paste0(out_data_int, '/subset_cosmic.txt'),
                   export_probabilities=TRUE,
                   make_plots=TRUE)
        
      } else if (signature_tool == 'BASCULE'){
        out_data_raw_new = paste0(out, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', 'SigProfiler', '_raw/')
        out_data_int_new = paste0(out, spn, '/', cov, 'x_', pur, 'p_', mut_caller, '_', cna_caller, '/', tool, '_', 'SigProfiler', '_int/')
        
        if (context == 'SBS96'){
          data_id_raw = paste0(out_data_raw_new, '/output/SBS/', spn, '.SBS96.all')
          data_id_int = paste0(out_data_int_new, '/output/SBS/', spn, '.SBS96.all')
          
          c_type = "96"
        }  else if (context == 'ID83'){
          data_id_raw = paste0(out_data_raw_new, '/output/ID/', spn, '.ID83.all')
          data_id_int = paste0(out_data_int_new, '/output/ID/', spn, '.ID83.all')
          
          c_type = "83"
        }
        
        if (file.exists(data_id_raw) & file.exists(data_id_int)){
          
          for (type in c('raw', 'int')){
            if ( type == 'raw'){
              data_id = data_id_raw
              out_new = out_raw
            } else {
              data_id = data_id_int
              out_new = out_int
            }
          
            sig_matrix <- read.table(data_id, header = T)
            rownames(sig_matrix) <- sig_matrix$MutationType
            sig_matrix <- sig_matrix %>% select(!MutationType)
            sig_counts <- t(sig_matrix) %>% as.data.frame()
            rownames(sig_counts) <- sub("^[^_]+_(.*)", "\\1", rownames(sig_counts))
            
            if (context == 'SBS96'){
              cat = list(bascule::COSMIC_sbs[signature,] %>% as.data.frame())
              cat[[1]] = cat[[1]][rownames(cat[[1]]) != 'NA',]
              names(cat) = c("SBS")
              sig_counts = sig_counts[rowSums(sig_counts) > 0,]
              input = list("SBS"=sig_counts)
            } else {
              cat = list(bascule::COSMIC_indels[signature,] %>% as.data.frame())
              cat[[1]] = cat[[1]][rownames(cat[[1]]) != 'NA',]
              names(cat) = c("ID")
              sig_counts = sig_counts[rowSums(sig_counts) > 0,]
              input = list("ID"=sig_counts)
            }
            
            
            x = bascule::fit(
              counts = input,
              k_list = 0,
              reference_cat = cat,
              keep_sigs = c("SBS1","SBS5"),
              hyperparameters = NULL,
              lr = 0.005,
              optim_gamma = 0.1,
              n_steps = 3000,
              py = NULL,
              enumer = "parallel",
              nonparametric = TRUE,
              autoguide = FALSE,
              filter_dn = FALSE,
              min_exposure = 0.1,
              CUDA = FALSE,
              compile = FALSE,
              store_parameters = FALSE,
              store_fits = TRUE,
              seed_list = 10
            )
            saveRDS(object = x, file = paste0(out_new,"/bascule_fit_", type, ".rds"))
  
            exp = plot_exposures(x=x, sample_name = T)
            ggplot2::ggsave(filename = paste0(out_new,"/bascule_exposure_", type, ".pdf"),plot=exp, width = 8, height = 8, units = 'in', dpi = 300)
          }
        }
      }
    }
}

