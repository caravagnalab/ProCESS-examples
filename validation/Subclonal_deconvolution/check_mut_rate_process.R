library(ggplot2)
library(tidyverse)
library(ProCESS)
library(CNAqc)
library(mobster)

spn = 'SPN03'
coverage = 100
purity = 0.9
vcf_caller = "mutect2"
cna_caller = "ascat"

tool = "mobster"

source('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/figures/figure3/utils_plot.R')

# SPN_colors = c("01"='steelblue', "02"='seagreen', "03"='goldenrod', 
               # "04"='coral', "05"="magenta4","06"='palevioletred', "07"='indianred3')


github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
setwd(main_path)

source(file.path(github_path, "getters/process_getters.R"))
source(file.path(github_path, "getters/tumourevo_getters.R"))
source(file.path(github_path, "validation/Subclonal_deconvolution/utils_tables.R"))


coverage_list = c(50, 100, 150)
purity_list = c(0.3, 0.6, 0.9)
vcf_caller_list = c("mutect2")
cna_caller_list = c("ascat")
spn_list = paste("SPN", 1:7, sep="0")
combs = expand.grid(coverage=coverage_list,
                    purity=purity_list,
                    vcf_caller=vcf_caller_list,
                    cna_caller=cna_caller_list,
                    spn=spn_list)
# 
# expected_M_list = list()
# M_list = list()

df = data.frame()

for(i in 1:nrow(combs)){
  print(i)
  coverage = combs[i, "coverage"]
  purity = combs[i, "purity"]
  vcf_caller = combs[i, "vcf_caller"]
  cna_caller = combs[i, "cna_caller"]
  spn = combs[i, "spn"]
  
  simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
  
  path_m = file.path(spn, "tumourevo", simulation_id, "subclonal_deconvolution", tool, "SCOUT", spn)
  samples = list.dirs(path_m, recursive=F, full.names=F)  # get_sample_names(spn, base_path=main_path)
  
  print(samples)
  
  for(sample_name in samples) {
    # sample_name = samples[[1]]
    print(file.path(path_m,
                    paste0(sample_name),
                    paste0("SCOUT_", spn, "_", sample_name, "_mobsterh_st_best_fit.rds")))
    
    obj_all = readRDS(file.path(path_m,
                                paste0(sample_name),
                                paste0("SCOUT_", spn, "_", sample_name, "_mobsterh_st_best_fit.rds")))
    obj = obj_all
    
    
    # Compute mutation rate #####
    mutation_rate = NA
    if (obj$fit.tail == TRUE) {
      lq = 0.05
      uq = 0.95
      ploidy = 2
      ncells = 2
      VAFvec <- dplyr::filter(obj$data, cluster == 'Tail', type == 'SNV') %>%
        dplyr::filter(VAF < quantile(VAF, uq) &
                        VAF > quantile(VAF, lq)) %>%
        dplyr::pull(VAF)
      min_vaf = min(VAFvec)
      max_vaf = max(VAFvec)
      M = length(VAFvec)
    
    # append(M_list,M)
    
    segments = obj$data$segment_id %>% unique()
    total_length = sum(sapply(strsplit(segments, ":"), function(x) {
      as.numeric(x[3]) - as.numeric(x[2])
    }))
    
    mut_rate_table = readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/muts_rate.rds')
    
    mu_process = mut_rate_table %>% 
      filter(sample_id == sample_name, coverage == coverage, purity == purity) %>% 
      pull(SNV_rate_process) %>% 
      unique()
    
    expected_M = mu_process * total_length * ((1/min_vaf) - (1/max_vaf))
    # append(expected_M_list, expected_M)
    df = bind_rows(df,
      data.frame(expected = expected_M,
      true = M,
      spn = spn,
      coverage = coverage,
      purity = purity,
      sample = sample_name
    ))
    }
  }
}

df

saveRDS(df, '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/metrics_tables/mu_df.rds')
df = readRDS('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/metrics_tables/mu_df.rds')
ggplot(data = df, aes(x = expected, y = true))+
  geom_point(aes(color=spn))+
  xlab('Number of expected tail mutations with process mu')+
  ylab('Number of mutations in the tail')+
  scale_color_manual(values=SPN_colors)+
  my_ggplot_theme()#+
  # geom_smooth(method = "lm")
 
