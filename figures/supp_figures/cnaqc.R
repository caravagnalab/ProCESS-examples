rm(list = ls())
options(bitmapType='cairo')
library(optparse)
library(tidyverse)
library(gt)

# source('/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/utils.R')

source('/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R')
source('/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples//getters/tumourevo_getters.R')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
df_all_combs_SPN_cna<- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/cna/cna_all.rds")
df_all_combs_SPN_cna <- df_all_combs_SPN_cna %>% 
  dplyr::rename(cna_caller=tool) %>% 
  mutate(cna_caller=case_when(cna_caller=="Sequenza"~"sequenza",
                              cna_caller=="ASCAT" ~ "ascat",
                              TRUE~cna_caller))

fga_df <-readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/fga_df.rds") %>% 
  dplyr::rename(sample=spn) %>% 
  mutate(fga_class=case_when(sample%in%c("SPN01_1.1","SPN01_1.3","SPN06_3.2")~"WGD",
                             TRUE ~ fga_class))
### select columns of interest

df_all_combs_SPN_cna <- df_all_combs_SPN_cna %>% 
  select(sample,spn,coverage,true_purity,class,cna_caller) %>% 
  mutate(true_purity=as.numeric(true_purity),
         coverage=as.numeric(coverage))
SPN <- paste0('SPN0', c(1,3,4,5,6,7))
coverages <- c(50, 100, 150)
purities <- c(0.3, 0.6, 0.9)
cna_callers <- c("sequenza","ascat")
variant_callers <- c("mutect2","strelka")
params_grid = expand.grid(coverages, purities,cna_callers,variant_callers)
colnames(params_grid) = c("coverage", "purity","cna_caller","snv_caller")
validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION/"
validate_spns <- list()
for (spn in SPN){
  print(spn)
  validate_tables <- list()
  for (i in 1:nrow(params_grid)){
    print(i)
    cov <- params_grid$coverage[i]
    pur <- params_grid$purity[i]
    cna_caller <- params_grid$cna_caller[i]
    snv_caller <- params_grid$snv_caller[i]
    
    cnaqc_file <- file.path(validation_dir,spn,"cnaqc",paste0(cov, 'x_', pur),paste0(snv_caller, "_",cna_caller),"cnaqc_validate.rds")
    if (file.exists(cnaqc_file)){
      validate_tables[[i]] <- readRDS(cnaqc_file)
    } else{
      print(paste0('Missing file for: ', spn, '-', cov, 'x-', pur, 'p', snv_caller, "-",cna_caller))
    }
    
  }
  validate_spns[[spn]] <- do.call("bind_rows",validate_tables)
}
cnaqc_df <-do.call("bind_rows",validate_spns)
cnaqc_df <- cnaqc_df %>% 
  left_join(y = df_all_combs_SPN_cna)
cnaqc_df_filt <- cnaqc_df %>%
  mutate(comb_te=paste0(vcf_caller,"_",cna_caller)) %>% 
  mutate(CNAqc = ifelse(is.na(CNAqc), 'NA', CNAqc)) %>% 
  group_by(sample, coverage, purity, true_karyo,comb_te) %>%
  filter(!(all(CNAqc == "Subclonal", na.rm = T))) %>%
  ungroup() %>% 
  left_join(fga_df,by = "sample") %>%
  mutate(comb_sarek=paste(coverage,true_purity,sep="_")) %>%
  mutate(
    CNAqc_superclass = case_when(
      str_detect(CNAqc, "CNAqc FAIL") ~ "CNAqc_FAIL",
      TRUE ~ "CNAqc_OK"
    )) 

cnaqc_df_filt_tot_counts <- cnaqc_df_filt %>%
  select(true_purity, class,cna_caller,sample) %>%
  unique() %>%
  group_by(true_purity, class,cna_caller) %>%
  # group_by(true_purity, class,cna_caller) %>%
  summarize(n = n(), .groups = "drop")



cnaqc_df_filt %>% 
  group_by(cna_caller,CNAqc_superclass,class,true_purity) %>%
  summarise(n = n())  %>% 
  ggplot() + 
  geom_col(aes(x = class, y=n, fill = CNAqc_superclass), position = "fill") +
  ggplot2::scale_fill_manual('', values = c(
    `CNAqc_OK` = 'seagreen3',
    `CNAqc_FAIL` = 'indianred3',
    `Subclonal` = 'gainsboro',
    `NA` = 'gray60'
  )) +

  geom_text(
    data = cnaqc_df_filt_tot_counts,
    aes(x = class, y = 0.9, label = paste0("n= ",n)),
    vjust = -0.3,
    inherit.aes = FALSE,
    size = 3,col="black"
  ) +
  ylab('% of segments') +
  # ggh4x::facet_nested(comb_te ~ sample+coverage)+
  facet_grid(cna_caller ~ true_purity, scales = "free") +
  my_ggplot_theme()+
  coord_cartesian(clip = "off")+
  theme(legend.position = "bottom")


##########################################
cnaqc_df_filt %>% 
  # filter(spn=="SPN04") %>% 
  group_by(cna_caller,CNAqc_superclass,class,true_purity,sample) %>%
  mutate(n_superclass = n())  %>% 
  group_by(cna_caller,class,true_purity,sample) %>% 
  mutate(n_total = n())  %>% 
  mutate(freq_superclass=n_superclass/n_total) %>% 
  select(c("spn","cna_caller","class","true_purity","sample","freq_superclass","CNAqc_superclass","coverage","fga_class")) %>% 
  unique() %>% 
  pivot_wider(
    names_from = CNAqc_superclass,
    values_from = freq_superclass,
    names_prefix = "freq_"
  ) %>% 
  unique() %>% 
  mutate(sample_cnaqc_class=case_when(freq_CNAqc_OK>=0.8~"CNAqc_sample_OK",
                                      TRUE~"CNAqc_sample_FAIL" ))  %>% 
  filter(spn=="SPN04") %>% 
  ggplot(aes(x = sample, y = cna_caller, fill = sample_cnaqc_class,color=class)) +
  geom_tile(linewidth = 1,alpha=0.4) +
  labs(x = "Sample", y = "CNA caller", fill = "Class",color="CN caller performance",caption = "A sample is CNAqc OK when more than 80% of the segments are CNAqc-OK") +
  scale_fill_manual(values=c("CNAqc_sample_OK"="seagreen3","CNAqc_sample_FAIL"="indianred3")) + 
  scale_color_manual(values=c("correct"="forestgreen","uncorrect purity"="orange","uncorrect ploidy"="purple")) + 
  ggh4x::facet_nested(~ true_purity+coverage)+
  my_ggplot_theme()+
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )







#ggsave(filename = 'validation_CNAqc.png', plot = plt, width = 18, height = 7, units = 'in', dpi = 400)
