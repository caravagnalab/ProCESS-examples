library(dplyr)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(ggh4x)
purities <- c(0.3,0.6,0.9)
coverages <- c(50, 100,150)
contexts <- c("SBS96","ID83")
all_combs <- list()
all_metrics <- list()
all_cosine  <- list()
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN05", "SPN06", 'SPN07')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION/"

cna_caller <- "ascat"
variant_caller <- "mutect2"

mut_counts_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/mutations_counts/"
params_grid = expand.grid(coverages, purities)
colnames(params_grid) = c("coverage", "purity")
df_all_SPN_muts_count <- list()
for (SPN in SPNS){
  spn <- SPN
  df_all_SPN_muts_count[[spn]] = lapply(1:nrow(params_grid), function(i) {
    coverage = params_grid[i,]$coverage
    purity = params_grid[i,]$purity
    comb <- paste0(coverage,"x_",purity,"p")
    mut_counts_dir_spn <- file.path(mut_counts_dir,spn)
    all_metrics_counts <- list()
    df <- readRDS(file.path(mut_counts_dir_spn,paste0("mutations_counts_",comb,".rds")))
  }) %>% bind_rows()
}
df_all_SPN_muts_count <- do.call("rbind",df_all_SPN_muts_count) %>% 
  mutate(context=case_when(type=="SNV"~"SBS96",
                           TRUE ~ "ID83"))


gt_signatures <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/signature_exposures/"


df_simultated_signatures <- lapply(SPNS, function(spn){
  
  exp_all <- tibble()
  
  for (pur in purities){
    for (cov in coverages){
      
      comb <- paste0(cov, "x_", pur, "p")
      
      ## SBS
      exp_sbs <- readRDS(
        file.path(gt_signatures, spn,
                  paste0("signature_exposure_", comb, "_SBS.rds"))
      ) %>%
        as.data.frame()
      
      colnames(exp_sbs) <- get_sample_names(spn)
      
      exp_sbs <- exp_sbs %>%
        rownames_to_column(var = "signature") %>%
        pivot_longer(
          cols = starts_with("SPN"),
          names_to = "sample",
          values_to = "exposure"
        ) %>%
        mutate(
          coverage = cov,
          purity   = pur,
          type     = "SNV",
          spn = spn
        )
      
      ## ID
      exp_id <- readRDS(
        file.path(gt_signatures, spn,
                  paste0("signature_exposure_", comb, "_ID.rds"))
      ) %>%
        as.data.frame()
      
      colnames(exp_id) <- get_sample_names(spn)
      
      exp_id <- exp_id %>%
        rownames_to_column(var = "signature") %>%
        pivot_longer(
          cols = starts_with("SPN"),
          names_to = "sample",
          values_to = "exposure"
        ) %>%
        mutate(
          coverage = cov,
          purity   = pur,
          type     = "INDEL",
          spn=spn
        )
      
      ## combine + accumulate
      exp_all <- bind_rows(exp_all, exp_sbs, exp_id)
    }
  }
  
  return(exp_all)
}) %>%
  bind_rows()



df_all_combs_SPN_signatures <- lapply(SPNS,function(spn){
  validation_dir_spn <- paste0(validation_dir,spn,"/signature/")
  all_combs <- list()
  all_metrics <- list()
  all_cosine  <- list()
  all_metrics_sample <-list()
  for (purity in purities){
    for (coverage in coverages){
      comb <- paste0(coverage,"x_",purity)
      metrics_spn <- list()
      cosine_spn  <- list()
      metrics_sample <- list()
      comb_caller <- paste0(variant_caller,"_",cna_caller)
      for (ctx in contexts){
        rds_metrics <- file.path(validation_dir_spn,comb,comb_caller,paste0("metrics_",ctx,"_spn.rds"))
        rds_metrics_sample <- file.path(validation_dir_spn,comb,comb_caller,paste0("metrics_",ctx,"_sample.rds"))
        rds_cosine  <- file.path(validation_dir_spn,comb,comb_caller,paste0("cosine_mse_",ctx,".rds"))
        
        metrics_spn[[ctx]] <- readRDS(rds_metrics) %>% 
          dplyr::mutate(context=ctx,
                        purity=purity,
                        coverage=coverage)
        metrics_sample[[ctx]] <- readRDS(rds_metrics_sample) %>% 
          dplyr::mutate(context=ctx,
                        purity=purity,
                        coverage=coverage) %>% 
          unique()
        cosine_spn[[ctx]] <- readRDS(rds_cosine) %>% 
          dplyr::mutate(context=ctx,
                        purity=purity,
                        coverage=coverage)
      }
      
      # bind SBS96 + ID83 for this purity/coverage
      metrics_df <- dplyr::bind_rows(metrics_spn)
      cosine_df  <- dplyr::bind_rows(cosine_spn)
      metrics_sample_df <- dplyr::bind_rows(metrics_sample)
      # save into global list
      all_metrics[[comb]] <- metrics_df
      all_cosine[[comb]]  <- cosine_df
      all_metrics_sample[[comb]] <- metrics_sample_df
    }
  }
  cosine_metrics <- all_cosine %>% bind_rows()
  # metrics <- all_metrics %>% bind_rows()
  metrics_ss <- all_metrics_sample %>% bind_rows()
  final_df <- inner_join(cosine_metrics,metrics_ss)
  # all <- list("cosine"=cosine_metrics,"general_metrics"=metrics, "sample_metrics"=metrics_ss)
}) %>% bind_rows()

df_all_combs_SPN_signatures <-df_all_combs_SPN_signatures %>% 
  mutate(spn=case_when(spn=="SPN07_last"~"SPN07",
                       TRUE ~spn)) %>% 
  mutate(sample = str_replace_all(sample, "last_", ""))

df_simultated_signatures <-df_simultated_signatures %>% 
  inner_join(df_all_SPN_muts_count) %>% 
  mutate(log_total_counts=log10(mutation_count)) %>% 
  mutate(log_counts_sign=log_total_counts*exposure) %>% 
  mutate(abs_counts_sign=mutation_count*exposure) %>% 
  group_by(sample,type) %>% 
  mutate(mean_mutation_count=mean(mutation_count)) 

df_simultated_signatures_sample_spec <- df_simultated_signatures %>%
  select(spn,sample,coverage,purity,signature,exposure,mutation_count,context) %>%
  filter(exposure>0) %>%
  group_by(sample,coverage,purity,context) %>%
  mutate(n_signatures=n()) %>%
  mutate(muts_per_signature=mutation_count/n_signatures) %>%
  mutate(log10_muts_per_signature=log10(muts_per_signature)) %>%
  select(spn,sample,coverage,purity,muts_per_signature) %>%
  distinct()

df_all_combs_SPN_signatures <- df_all_combs_SPN_signatures %>% 
  left_join(y = df_simultated_signatures_sample_spec,
            by = c("spn","sample","coverage","purity","context"),
            relationship = "many-to-many") 

plt_signatures <- df_all_combs_SPN_signatures %>% 
  group_by(context) %>% 
  mutate(mean_cosine=mean(cosine)) %>% 
  mutate(purity = as.character(purity)) %>% 
  arrange(desc(muts_per_signature)) %>% 
  mutate(log10_muts_count = log10(muts_per_signature)) %>% 
  ggplot(aes(x = log10_muts_count,
             y = cosine,
             group = interaction(caller, context))) +
  geom_point(alpha=0.4,aes(color=spn)) +
  scale_color_manual('SPN', values = SPN_colors)+
  ggnewscale::new_scale_color() +
  geom_smooth(span = 0.1,se = F,aes(fill=caller,color=caller),alpha = 0.1,method = "loess") +
  scale_color_manual(values = c('SigProfiler'='sienna1',
                                'BASCULE'='dodgerblue4'))+
  geom_hline(aes(yintercept = mean_cosine),linetype = "dashed",color="grey40")+
  facet_wrap(~context,scales = "free",nrow = 2,strip.position = "right")+
  my_ggplot_theme()+
  xlab("Log10(Median Mutation count per signature)")+
  ylab("Cosine Similarity")  +
  guides(col  = guide_legend(alpha = 1)) 
plt_signatures
