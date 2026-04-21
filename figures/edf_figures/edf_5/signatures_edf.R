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




abs_counts_exp_plot <- df_simultated_signatures %>% 
  group_by(sample,signature) %>%
  mutate(median_counts=median(log_counts_sign)) %>%
  select(spn,sample,signature,median_counts,context) %>% 
  unique() %>% 
  # mutate(log_counts=log10(median_counts)) %>%
  ggplot(aes(x = sample, y = median_counts, fill = signature)) +
  geom_col()+
  # scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  labs(
    x = "Name",
    y = "Log10(Mutation Counts)",
    fill = "Signature"
  ) +
  my_ggplot_theme()+
  scale_fill_manual(values = c(sbs_colors,id_colors))+
  facet_nested(
    ~ context + spn,
    scales = "free_x",
    space = "free_x",
  ) +
  xlab("")+
  theme(axis.text.x = element_blank(),
        # strip.text.x.top = element_blank(),
        panel.spacing.x = unit(0.2, "lines"),
        legend.box = "horizontal",
        panel.grid.minor = element_blank(),
        # panel.grid = 
        axis.ticks.x = element_blank(),axisticks = element_blank())+
  guides(col  = guide_legend(nrow = 3))

abs_counts_exp_plot



cosine_samples <-df_all_combs_SPN_signatures %>% 
  separate(col = sample,into = c("spn_id","sample_small"),sep = "_",remove = F) %>% 
  ggplot(aes(x = sample_small, y = cosine, color = caller)) +
  geom_boxplot() +
  facet_nested(
    ~ context + spn,
    scales = "free",
    space = "free_x"
  ) +
  scale_color_manual("Tool",
    values = c(
      'SigProfiler' = 'sienna1',
      'BASCULE'     = 'dodgerblue4'
    )
  ) +
  my_ggplot_theme() +
  labs(x = "Sample ID", y = "Cosine similarity") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.x.top = element_blank(),panel.spacing.x = unit(0.2, "lines")
  )+
  guides(
    col  = guide_legend(nrow = 2)
  )

subclonal_arch<-readRDS(file="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure5/signatures_cohort/sample_classification_clonal_heterogeneity.rds")
jaccard_plot <-subclonal_arch %>% 
  select(sample,context,median_jaccard_similary_clusters,sample_class) %>% 
  distinct() %>% 
  ggplot(aes(x=sample,y="",fill=median_jaccard_similary_clusters))+
  geom_tile()+
  scale_fill_continuous("Jaccard similarity",palette = rev(c("#FEE0D2", "#FC9272", "#DE2D26")))+
  my_ggplot_theme()+
  facet_wrap(~context)+
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.text.x.top = element_blank(),panel.spacing.x = unit(0.2, "lines")
  )

cosine_plot <-subclonal_arch %>% 
  select(sample,context,mean_cosine_similary_clusters,sample_class) %>% 
  distinct() %>% 
  ggplot(aes(x=sample,y="",fill=mean_cosine_similary_clusters))+
  geom_tile()+
  scale_fill_continuous("Cosine similarity",palette = c("#E5F5E0", "#A1D99B", "#31A354"))+
  my_ggplot_theme()+
  facet_wrap(~context)+
  theme(
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_blank(),
    strip.text.x.top = element_blank(),panel.spacing.x = unit(0.2, "lines")
  )
clonal_subcl_exp <- subclonal_arch %>% 
  filter(coverage==150, purity==0.9) %>% 
  ggplot(aes(
    x = interaction(cluster_id_process,sample),  # ← key change
    y = exposure, 
    fill = causes
  )) +
  geom_col() +
  labs(
    x = "Name",
    y = "Exposure",
    fill = "Signature"
  ) +
  my_ggplot_theme() +
  scale_fill_manual(values = c(sbs_colors, id_colors)) +
  facet_nested(
    ~ context + spn ,
    scales = "free_x",
    space = "free_x"
  ) +
  xlab("") +
  theme(
    axis.text.x = element_blank(),
    panel.spacing.x = unit(0.2, "lines"),
    legend.box = "horizontal",
    panel.grid.minor = element_blank(),
    axis.ticks.x = element_blank(),
    axisticks = element_blank()
  ) +
  guides(col = guide_legend(nrow = 3))

cosine_samples <- cosine_samples + theme(plot.margin = margin(0, 0, 0, 0),panel.spacing = unit(0.05, "lines"))
abs_counts_exp_plot <- abs_counts_exp_plot + theme(plot.margin = margin(0, 0, 0, 0),panel.spacing = unit(0.05, "lines"))
plt_signatures_exp <- wrap_plots(list(clonal_subcl_exp,cosine_samples,jaccard_plot,cosine_plot),
           design = "AAAA\nCCCC\nDDDD\nBBBB\nBBBB\nBBBB",guides="collect")+
  plot_layout(heights = c(1, 0.6, 0.6, 1)) &
  theme(legend.position = "bottom")
ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/edf_figures/edf_5/edf5.pdf",plot = plt_signatures_exp,
       height = 5,width = 10)
