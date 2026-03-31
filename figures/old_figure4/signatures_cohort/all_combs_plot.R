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
# plot_all_spns <- lapply(spns,function(spn){
#   validation_dir_spn <- paste0(validation_dir,spn,"/signature/")
#   all_combs <- list()
#   all_metrics <- list()
#   all_cosine  <- list()
#   for (purity in purities){
#     for (coverage in coverages){
#       comb <- paste0(coverage,"x_",purity)
#       metrics_spn <- list()
#       cosine_spn  <- list()
#       comb_caller <- paste0(variant_caller,"_",cna_caller)
#       for (ctx in contexts){
#         rds_metrics <- file.path(validation_dir_spn,comb,comb_caller,paste0("metrics_",ctx,"_spn.rds"))
#         rds_cosine  <- file.path(validation_dir_spn,comb,comb_caller,paste0("cosine_mse_",ctx,".rds"))
#         
#         metrics_spn[[ctx]] <- readRDS(rds_metrics) %>% 
#           dplyr::mutate(context=ctx,
#                         purity=purity,
#                         coverage=coverage)
#         
#         cosine_spn[[ctx]] <- readRDS(rds_cosine) %>% 
#           dplyr::mutate(context=ctx,
#                         purity=purity,
#                         coverage=coverage)
#       }
#       
#       # bind SBS96 + ID83 for this purity/coverage
#       metrics_df <- dplyr::bind_rows(metrics_spn)
#       cosine_df  <- dplyr::bind_rows(cosine_spn)
#       
#       # save into global list
#       all_metrics[[comb]] <- metrics_df
#       all_cosine[[comb]]  <- cosine_df
#     }
#   }
#   
#   # final dataframes
#   final_metrics <- dplyr::bind_rows(all_metrics)
#   final_cosine  <- dplyr::bind_rows(all_cosine)
  color_caller <- list("ID83"=c('BASCULE' = '#E1BFF8','SigProfiler' = '#F6AF92'),
                       "SBS96"=c('BASCULE' = 'purple1','SigProfiler' = 'orange2'))
#   
#   plots <- list()
#   for (ctx in contexts){
#     
#     cosine_plot <-final_cosine %>%
#       filter(context==ctx) %>% 
#       ggplot(aes(x = sample, y = cosine, fill = caller)) +
#       geom_col(position = position_dodge()) +
#       geom_errorbar(
#         aes(ymin = cosine - mse, ymax = cosine + mse),
#         position = position_dodge(width = 0.9),
#         width = 0.3
#       ) +
#       scale_fill_manual(values = color_caller[[ctx]]) +
#       facet_grid(as.numeric(coverage) ~ as.numeric(purity)) +
#       theme_bw() +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#       ylab(label = paste0("cosine ",ctx))
#     
#     
#     metrics_plot <- final_metrics %>%
#       dplyr::filter(context==ctx) %>% 
#       ggplot(aes(x = metric, y = mean, fill = caller)) +
#       geom_col(position = position_dodge()) +
#       geom_errorbar(
#         aes(x=metric,ymin = mean - sd, ymax = mean + sd),
#         position = position_dodge(width = 0.9),
#         width = 0.3
#       ) +
#       scale_fill_manual(values = color_caller[[ctx]]) +
#       facet_grid(as.numeric(coverage) ~ as.numeric(purity)) +
#       theme_bw() +
#       theme(axis.text.x = element_text(angle = 45, hjust = 1))+
#       ylab(label = paste0("mean ",ctx))
#     
#     plots[[ctx]] <- wrap_plots(list(cosine_plot,metrics_plot))+plot_layout(guides = "collect")
#   }
# 
#    wrap_plots(plots,nrow = 2)+plot_annotation(title = spn)
# 
# })
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

#df_all_combs_SPN_signatures_cosine <- all$cosine %>% bind_rows()
#df_all_combs_SPN_signatures_metrics_ss <- all$sample_metrics %>% bind_rows()

# df_all_combs_SPN_signature_exposures<- lapply(SPNS,function(spn){
#   validation_dir_spn <- paste0(validation_dir,spn,"/signature/")
#   all_combs <- list()
#   for (purity in purities){
#     for (coverage in coverages){
#       comb <- paste0(coverage,"x_",purity,"p")
#       metrics_spn <- list()
#       
#       comb_caller <- paste0(variant_caller,"_",cna_caller)
#       contexts_all <- c("SBS96","ID83")
#       context_classes <-  gsub('[[:digit:]]+', '', contexts_all) 
#       gt_exposure <- lapply(context_classes, function(c){
#         d=readRDS(file.path("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/signature_exposures/",spn,
#                             paste0("signature_exposure_",coverage,"x_",purity,"p_",c,".rds")))
#         colnames(d) <- get_sample_names(spn)
#         td=t(d)
#       })
#       names(gt_exposure)<-context_classes
#       for (ctx in contexts_all){
#         context_classes <- gsub('[[:digit:]]+', '', ctx)
#         # ground_truth_nested <- gt_exposure[[context_classes]]  %>%
#         #   tibble::column_to_rownames("Sample_ID") %>% 
#         #   as.matrix()
#         # ground_truth_nested[is.na(ground_truth_nested)] <- 0
#         ground_truth_nested <- gt_exposure[[context_class]]
#         result <- tryCatch({
#           
#           # Get SigProfiler paths
#           sigprofiler <- get_tumourevo_signatures(
#             spn = spn,
#             coverage = coverage,
#             purity = purity,
#             vcf_caller = vcf_caller,
#             cna_caller = cna_caller,
#             tool = "SigProfiler",
#             context = context
#           )
#           
#           bascule <- get_tumourevo_signatures(
#             spn = spn,
#             coverage = coverage,
#             purity = purity,
#             vcf_caller = vcf_caller,
#             cna_caller = cna_caller,
#             tool = "BASCULE",
#             context = context_classes
#           )
#           # Combine and load paths
#           paths <- c(
#             
#             sigprofiler$COSMIC_exposure,
#             sigprofiler$COSMIC_signatures,
#             sigprofiler$denovo_exposure,
#             sigprofiler$denovo_signatures,
#             bascule$refined_fit
#           )
#           
#           data <- load_signature_data(paths)
#           
#           names(data) <- c(
#             "SigProfiler_COSMIC_exposure",
#             "SigProfiler_COSMIC_signatures",
#             "SigProfiler_denovo_exposure",
#             "SigProfiler_denovo_signatures",
#             "BASCULE_refined_fit"
#           )
#           data
#           
#         })
#         
#         if (!is.null(result)) {
#           tumourevo_signature_res <- result
#         }
#         ground_truth_nested <- as.data.frame(ground_truth_nested)
#         ground_truth_nested %>% pivot_longer(cols = starts_with(context_classes),
#                                              names_to = "Signature",
#                                              values_to = "Exposure")
#         sigprof_aligned <- align_callers(tumourevo_signature_res = tumourevo_signature_res,tool = "SigProfiler",spn = spn) %>% as.data.frame()
#         bascule_aligned <- align_callers(tumourevo_signature_res = tumourevo_signature_res,tool = "BASCULE",spn = spn) %>% as.data.frame()
#         
#         df_long <- bind_rows(
#           ground_truth_nested %>% mutate(sample=row.names(.), Method="ProCESS"),
#           bascule_aligned %>% mutate(sample=row.names(.), Method="BASCULE"),
#           sigprof_aligned %>% mutate(sample=row.names(.), Method="SigProfiler")
#         ) %>%
#           pivot_longer(cols = starts_with(context_classes),
#                        names_to = "Signature",
#                        values_to = "Exposure")
#         metrics_spn[[ctx]] <- df_long
#       }
#       
#       # bind SBS96 + ID83 for this purity/coverage
#       metrics_df <- dplyr::bind_rows(metrics_spn)
#       
#       # save into global list
#       all_metrics[[comb]] <- metrics_df
#     }
#   }
#   
#   metrics <- all_metrics %>% bind_rows()
# }) %>% bind_rows()


df_all_combs_SPN_signatures <-df_all_combs_SPN_signatures %>% 
  mutate(spn=case_when(spn=="SPN07_last"~"SPN07",
                       TRUE ~spn)) %>% 
  mutate(sample = str_replace_all(sample, "last_", "")) 
  # left_join(df_all_SPN_muts_count)
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
  select(spn,sample,coverage,purity,muts_per_signature) %>%
  distinct()

# df_simultated_signatures_sample_spec <- df_simultated_signatures %>% 
#   # select(spn,sample,coverage,purity,signature,exposure,mutation_count,context) %>% 
#   filter(exposure>0) %>% 
#   group_by(sample,coverage,purity,context) %>% 
#   mutate(min_mut_signature=min(exposure))


##### APPLY CLASSIFICATION ###
df_simultated_signatures <-df_simultated_signatures %>% 
  filter(exposure>0) %>% 
  mutate(exposure_rounded=round(exposure,1)) %>% 
  mutate(sample_class=case_when(context=="SBS96" & mean_mutation_count<4000 ~"low mutations",
                                 context=="ID83" & mean_mutation_count<1000 ~"low mutations",
                                 TRUE~"high mutations"))
  # mutate(sample_class2=case_when(sample_class1=="high mutations" & context=="SBS96" & abs_counts_sign<1000 ~ "low exposure",
  #                                sample_class1=="high mutations" & context=="ID83" & abs_counts_sign<500 ~ "low exposure",
  #                                TRUE ~ "high exposure")) 
  # mutate(exposure_class=case_when(context=="SBS96" & abs_counts_sign<400 ~"low counts",
  #                                 context=="ID83" & abs_counts_sign<100 ~"low counts",
  #                                 TRUE~"high counts")) %>%
  # group_by(sample, context) %>%
  # mutate(
  #   n_low_counts = sum(exposure_class == "low counts")
  # ) %>%
  # ungroup() %>% 
  # mutate(sample_exposure_class=case_when(n_low_counts>2~"low exposure",
  #                                        TRUE~"normal exposure"))
  
  # mutate(sample_exposure_class_1=case_when(sample%in%c("SPN03_1.1","SPN03_2.1","SPN03_3.1",
  #                                                      "SPN03_4.1","SPN04_2.1","SPN06_3.1",
  #                                                      "SPN06_3.2","SPN07_3.1","SPN07_3.2","SPN07_3.3")~"low exposure",
  #                                        TRUE~"normal exposure"))




df_all_combs_SPN_signatures <- df_all_combs_SPN_signatures %>% 
  left_join(y = df_simultated_signatures_sample_spec,
            by = c("spn","sample","coverage","purity","context"),
            relationship = "many-to-many") 
# df_all_combs_SPN_signatures %>% 
#   ggplot(aes(x=log10(mutation_count),y=cosine,fill=caller,color=caller,shape=context)) +
#   geom_point(size=2) +
#   scale_color_manual(values = c('BASCULE' = 'purple1','SigProfiler' = 'orange2'))+
#   my_ggplot_theme()
# 
# 
# df_all_combs_SPN_signatures %>% 
#   mutate(number_signatures=as.factor(TP+FN)) %>% 
#   # ggplot(aes(x=number_signatures,y=cosine,fill=caller)) +
#   ggplot(aes(x=cosine,fill=caller,color=caller)) +
#   # geom_boxplot()+
#   geom_histogram(alpha=0.5,binwidth = 0.01,position = 'identity')+
#   # geom_density(alpha=0.1)+
#   scale_fill_manual(values = c('BASCULE' = 'purple1','SigProfiler' = 'orange2'))+
#   scale_color_manual(values = c('BASCULE' = 'purple1','SigProfiler' = 'orange2'))+
#   my_ggplot_theme()+
#   facet_wrap(~context)
# 
# df_all_combs_SPN_signatures %>% 
#   mutate(number_signatures=as.factor(TP+FN)) %>% 
#   ggplot(aes(x=number_signatures,y=mse,fill=caller)) +
#   # ggplot(aes(x=spn,y=cosine,fill=caller)) +
#   geom_boxplot()+
#   scale_fill_manual(values = c('BASCULE' = 'purple1','SigProfiler' = 'orange2'))+
#   my_ggplot_theme()+
#   facet_wrap(~context)
# 
# df_all_combs_SPN_signatures %>% 
#   filter(context=="SBS96") %>% 
#   mutate(number_signatures=(TP+FN)) %>%
#   ggplot(aes(y=cosine,x=spn,fill=spn)) +
#   geom_col()+
#   scale_fill_manual(values =SPN_colors)+
#   my_ggplot_theme()




df_all_combs_SPN_signatures %>% 
  group_by(context) %>% 
  mutate(mean_cosine=mean(cosine)) %>% 
  mutate(purity = as.character(purity)) %>% 
  arrange(desc(muts_per_signature)) %>% 
  mutate(log10_muts_count = log10(muts_per_signature)) %>% 
  ggplot(aes(x = log10_muts_count,
             y = cosine,
            color = caller,
             # linetype = context,
             group = interaction(caller, context))) +
  # geom_line() 
  scale_color_manual(values = c('SigProfiler'='sienna1',
  'BASCULE'='dodgerblue4'))+
  scale_fill_manual(values = c('SigProfiler'='sienna1',
                                'BASCULE'='dodgerblue4'))+
  geom_smooth(span = 0.1,se = F,aes(fill=caller),alpha = 0.1,method = "loess") +
  ggnewscale::new_scale_color() +
  geom_point(alpha=0.4,aes(color=spn)) +
  scale_color_manual(values = SPN_colors)+
  geom_hline(aes(yintercept = mean_cosine),linetype = "dashed",color="grey40")+
  facet_wrap(~context,scales = "free",nrow = 2)+
  my_ggplot_theme()

df_all_combs_SPN_signatures %>% 
  group_by(context) %>% 
  mutate(mean_cosine=mean(cosine)) %>% 
  mutate(purity = as.character(purity)) %>% 
  arrange(desc(muts_per_signature)) %>% 
  mutate(log10_muts_count = log10(muts_per_signature)) %>% 
  ggplot(aes(x = spn,
             y = cosine,
             color = caller))+
             # linetype = context,
             # group = interaction(caller, context))) +
  geom_boxplot()+
  scale_color_manual(values = c('SigProfiler'='sienna1',
                                'BASCULE'='dodgerblue4'))+
  facet_wrap(~context,scales = "free")+
  my_ggplot_theme()

df_all_combs_SPN_signatures %>% 
  group_by(context) %>% 
  mutate(mean_cosine=mean(cosine)) %>% 
  mutate(purity = as.character(purity)) %>% 
  arrange(desc(muts_per_signature)) %>% 
  mutate(log10_muts_count = log10(muts_per_signature)) %>% 
  ggplot(aes(x = context,
             y = cosine))+
  geom_boxplot(outliers = F)+
  my_ggplot_theme()

df_all_combs_SPN_signatures %>% 
  mutate(purity = as.character(purity)) %>% 
  # filter(spn!="SPN07") %>%
  arrange(desc(min_mut_signature)) %>%
  # mutate(log10_muts_count = log10(muts_per_signature)) %>% 
  ggplot(aes(x = min_mut_signature,
             y = cosine,
             color = caller,
             group = caller)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c('SigProfiler'='sienna1',
                                'BASCULE'='dodgerblue4'))+
  # geom_smooth(span = 0.1,se = F,aes(fill=caller),alpha = 0.1) +
  scale_fill_manual(values = c('SigProfiler'='sienna1',
                               'BASCULE'='dodgerblue4'))+
  facet_wrap(~context,scales = "free_x",nrow = 2)+
  
  my_ggplot_theme()

df_all_combs_SPN_signatures %>% 
  mutate(purity = as.character(purity)) %>% 
  filter(spn!="SPN07") %>%
  arrange(desc(muts_per_signature)) %>% 
  mutate(log10_muts_count = log10(muts_per_signature)) %>% 
  ggplot(aes(x = log10_muts_count,
             y = cosine,
             color = caller,
             group = caller)) +
  # geom_line() +
  # geom_point() +
  scale_color_manual(values = c('SigProfiler'='sienna1',
                                'BASCULE'='dodgerblue4'))+
  geom_smooth(span = 0.1,se = F,aes(fill=caller),alpha = 0.1) +
  scale_fill_manual(values = c('SigProfiler'='sienna1',
                               'BASCULE'='dodgerblue4'))+
  facet_wrap(~context,scales = "free_x",nrow = 2)+
  
  my_ggplot_theme()
  
library(ggpubr)

df_all_combs_SPN_signatures %>% 
  mutate(purity = as.character(purity)) %>% 
  ggplot(aes(x = sample_class, y = cosine, color = caller)) +
  geom_boxplot(outliers = FALSE) +
  facet_grid(.~context, scales = "free_x") +
  stat_compare_means(
    comparisons = list(c("low mutations", "high mutations")), 
    method = "wilcox.test",                # or "t.test" if appropriate
    label = "p.format"
  ) +
  scale_color_manual(values = c('SigProfiler'='sienna1',
                                'BASCULE'='dodgerblue4')) +
  my_ggplot_theme() +
  xlab("") 
  # theme(axis.text.x = element_text(angle = 45, hjust = 1))







df_all_combs_SPN_signatures %>% 
  select(spn,sample,sample_exposure_class) %>% View()


df_all_combs_SPN_signatures %>% 
  ggplot(aes(x=sample,y=cosine,color=sample_exposure_class))+
  geom_point()+
  facet_grid(context~caller)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1))





  
exp_plot <- df_simultated_signatures %>% 
  group_by(sample,signature) %>%
  mutate(median_exposure=median(exposure)) %>%
  mutate(context = sub("[0-9].*$", "", signature)) %>% 
  ggplot(aes(x = sample, y = median_exposure, fill = signature)) +
  geom_bar(stat = "identity", position = "fill") +
  # scale_y_continuous(labels = scales::percent_format(accuracy = 0.01)) +
  labs(
    x = "Name",
    y = "Signature contribution (%)",
    fill = "Signature"
  ) +
  my_ggplot_theme()+
  scale_fill_manual(values = c(sbs_colors,id_colors))+
  facet_wrap(~context)+
  # facet_grid(context~spn,scales = "free")+
  xlab("")+
  theme(axis.text.x = element_blank(),
        strip.text.x.top = element_blank(),
        axis.ticks.x = element_blank(),axis.ticks = element_blank())

abs_counts_exp_plot <- df_simultated_signatures %>% 
  group_by(sample,signature) %>%
  mutate(median_counts=median(log_counts_sign)) %>%
  mutate(context = sub("[0-9].*$", "", signature)) %>% 
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
    space = "free_x"
  ) +
  # facet_wrap(~context,scales = "free_y")+
  # facet_grid(context~spn,scales = "free")+
  xlab("")+
  theme(axis.text.x = element_blank(),
        # strip.text.x.top = element_blank(),
        panel.spacing.x = unit(0.2, "lines"),
        legend.box = "horizontal",
        axis.ticks.x = element_blank(),axisticks = element_blank())+
  guides(col  = guide_legend(ncol = 6))

abs_counts_exp_plot
cosine_samples <- df_all_combs_SPN_signatures %>% 
  ggplot(aes(x=sample,y=cosine,color=caller))+
  geom_boxplot()+
  facet_wrap(~context,scales = "free_x")+
  #ggh4x::facet_nested(~spn+context,scales = "free")+
  scale_color_manual(values=c('SigProfiler' = 'sienna1','BASCULE' = 'dodgerblue4'))+
  # scale_fill_manual(values=SPN_colors)+
  my_ggplot_theme()+
  xlab("")+
  theme(axis.text.x = element_text(angle = 45,hjust = 1),strip.text.x.top = element_blank())

library(ggh4x)

cosine_samples <-df_all_combs_SPN_signatures %>% 
  separate(col = sample,into = c("spn_id","sample_small"),sep = "_",remove = F) %>% 
  ggplot(aes(x = sample_small, y = cosine, color = caller)) +
  geom_boxplot() +
  facet_nested(
    ~ context + spn,
    scales = "free",
    space = "free_x"
  ) +
  scale_color_manual(
    values = c(
      'SigProfiler' = 'sienna1',
      'BASCULE'     = 'dodgerblue4'
    )
  ) +
  my_ggplot_theme() +
  labs(x = NULL, y = "Cosine similarity") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x.top = element_blank(),panel.spacing.x = unit(0.2, "lines")
  )+
  guides(
    col  = guide_legend(nrow = 2)
  )
cosine_samples
#### geom_line
mut_counts <- df_all_combs_SPN_signatures %>% 
  group_by(sample,context) %>%
  mutate(mean_mut_counts=mean(mutation_count)) %>% 
  ggplot(aes(x=sample,y=log10(mean_mut_counts), group = spn,colour = spn))+
  geom_point()+
  geom_line()+
  my_ggplot_theme()+
  scale_color_manual(values=SPN_colors)+
  facet_wrap(~context)+
  xlab("")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())
## geom_lollipop
mut_counts <- df_all_combs_SPN_signatures %>% 
  group_by(sample,context) %>%
  mutate(mean_mut_counts=mean(mutation_count)) %>% 
  ggplot(aes(x=sample,y=log10(mean_mut_counts), group = spn,colour = spn))+
  geom_segment(
    aes(
      x = sample,
      xend = sample,
      y = 0,
      yend = log10(mean_mut_counts),
      colour = spn
    ),
    linewidth = 0.3,
    linetype = "solid"
  ) +
  geom_point(size = 2)+
  my_ggplot_theme()+
  scale_color_manual(values=SPN_colors)+
  facet_wrap(~context,scales = "free")+
  xlab("")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

mut_counts <- df_all_combs_SPN_signatures %>% 
  group_by(sample,context) %>%
  mutate(mean_mut_counts=mean(mutation_count)) %>% 
  ggplot(aes(x = sample, y="",fill = log10(mean_mut_counts))) +
  geom_tile(color = "white", linewidth = 0.3) +
  # scale_fill_viridis_c(name = "Mutation count") +
  scale_fill_stepsn(
    colors = c("#f7fcf0", "#ccebc5", "#7bccc4", "#2b8cbe", "#084081"),
    name = "Mutation count"
  )+
  labs(x = "", y = "Log10(Mutation counts)") +
  my_ggplot_theme()+
  facet_wrap(~context,scales="free")+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank())

#exp_plot/cosine_samples/mut_counts
mut_counts     <- mut_counts     + theme(plot.margin = margin(0, 0, 0, 0),panel.spacing = unit(0.05, "lines"))
exp_plot       <- exp_plot       + theme(plot.margin = margin(0, 0, 0, 0),panel.spacing = unit(0.05, "lines"))
cosine_samples <- cosine_samples + theme(plot.margin = margin(0, 0, 0, 0),panel.spacing = unit(0.05, "lines"))
abs_counts_exp_plot <- abs_counts_exp_plot + theme(plot.margin = margin(0, 0, 0, 0),panel.spacing = unit(0.05, "lines"))
wrap_plots(list(abs_counts_exp_plot,cosine_samples),
           design = "AAAA\nBBBB\nBBBB\nBBBB",guides="collect")+
  plot_layout(heights = c(1, 1, 0.6, 0.6)) &
  theme(legend.position = "bottom")
  # theme(legend.position = "bottom",plot.margin = margin(0, 0, 0, 0))


### Mean summacontext### Mean summary stats
df_all_combs_SPN_signatures <- df_all_combs_SPN_signatures %>% 
  group_by(sample,purity,caller,context) %>% 
  mutate(mean_cosine=mean(cosine)) %>% 
  mutate(mean_mse=mean(mse)) %>% 
  mutate(tumour_type=case_when(spn=="SPN04" ~ "AML", 
                               spn=="SPN03" ~ "CLL",
                               spn=="SPN07" ~ "GBM",
                               spn=="SPN06" ~ "LUAD",
                               TRUE ~ "COADREAD")) %>% 
  mutate(mutation_class=case_when((mutation_count>=1e4 & type=="INDEL") | (mutation_count>=1e5 & type=="SNV")~"high TMB",
                                  TRUE ~ "low TMB"))
plt_singatures <-df_all_combs_SPN_signatures %>%  ggplot(aes(x = as.factor(purity), y = cosine, fill = caller)) +
  geom_boxplot(position = position_dodge(width = 0.8),
               outlier.alpha = 0.2)+
  facet_wrap(~context)+
  scale_fill_manual(values=c('BASCULE' = 'purple1','SigProfiler' = 'orange2'))+
  my_ggplot_theme()  
plt_singatures_ann <- ggplot(df_all_combs_SPN_signatures, aes(
  x = as.factor(purity),
  y = cosine,
  fill=caller,
  color=caller
)) +
  
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  scale_color_manual(name = "Signature caller", values = c('BASCULE' = 'purple1','SigProfiler' = 'orange2')) +
  scale_fill_manual(name = "Signature caller", values = c('BASCULE' = 'purple1','SigProfiler' = 'orange2')) +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  # Jitter points colored by FGA class
  geom_jitter(aes(y=mean_cosine,x=as.factor(purity),color = mutation_class, fill = mutation_class),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.8) +
  # scale_color_manual(name = "Tumour type", values = c("COADREAD"="#958dbaff","AML"="#5d8964ff",
  #                                                     "CLL"="#f35d74ff","GBM"="#5c9f9bff","LUAD"="#bf624cff")) +
  # scale_fill_manual(name = "Tumour type", values = c("COADREAD"="#958dbaff","AML"="#5d8964ff",
  #                                                    "CLL"="#f35d74ff","GBM"="#5c9f9bff","LUAD"="#bf624cff")) +
  my_ggplot_theme()+
  labs(x="Purity")+
  facet_wrap(~context)

ggplot(df_all_combs_SPN_signatures, aes(
  x = as.factor(purity),
  y = cosine,
  fill=caller,
  color=caller
)) +
  
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  scale_color_manual(name = "Signature caller", values = c('BASCULE' = 'purple1','SigProfiler' = 'orange2')) +
  scale_fill_manual(name = "Signature caller", values = c('BASCULE' = 'purple1','SigProfiler' = 'orange2')) +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  # Jitter points colored by FGA class
  geom_jitter(aes(y=mean_cosine,x=as.factor(purity),color = spn, fill = spn),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.4) +
  scale_color_manual(name = "SPN", values = SPN_colors) +
  scale_fill_manual(name = "SPN", values =SPN_colors)+

  # scale_color_manual(name = "Tumour type", values = c("COADREAD"="#958dbaff","AML"="#5d8964ff",
  #                                                     "CLL"="#f35d74ff","GBM"="#5c9f9bff","LUAD"="#bf624cff")) +
  # scale_fill_manual(name = "Tumour type", values = c("COADREAD"="#958dbaff","AML"="#5d8964ff",
  #                                                    "CLL"="#f35d74ff","GBM"="#5c9f9bff","LUAD"="#bf624cff")) +
  my_ggplot_theme()+
  labs(x="Purity")+
  facet_wrap(~context)


ggplot(df_all_combs_SPN_signatures, aes(
  x = as.factor(purity),
  y = mse,
  fill=caller,
  color=caller
)) +
  
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +
  scale_color_manual(name = "Signature caller", values = c('BASCULE' = 'purple1','SigProfiler' = 'orange2')) +
  scale_fill_manual(name = "Signature caller", values = c('BASCULE' = 'purple1','SigProfiler' = 'orange2')) +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  # Jitter points colored by FGA class
  geom_jitter(aes(y=mean_mse,x=as.factor(purity),color = spn, fill = spn),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              alpha = 0.4) +
  scale_color_manual(name = "SPN", values = SPN_colors) +
  scale_fill_manual(name = "SPN", values =SPN_colors)+
  
  # scale_color_manual(name = "Tumour type", values = c("COADREAD"="#958dbaff","AML"="#5d8964ff",
  #                                                     "CLL"="#f35d74ff","GBM"="#5c9f9bff","LUAD"="#bf624cff")) +
  # scale_fill_manual(name = "Tumour type", values = c("COADREAD"="#958dbaff","AML"="#5d8964ff",
  #                                                    "CLL"="#f35d74ff","GBM"="#5c9f9bff","LUAD"="#bf624cff")) +
  my_ggplot_theme()+
  labs(x="Purity")+
  facet_wrap(~context)
#ggsave("fig4_panelA.pdf",width = 10,height = 4)

