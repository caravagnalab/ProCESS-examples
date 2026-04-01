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
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN05", "SPN06",'SPN07',"SPN02")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION/"


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


cna_callers <- c("sequenza","ascat")
variant_callers <- c("mutect2","strelka")
contexts <- c("SBS96","ID83")
params_grid = expand.grid(coverages, purities,cna_callers,variant_callers,contexts)
colnames(params_grid) = c("coverage", "purity","cna_caller","snv_caller","context")

df_all_combs_SPN_signatures <- lapply(SPNS,function(spn){
  validation_dir_spn <- paste0(validation_dir,spn,"/signature/")
  all_combs <- list()
  all_metrics <- list()
  all_cosine  <- list()
  all_metrics_sample <-list()
  for (i in 1:nrow(params_grid)){
    coverage <- params_grid$coverage[i]
    purity <- params_grid$purity[i]
    comb <- paste0(coverage,"x_",purity)
    ctx <- params_grid$context[i]
    cna_caller <- params_grid$cna_caller[i]
    variant_caller <- params_grid$snv_caller[i]
    comb_caller <- paste0(variant_caller,"_",cna_caller)
    print(comb_caller)
    rds_metrics <- file.path(validation_dir_spn,comb,comb_caller,paste0("metrics_",ctx,"_spn.rds"))
    rds_metrics_sample <- file.path(validation_dir_spn,comb,comb_caller,paste0("metrics_",ctx,"_sample.rds"))
    rds_cosine  <- file.path(validation_dir_spn,comb,comb_caller,paste0("cosine_mse_",ctx,".rds"))
    if (file.exists(rds_metrics) & file.exists(rds_metrics_sample) & file.exists(rds_cosine)){
      all_metrics[[i]] <- readRDS(rds_metrics) %>% 
        dplyr::mutate(context=ctx,
                      purity=purity,
                      cna_caller=cna_caller,
                      snv_caller=variant_caller,
                      coverage=coverage)
      all_cosine[[i]] <- readRDS(rds_cosine) %>% 
        dplyr::mutate(context=ctx,
                      purity=purity,
                      cna_caller=cna_caller,
                      snv_caller=variant_caller,
                      coverage=coverage)
      all_metrics_sample[[i]] <-readRDS(rds_metrics_sample) %>% 
        dplyr::mutate(context=ctx,
                      purity=purity,
                      cna_caller=cna_caller,
                      snv_caller=variant_caller,
                      coverage=coverage) %>% 
        unique()
    } else{
      print(paste0("Missing data for ",comb_caller))
    }

  }
  cosine_metrics <- all_cosine %>% bind_rows()
  metrics_ss <- all_metrics_sample %>% bind_rows()
  final_df <- inner_join(cosine_metrics,metrics_ss)
}) %>% bind_rows()


df_all_combs_SPN_signatures <-df_all_combs_SPN_signatures %>% 
  mutate(spn=case_when(spn=="SPN07_last"~"SPN07",
                       TRUE ~spn)) %>% 
  mutate(sample = str_replace_all(sample, "last_", ""))




indir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assing_signature/"
jaccard_similarity <- function(x, y) {
  length(intersect(x, y)) / length(union(x, y))
}
subclonal_arch <- lapply(SPNS, function(spn) {
  
  all_combs <- list()
  all_combs_final <- list()
  for (pur in purities) {
    for (cov in coverages) {
      
      all_ctx <- list()
      
      for (ctx in c("SBS", "ID")) {
        
        comb <- paste0(cov, "x_", pur, "p")
        samples <- get_sample_names(spn = spn)
        all_samples <- list()
        for (s in samples){
          print(s)
          df <- readRDS(
            file.path(indir, spn, "process_univariate", comb,
                      paste0(s,"_exposure_", ctx, ".rds"))
          ) %>%
            mutate(spn = spn) %>% 
            mutate(context=ctx) %>% 
            mutate(cluster_id_process=case_when(is.na(cluster_id_process)~"Subclonal",
                   TRUE ~ cluster_id_process)) %>% 
            filter(cluster_id_process%in%c("Clonal","Subclonal")) %>% 
            filter(nmuts_cause>=20)
          c_df <- df %>% 
            select(cluster_id_process,causes) %>%
            distinct(cluster_id_process, causes) %>%
            group_by(cluster_id_process) %>%
            summarise(causes = list(causes)) %>%
            deframe()
          jaccard_score <-jaccard_similarity(c_df$Clonal,c_df$Subclonal)
          ex_df <- df %>%
            select(cluster_id_process,causes,exposure) %>%
            complete(cluster_id_process, causes, fill = list(exposure = 0)) %>%
            arrange(causes) %>%
            group_by(cluster_id_process) %>%
            summarise(causes = list(exposure)) %>%
            deframe()
          cosine_similarity <- as.numeric(lsa::cosine(x = ex_df$Clonal,y=ex_df$Subclonal))
          df <- df %>% 
            mutate(cosine_similary_clusters=cosine_similarity) %>%
            mutate(jaccard_similary_clusters=jaccard_score)
          all_samples[[s]] <- df
        }
        all_ctx[[ctx]] <- bind_rows(all_samples)
      }
      
      all_combs[[comb]] <- bind_rows(all_ctx)
    }
  }
  all_combs_final <- bind_rows(all_combs)
  return(all_combs_final)   # ← THIS WAS MISSING
}) %>% bind_rows()

subclonal_arch <- subclonal_arch %>% 
  group_by(sample_id,context) %>% 
  mutate(median_jaccard_similary_clusters=median(jaccard_similary_clusters)) %>% 
  mutate(sample_class=case_when(median_jaccard_similary_clusters>=0.9~"low complexity",
                                TRUE~"high complexity")) %>% 
  mutate(context=case_when(context=="SBS"~"SBS96",
                           context=="ID"~"ID83")) %>% 
  mutate(mean_cosine_similary_clusters=mean(cosine_similary_clusters)) %>%
  # mutate(sample_class=case_when(median_jaccard_similary_clusters>=0.9 & mean_cosine_similary_clusters>=median(.data$mean_cosine_similary_clusters) ~ "low complexity",
  #                               median_jaccard_similary_clusters>=0.9 & mean_cosine_similary_clusters<median(.data$mean_cosine_similary_clusters) ~ "high complexity",
  #                               median_jaccard_similary_clusters<0.9 ~ "high complexity")) %>% 
  dplyr::rename(sample=sample_id) %>% 
  ungroup()

#saveRDS(object = subclonal_arch,file="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure5/signatures_cohort/sample_classification_clonal_heterogeneity.rds")
#subclonal_arch <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples-my/figures/figure5/signatures_cohort/sample_classification_clonal_heterogeneity.rds")
df_all_combs_SPN_signatures <- df_all_combs_SPN_signatures %>% 
  left_join(subclonal_arch %>% select(sample,sample_class,context,median_jaccard_similary_clusters) %>% distinct(),by=c("sample","context")) 

# table <- df_all_combs_SPN_signatures %>% 
#   group_by(sample_class, caller,context ) %>% 
#   summarize(mean_cs = mean(cosine),
#             median_cs = median(cosine),
#             sd_cs = sd(cosine))
# write.table(x = table, file = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure5/panelB_cohort_signature.tsv', quote = F, sep = '\t', row.names = F)

my_comparisons <- list(c("High\nComplexity", "Low\nComplexity"))
plt_signatures <- df_all_combs_SPN_signatures %>%
  mutate(comb_te=paste0(cna_caller,"_",snv_caller)) %>% 
  mutate(sample_class = ifelse(sample_class == "high complexity", "High\nComplexity", "Low\nComplexity")) %>% 
  ggplot(aes(x=sample_class,y=cosine))+
  geom_boxplot(aes(color=caller, fill = caller),outliers = T,width = 0.5,alpha=0.3)+

  scale_color_manual('Caller', values = c('SigProfiler'='sienna1',
                                'BASCULE'='dodgerblue4'))+
  scale_fill_manual('Caller',values = c('SigProfiler'='sienna1',
                                'BASCULE'='dodgerblue4'))+
  facet_grid(snv_caller~context)+
  stat_compare_means(comparisons = my_comparisons,label = "p.signif",vjust = 0.5)+
  my_ggplot_theme() +
  xlab('Sample Class') + 
  ylab('Exposure Accuracy\nCosine Similarity')

df_all_combs_SPN_signatures %>%
  mutate(comb_te=paste0(cna_caller,"_",snv_caller)) %>% 
  mutate(sample_class = ifelse(sample_class == "high complexity", "High\nComplexity", "Low\nComplexity")) %>% 
  ggplot(aes(x=spn,y=cosine))+
  geom_boxplot(aes(color=caller, fill = caller),outliers = T,width = 0.5,alpha=0.3)+
  
  scale_color_manual('Caller', values = c('SigProfiler'='sienna1',
                                          'BASCULE'='dodgerblue4'))+
  scale_fill_manual('Caller',values = c('SigProfiler'='sienna1',
                                        'BASCULE'='dodgerblue4'))+
  facet_grid(sample_class~context)+my_ggplot_theme()

# 
# plt_spn03_1 <-df_all_combs_SPN_signatures %>% 
#   filter(spn=="SPN03") %>% 
#   group_by(sample,context) %>% 
#   mutate(mean_cosine=mean(cosine)) %>% 
#   mutate(upper_cosine=max(cosine)) %>% 
#   mutate(lower_cosine=min(cosine)) %>% 
#   ggplot(aes(x=sample,y = mean_cosine))+
#     geom_pointrange(aes(ymin = lower_cosine, 
#                         ymax = upper_cosine,color =  context),
#                     position = position_dodge(0.3))+
#   xlab("")+
#   ylab("Cosine similarity")+
#   my_ggplot_theme()
# 
# plt_spn03_1.1 <-df_all_combs_SPN_signatures %>% 
#   filter(spn=="SPN03") %>% 
#   ggplot(aes(x=sample,y = cosine))+
#   geom_boxplot(aes(color=context),outliers = F,width = 0.5)+
#   xlab("")+
#   ylab("Cosine similarity")+
#   my_ggplot_theme()
# 
# plt_spn03_2 <- subclonal_arch %>% 
#   filter(spn=="SPN03") %>% 
#   select(sample,median_jaccard_similary_clusters,sample_class,context) %>% 
#   distinct() %>% 
#   ggplot(aes(x=sample,y=median_jaccard_similary_clusters,fill=context))+
#     geom_col(width = 0.4,aes(fill=context),alpha=0.4,position = position_dodge())+
#     # scale_fill_manual(values=list("high complexity"="firebrick4","low complexity"="lightpink1"))+
#   my_ggplot_theme()+
#   xlab("")+
#   ylab("Jaccard Index")+
#   theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
#         axis.minor.ticks=element_blank(),axis.ticks.x = element_blank())
# 
# plt_spn03_3 <- subclonal_arch %>% 
#   filter(spn=="SPN03") %>% 
#   ggplot(aes(x=sample,y=cosine_similary_clusters,fill=context))+
#   geom_boxplot(width = 0.4,aes(fill=context),alpha=0.4)+
#   scale_fill_manual(values=list("high complexity"="firebrick4","low complexity"="lightpink1"))+
#   my_ggplot_theme()+
#   xlab("")+
#   ylab("Clonal Heterogeneity")+
#   theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
#         axis.minor.ticks=element_blank(),axis.ticks.x = element_blank())+
#   ggtitle(label = "SPN03")
# 
# 
# plt_spn03_4 <- subclonal_arch %>% 
#   filter(spn == "SPN03") %>% 
#   ggplot(aes(y = mean_cosine_similary_clusters,
#              x = sample)) +
#   geom_linerange(aes(ymin = 0,
#                    ymax = mean_cosine_similary_clusters,
#                    color = context),
#                  # height = 0,
#                position = position_dodge(width = .3)) +
#   geom_point(aes(size = median_jaccard_similary_clusters,
#                  color = context),
#              position = position_dodge(width = .3)) +
#   ylim(0,0.75)+
#   my_ggplot_theme()+
#   theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
#         axis.minor.ticks=element_blank(),axis.ticks.x = element_blank())
#     
# wrap_plots(list(plt_spn03_4,plt_spn03_1.1),design = "AAA\nBBB\nBBB",guides = "collect") &
#   theme(legend.position = "bottom")
# wrap_plots(list(plt_spn03_3,plt_spn03_2,plt_spn03_1),design = "AAA\nBBB\nCCC\nCCC",guides = "collect") &
#   theme(legend.position = "bottom")
purity_effect <-df_all_combs_SPN_signatures %>%
  mutate(comb_te=paste0(cna_caller,"_",snv_caller)) %>% 
  mutate(sample_class = ifelse(sample_class == "high complexity", "High\nComplexity", "Low\nComplexity"))  %>%

  ggplot(aes(
    x = sample_class,
    y=cosine,
    fill = as.factor(purity)
  )) +
  geom_boxplot(outlier.alpha = 0.4)+
  facet_grid(snv_caller~context)+
  scale_fill_manual("Purity",values = purity_colors) +
  my_ggplot_theme() +
  labs(
    y = "Cosine similarity",
    x = "Purity"
  )
# 
# 
# 
# 
# wrap_plots(list(ggplot(),plt_caller,purity_effect,plt_spn03_4,plt_spn03_1.1),design = "AAABB\nAAACC\nDD###\nEE###",
#            guides="collect") & theme(legend.position = "bottom")



# plt_signatures <- df_all_combs_SPN_signatures %>% 
#   group_by(context) %>% 
#   mutate(mean_cosine=mean(cosine)) %>% 
#   mutate(purity = as.character(purity)) %>% 
#   arrange(desc(muts_per_signature)) %>% 
#   mutate(log10_muts_count = log10(muts_per_signature)) %>% 
#   ggplot(aes(x = log10_muts_count,
#              y = cosine,
#              group = interaction(caller, context))) +
#   geom_point(alpha=0.4,aes(color=spn)) +
#   scale_color_manual(values = SPN_colors)+
#   ggnewscale::new_scale_color() +
#   geom_smooth(span = 0.1,se = F,aes(fill=caller,color=caller),alpha = 0.1,method = "loess") +
#   scale_color_manual(values = c('SigProfiler'='sienna1',
#                                 'BASCULE'='dodgerblue4'))+
#   geom_hline(aes(yintercept = mean_cosine),linetype = "dashed",color="grey40")+
#   facet_wrap(~context,scales = "free",nrow = 2,strip.position = "right")+
#   my_ggplot_theme()+
#   xlab("Log10(Median Mutation count per signature)")+
#   ylab("Cosine Similarity")

