library(ggplot2)
library(dplyr)
library(VIBER)
library(patchwork)
library(gghighlight)
library(optparse)
setwd('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assing_signature/')

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/tumourevo_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN03'),
                    make_option(c("--purity"), type = "double", default = 0.9),
                    make_option(c("--coverage"), type = "integer", default = 50),
                    make_option(c("--cna_caller"), type = "character", default = 'ascat'),
                    make_option(c("--vcf_caller"), type = "character", default = 'mutect2'),
                    make_option(c("--signature_caller"), type = "character", default = 'SigProfiler'),
                    make_option(c("--subclonal_caller"), type = "character", default = 'pyclonevi')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
spn <- opt$spn_id
cov <- opt$coverage
pur <- opt$purity
cna_caller <- opt$cna_caller
vcf_caller <- opt$vcf_caller

sig_dec_tool <- opt$signature_caller
sub_dec_tool <- opt$subclonal_caller

base_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assing_signature/"
data_subclonal <- get_tumourevo_subclonal(spn = spn, 
                                          coverage = cov, 
                                          purity = pur, 
                                          tool = sub_dec_tool,
                                          vcf_caller = vcf_caller,
                                          cna_caller = cna_caller)


if (sub_dec_tool=="viber"){
  fit <- readRDS(data_subclonal$viber_best_st_fit_rds)
  data_table <- bind_cols(cluster = fit$x$cluster.Binomial, fit$data) %>% 
    select(contains("VAF"),cluster) %>% 
    distinct()
} else if (sub_dec_tool=="pyclonevi"){
  fit <- read.table(file = data_subclonal$best_fit_txt,header = T,sep="\t")
  input <- read.table(file = data_subclonal$pyclone_input_all_samples_tsv,header = T,sep="\t")
  data_table <- inner_join(fit,input,by = c("mutation_id","sample_id")) %>% 
    mutate(VAF=alt_counts/(alt_counts+ref_counts)) %>% 
    pivot_wider(
      id_cols = c(mutation_id, cluster_id, driver_label, patient_id),  # identifiers to keep
      names_from = sample_id,                                      # sample names as new columns
      values_from = VAF,                                              # values to spread
      names_prefix = "VAF."                                           # prefix for new columns
    ) %>% 
    mutate(cluster=paste0("C",cluster_id)) %>% 
    select(contains("VAF"),cluster) %>% 
    distinct()
}



stats_cluster <- data_table %>% group_by(cluster) %>% summarize(N = n()) 
samples <- data_table %>% 
  select(contains("VAF")) %>% 
  colnames() %>% 
  str_replace(pattern = "VAF.",replacement = "")

pairs <- combn(samples, 2, simplify = FALSE)

res_assignment_dir <- file.path(base_dir,spn,paste0(cov,"x_",pur,"p_",vcf_caller,"_",cna_caller),paste0(sub_dec_tool,"_",sig_dec_tool))
contexts <- c("ID83","SBS96")
activity_data_list <- list()
for (context in contexts){
  activity_data_list[[context]] <- read.table(file = paste0(res_assignment_dir,"/",
                                                            context,
                                                            "/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"),header = T)
}


df_sbs <-activity_data_list$SBS96 %>% 
  pivot_longer(
               cols = starts_with("SBS"),
               names_to = "Signature",
               values_to = "Value")

df_id <-activity_data_list$ID83 %>% 
  pivot_longer(
    cols = starts_with("ID"),
    names_to = "Signature",
    values_to = "Value")
exp_id <- df_id %>% 
  ggplot(aes(fill=Signature, y=Value, x=as.factor(Samples))) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()+
  scale_fill_manual(values=id_colors)+
  CNAqc:::my_ggplot_theme()+
  ylab("Cluster")+
  xlab("Exposures")

exp_sbs <- df_sbs %>% 
  ggplot(aes(fill=Signature, y=Value, x=as.factor(Samples))) + 
  geom_bar(position="fill", stat="identity")+
  coord_flip()+
  scale_fill_manual(values=sbs_colors)+
  CNAqc:::my_ggplot_theme()+
  ylab("Cluster")+
  xlab("Exposures")

stats_cluster_plt <- stats_cluster %>% ggplot(aes(x=cluster,y=N,fill=cluster))+geom_col()+scale_fill_manual(values=cluster_palette)+
  CNAqc:::my_ggplot_theme()+
  theme(legend.position = "bottom")+
  ylab("Number of mutations")+
  coord_flip()
cluster_plots_pair <- list()

for (i in 1:length(pairs)){
  s1 <- paste0("VAF.",pairs[[i]][1])
  s2 <- paste0("VAF.",pairs[[i]][2])
  
  cluster_plots_pair[[i]] = data_table %>%
    ggplot(aes(x =.data[[s1]], y = .data[[s2]],color=cluster))+
    geom_point()+
    scale_color_manual(values = cluster_palette)+
    xlim(0,1)+
    ylim(0,1)+
    xlab(pairs[[i]][1]) +
    ylab(pairs[[i]][2])+
    CNAqc:::my_ggplot_theme()+
    theme(legend.position = "bottom")
}

plt_marg <- wrap_plots(cluster_plots_pair)

vaf_histo <- list()
for (s in samples){
  s1 <- paste0("VAF.",s)
  #   s2 <- paste0("VAF.",pairs[[i]][2])
  vaf_histo[[s]] = data_table %>%
    ggplot(aes(x = .data[[s1]], fill = cluster)) +
    xlim(0.01,1)+
    geom_histogram(binwidth = 0.01, na.rm = TRUE) +
    # gghighlight(cluster == c,
    #             keep_scales = TRUE,unhighlighted_params = list( fill = alpha("grey", 0.2)))+
    scale_fill_manual(values=cluster_palette)+
    CNAqc:::my_ggplot_theme() +
    theme(legend.position = "bottom")
  
}
vaf_histo_plot <- wrap_plots(vaf_histo,nrow = 1)

wrap_plots(list(exp_sbs,exp_id,stats_cluster_plt,plt_marg,vaf_histo_plot),guides="collect",design = "ACDDD\nBCDDD\nEEEEE") & theme(legend.position = "bottom")


# clusters <- data_table %>% pull(cluster) %>% unique()
# cluster_palette <- hue_pal()(length(clusters))
# names(cluster_palette) <- clusters
# final_clusters_plot <- list()
# pairs <- combn(samples, 2, simplify = FALSE)
# for (c in clusters){
#   exp <- df %>% 
#     filter(Samples==c) %>% 
#     filter(Value!=0) %>% 
#     ggplot(aes(fill=Signature, y=Value, x=as.factor(Samples))) + 
#     geom_bar(position="fill", stat="identity")+
#     coord_flip()+
#     scale_fill_manual(values=sbs_colors)+
#     theme_bw()
#   cluster_plots_pair <- list()
#   vaf_histo <- list()
#   for (i in 1:length(pairs)){
#     s1 <- paste0("VAF.",pairs[[i]][1])
#     s2 <- paste0("VAF.",pairs[[i]][2])
# 
#     cluster_plots_pair[[i]] = data_table %>%
#       ggplot(aes(x =.data[[s1]], y = .data[[s2]],color=cluster))+
#       geom_point()+
#       scale_color_manual(values = cluster_palette)+
#       xlab(s1) +
#       ylab(s2)+
#       CNAqc:::my_ggplot_theme()+
#       theme(legend.position = "bottom")
#   }
#   for (s in samples){
#     s1 <- paste0("VAF.",s)
#     #   s2 <- paste0("VAF.",pairs[[i]][2])
#     vaf_histo[[s]] = data_table %>%
#         ggplot(aes(x = .data[[s1]], fill = cluster)) +
#         geom_histogram(binwidth = 0.01, na.rm = TRUE) +
#         gghighlight(cluster == c,
#                     keep_scales = TRUE,unhighlighted_params = list( fill = alpha("grey", 0.2)))+
#         CNAqc:::my_ggplot_theme() +
#         theme(legend.position = "bottom")
#     
#   }
#   cluster_plot <- wrap_plots(cluster_plots_pair,nrow = 3)
#   vaf_histo_plot <- wrap_plots(vaf_histo,nrow = 1)
#   final_clusters_plot[[c]] <- wrap_plots(exp,vaf_histo_plot,design = "AAAAAA
#                                                  BBBBBB")
# 
# }
# 
# wrap_plots(final_clusters_plot,ncol = 2,guides="collect")
# 
