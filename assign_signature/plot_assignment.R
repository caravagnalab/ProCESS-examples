rm(list=ls())
library(ggplot2)
library(dplyr)
library(ctree)
library(VIBER)
library(patchwork)
library(tidyr)
library(gghighlight)
library(optparse)
library(stringr)
library(forcats)
library(scales)
library(ggrepel)
base_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assing_signature/"
setwd(base_dir)

source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/tumourevo_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN06'),
                    make_option(c("--purity"), type = "double", default = 0.9),
                    make_option(c("--coverage"), type = "integer", default = 100),
                    make_option(c("--cna_caller"), type = "character", default = 'ascat'),
                    make_option(c("--vcf_caller"), type = "character", default = 'mutect2'),
                    make_option(c("--signature_caller"), type = "character", default = 'SigProfiler'),
                    make_option(c("--subclonal_caller"), type = "character", default = 'pyclonevi')
)
spns <- c("SPN06")
coverages <- c(100)
purities <- c(0.9)
sig_dec_tools <- c("SigProfiler")
sub_dec_tools <- c("viber_heuristic")
params_grid = expand.grid(spns, coverages,purities,sig_dec_tools,sub_dec_tools)
colnames(params_grid) = c("spn","coverage", "purity","signature_caller","subclonal_caller")
for (i in 1:nrow(params_grid)){
  spn <- params_grid[i,]$spn
  cov <- params_grid[i,]$coverage
  pur <- params_grid[i,]$purity
  cna_caller <- "ascat"
  vcf_caller <- "mutect2"

  sig_dec_tool <- params_grid[i,]$signature_caller
  sub_dec_tool <- params_grid[i,]$subclonal_caller
  if (grepl(pattern = "viber",sub_dec_tool)){
    data_subclonal <- get_tumourevo_subclonal(spn = spn, 
                                              coverage = cov, 
                                              purity = pur, 
                                              tool = "viber",
                                              vcf_caller = vcf_caller,
                                              cna_caller = cna_caller)  
  } else {
    data_subclonal <- get_tumourevo_subclonal(spn = spn, 
                                              coverage = cov, 
                                              purity = pur, 
                                              tool = sub_dec_tool,
                                              vcf_caller = vcf_caller,
                                              cna_caller = cna_caller)  
  }
  
  
  data_subclonal_ctree <- get_tumourevo_subclonal(spn = spn, 
                                                  coverage = cov, 
                                                  purity = pur, 
                                                  tool = "ctree",
                                                  vcf_caller = vcf_caller,
                                                  cna_caller = cna_caller,sample = spn)
  
  if (sub_dec_tool=="viber"){
    fit <- readRDS(data_subclonal$viber_best_st_fit_rds)
    data_table <- bind_cols(cluster = fit$x$cluster.Binomial, fit$data) %>% 
      select(contains("VAF"),cluster,driver,gene) %>% 
      distinct()
    ctree_ccf <- readRDS(data_subclonal_ctree[["ctree_VIBER_rds"]])[[1]]$CCF %>% 
      rename("Driver cluster"=is.driver) %>% 
      rename("Clonal cluster"=is.clonal)
    ctree_ccf_long <- pivot_longer(ctree_ccf, cols = c("Driver cluster","Clonal cluster"),
                                   names_to = "feature", values_to = "value")
    ctree_drivers <- readRDS(data_subclonal_ctree[["ctree_VIBER_rds"]])[[1]]$drivers
    ctree_plot <-  plot(readRDS(data_subclonal_ctree[["ctree_VIBER_rds"]])[[1]])
    
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
      select(contains("VAF"),cluster,driver_label) %>% 
      rename(variantID=driver_label) %>% 
      distinct()
    ctree_ccf <- readRDS(data_subclonal_ctree[["ctree_pyclonevi_rds"]])[[1]]$CCF %>% 
      rename("Driver cluster"=is.driver) %>% 
      rename("Clonal cluster"=is.clonal)
    ctree_ccf_long <- pivot_longer(ctree_ccf, cols = c("Driver cluster","Clonal cluster"),
                                   names_to = "feature", values_to = "value")
    ctree_drivers <- readRDS(data_subclonal_ctree[["ctree_pyclonevi_rds"]])[[1]]$drivers %>% 
      inner_join(y = data_table,by = c("variantID","cluster"))
    ctree_plot <-  plot(readRDS(data_subclonal_ctree[["ctree_pyclonevi_rds"]])[[1]])
  } else if (sub_dec_tool=="viber_heuristic"){
    fit <- readRDS(data_subclonal$viber_best_st_heuristic_fit_rds)
    fit[["theta_k"]]<-round(fit[["theta_k"]],2)
    data_table <- bind_cols(cluster = fit$x$cluster.Binomial, fit$data) %>% 
      select(contains("VAF"),cluster) %>% 
      distinct()
    trees <- get_clone_trees(fit)
    ctree_ccf <- trees[[1]]$CCF %>%
      rename("Driver cluster"=is.driver) %>%
      rename("Clonal cluster"=is.clonal)
    ctree_ccf_long <- pivot_longer(ctree_ccf, cols = c("Driver cluster","Clonal cluster"),
                                   names_to = "feature", values_to = "value")
    ctree_drivers <- trees[[1]]$drivers
    ctree_plot <-  plot(trees[[1]])
    ctree_plot <- ctree_plot +theme(legend.position = "none")
    
  }
  
  
  
  stats_cluster <- data_table %>% filter(!is.na(cluster)) %>% group_by(cluster) %>% summarize(N = n()) 
  samples <- data_table %>% 
    select(contains("VAF")) %>% 
    colnames() %>% 
    str_replace(pattern = "VAF.",replacement = "") %>% 
    sort()
  
  contexts <- c("ID83","SBS96")
  res_assignment_dir <- file.path(base_dir,spn,paste0(cov,"x_",pur,"p_",vcf_caller,"_",cna_caller),paste0(sub_dec_tool,"_",sig_dec_tool))
  if (sig_dec_tool=="SigProfiler"){
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
        values_to = "Exposure")
    
    df_id <-activity_data_list$ID83 %>% 
      pivot_longer(
        cols = starts_with("ID"),
        names_to = "Signature",
        values_to = "Exposure")
    
  } else {
    activity_data_list <- list()
    for (context in contexts){
      activity_data_list[[context]] <-readRDS(file = paste0(res_assignment_dir,"/",
                                                            context,
                                                            "/bascule_fit.rds"))
    }
    df_sbs <-activity_data_list$SBS96$nmf$SBS$exposure %>% 
      rename(Signature=sigs) %>% 
      rename(Exposure=value) %>% 
      rename(Samples=samples)
    df_id <-activity_data_list$ID83$nmf$ID$exposure %>% 
      rename(Signature=sigs) %>% 
      rename(Exposure=value) %>% 
      rename(Samples=samples)
  }
  
  cluster_order <- stats_cluster %>% arrange((N)) %>% pull(cluster)
  exp_id <- df_id %>% 
    filter(Samples!="NotAssigned") %>% 
    mutate(Samples = factor(Samples, levels = cluster_order)) %>% 
    ggplot(aes(fill=Signature, y=Exposure, x=as.factor(Samples))) + 
    geom_bar(position="fill", stat="identity")+
    coord_flip()+
    scale_fill_manual(values=id_colors)+
    # CNAqc:::my_ggplot_theme()+
    theme_bw()+
    theme(legend.position = "bottom")+
    ylab("Cluster")+
    xlab("Exposures")+
    ggtitle(label = "ID83 Exposures")
  
  
  exp_sbs <- df_sbs %>% 
    filter(Samples!="NotAssigned") %>% 
    mutate(Samples = factor(Samples, levels = cluster_order)) %>% 
    ggplot(aes(fill=Signature, y=Exposure, x=as.factor(Samples))) + 
    geom_bar(position="fill", stat="identity")+
    coord_flip()+
    scale_fill_manual(values=sbs_colors)+
    theme_bw()+
    theme(legend.position = "bottom")+
    ylab("Cluster")+
    xlab("Exposures")+
    ggtitle(label = "SBS96 Exposures")
  clusters <- data_table %>% pull(cluster) %>% unique()
  cluster_palette <- hue_pal()(length(clusters))
  names(cluster_palette) <- clusters
  
  stats_cluster_plt <- stats_cluster %>%
    mutate(cluster = fct_reorder(cluster, N)) %>% 
    ggplot(aes(x=cluster,y=N,fill=cluster))+geom_col()+scale_fill_manual(values=cluster_palette)+
    theme_bw()+
    theme(legend.position = "bottom")+
    ylab("Number of mutations")+
    xlab("Cluster")+
    coord_flip()
  
  # plot
  cluster_features_plt <- ctree_ccf_long %>% 
    full_join(ctree_drivers,by = c("cluster"),relationship = "many-to-many") %>% 
    group_by(cluster,feature) %>% 
    mutate(cluster = factor(cluster, levels = cluster_order)) %>% 
    mutate(driver = paste0(variantID, collapse = "\n"))  %>% 
    select(cluster,driver,feature,value,nMuts) %>% 
    unique() %>% 
    # arrange(desc(nMuts)) %>% 
    mutate(cluster = fct_reorder(cluster, nMuts)) %>% 
    ggplot(aes(x = feature, y = cluster, fill = value)) +
    geom_tile(color = "white") +
    geom_text(aes(label = driver), color = "white", size = 2)+
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "grey")) +
    theme_bw()+
    theme(legend.position = "none")+
    labs(x = "", y = "", fill = "State")
  
  
  pairs <- combn(samples, 2, simplify = FALSE)
  cluster_plots_pair <- list()
  
  for (i in 1:length(pairs)){
    s1 <- paste0("VAF.",pairs[[i]][1])
    s2 <- paste0("VAF.",pairs[[i]][2])
    
    # cluster_plots_pair[[i]] = 
    p=data_table %>%
      ggplot(aes(x =.data[[s1]], y = .data[[s2]],color=cluster,alpha=0.4))+
      geom_point()+
      scale_color_manual(values = cluster_palette)+
      xlim(0,1)+
      ylim(0,1)+
      xlab(pairs[[i]][1]) +
      ylab(pairs[[i]][2])+
      theme_bw()+
      theme(legend.position = "none")
    cluster_plots_pair[[i]] =p+
      geom_label_repel(
        data = ctree_drivers,
        aes(
          x = .data[[s1]],
          y = .data[[s2]],
          label = variantID,colour = cluster
        ),
        inherit.aes = FALSE,
        # color = "red",
        fontface = "bold",
        size = 3,
        min.segment.length = 0,
        box.padding = 1)+
      theme(legend.position = "none")
  }
  ncol=3
  nrow=ceiling(length(pairs)/ncol)
  plt_marg <- wrap_plots(cluster_plots_pair,ncol = ncol)
  
  # vaf_histo <- list()
  # for (s in samples){
  #   s1 <- paste0("VAF.",s)
  #   vaf_histo[[s]] <-data_table %>%
  #     ggplot(aes(x = .data[[s1]], fill = cluster)) +
  #     xlim(0,1)+
  #     xlab(s) +
  #     geom_histogram(binwidth = 0.01, na.rm = TRUE) +
  #     # gghighlight(cluster == c,
  #     #             keep_scales = TRUE,unhighlighted_params = list( fill = alpha("grey", 0.2)))+
  #     scale_fill_manual(values=cluster_palette)+
  #     CNAqc:::my_ggplot_theme() +
  #     theme(legend.position = "none")
  # 
  # 
  # 
  #   
  # }
  # vaf_histo_plot <- wrap_plots(vaf_histo,nrow = 1)
  
  final_plot <- wrap_plots(
    list(exp_sbs, exp_id, stats_cluster_plt, cluster_features_plt, ctree_plot, plt_marg),
    design = "ACDE\nBCDE\nFFFF\nFFFF\nFFFF",
    guides = "collect"
  ) +
    plot_annotation(
      title = paste0("Clone-specific mutational signatures for ", spn),
      subtitle = paste0("Clones identified with ", sub_dec_tool,
                        " and mutational signatures defined with ", sig_dec_tool),
      caption = paste0("Coverage: ", cov, ", Purity: ", pur,
                       ", Mutation caller: ", vcf_caller,
                       ", CNA caller: ", cna_caller),
      theme = theme(legend.position = "bottom")
    )
  width = 4*ncol
  height = 2+(4*nrow)
  # max_height =12/length(pairs)/ncol
  
  ggsave(filename = file.path(base_dir,spn,paste0(cov,"x_",pur,"p_",vcf_caller,"_",cna_caller),paste0(sub_dec_tool,"_",sig_dec_tool),"report.pdf"),
         plot = final_plot,height = height,width = width,device = "pdf",dpi = 300)
  
}

# opt_parser <- OptionParser(option_list = option_list)
# opt <- parse_args(opt_parser)
# spn <- opt$spn_id
# cov <- opt$coverage
# pur <- opt$purity
# cna_caller <- opt$cna_caller
# vcf_caller <- opt$vcf_caller
# 
# sig_dec_tool <- opt$signature_caller
# sub_dec_tool <- opt$subclonal_caller

