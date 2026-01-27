rm(list=ls())
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggbreak)
library(tidyr)
library(ggh4x)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure2")
mut_counts_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/mutations_counts"


spns <- c("SPN01","SPN02","SPN03","SPN04","SPN05","SPN06","SPN07")
n_samples <- c(3,2,4,2,3,5,5)
coverages <- c(50,100,150,200)
lots_groups <- c(10,20,30,40)
purities <- c(0.3,0.6,0.9)
params_grid = expand.grid(spns, coverages,purities)
colnames(params_grid) = c("spn","coverage", "purity")

muts_counts_stats <-list()
real_muts_counts_stats <- list()
fastq_sizes <- list()
for (i in 1:nrow(params_grid)){
  spn <- params_grid[i,]$spn
  cov <- params_grid[i,]$coverage
  pur <- params_grid[i,]$purity
  muts_counts_stats[[i]] <- readRDS(file = file.path(mut_counts_dir,spn,paste0("mutations_counts_",cov,"x_",pur,"p.rds"))) %>% 
    mutate(spn=spn)
}
stats <- do.call("rbind",muts_counts_stats)


# for (spn in spns){
#   real_muts_counts_stats[[spn]] <- readRDS(file = file.path(mut_counts_dir,spn,"real_mutations_counts.rds")) %>% 
#     mutate(spn=spn) %>% 
#     mutate(type=case_when(type=="indel"~"INDEL",
#                           TRUE~type)) %>% 
#     rename(process_muts=n) %>% 
#     select(type,process_muts,spn)
# }

# real_stats <- do.call("rbind",real_muts_counts_stats)
# # stats <- stats %>% 
# #   group_by(type,coverage,purity,spn) %>% 
# #   mutate(mean_mutation_count=mean(mutation_count))
# stats_total <- left_join(x = stats,y = real_stats) %>% 
#   mutate(percentage_sequenced=(mutation_count/process_muts)*100) %>% 
#   group_by(type,coverage,purity) %>% 
#   mutate(mean_percentage_sequenced=mean(percentage_sequenced)) %>% 
#   mutate(mean_mutation_counts=mean(mutation_count)) %>% 
#   select(type,coverage,purity,mean_percentage_sequenced,mean_mutation_counts) %>% 
#   unique()



# Compute scale factor
# stats_total_indel <- stats_total %>% 
#   filter(type=="INDEL")
# scale_factor_indel <- max(stats_total_indel$mean_mutation_counts) / max(stats_total_indel$mean_percentage_sequenced)
# 
# plt_muts_indel <- stats_total_indel %>%
#   ggplot(aes(x = coverage, color = as.factor(purity), group = purity)) +
#   # Main variable (solid line)
#   geom_point(aes(y = mean_mutation_counts)) +
#   geom_line(aes(y = mean_mutation_counts)) +
#   
#   # Secondary variable (dashed line)
#   geom_point(aes(y = mean_percentage_sequenced * scale_factor_indel), shape = 1) +
#   geom_line(aes(y = mean_percentage_sequenced * scale_factor_indel), linetype = "dashed") +
#   
#   # Colors and axes
#   scale_color_manual(values = purity_colors) +
#   scale_y_continuous(
#     name = "Mean INDELs counts",
#     sec.axis = sec_axis(~ . / scale_factor_indel, name = "Mean % sequenced INDELs")
#   ) +
#   
#   # facet_wrap(~ type, scales = "free_y") +
#   labs(
#     x = "Coverage",
#     color = "Purity"
#   ) +
#   theme_bw()
# 
# stats_total_snv <- stats_total %>% 
#   filter(type=="SNV")
# scale_factor_snv <- max(stats_total_snv$mean_mutation_counts) / max(stats_total_snv$mean_percentage_sequenced)
# 
# plt_muts_snv <-stats_total_snv %>%
#   ggplot(aes(x = coverage, color = as.factor(purity), group = purity)) +
#   # Main variable (solid line)
#   geom_point(aes(y = mean_mutation_counts)) +
#   geom_line(aes(y = mean_mutation_counts)) +
#   
#   # Secondary variable (dashed line)
#   geom_point(aes(y = mean_percentage_sequenced * scale_factor_snv), shape = 1) +
#   geom_line(aes(y = mean_percentage_sequenced * scale_factor_snv), linetype = "dashed") +
#   
#   # Colors and axes
#   scale_color_manual(values = purity_colors) +
#   scale_y_continuous(
#     name = "Mean SNVs counts",
#     sec.axis = sec_axis(~ . / scale_factor_snv, name = "Mean % sequenced SNVs")
#   ) +
#   # 
#   # facet_wrap(~ type, scales = "free_y") +
#   labs(
#     x = "Coverage",
#     color = "Purity"
#   ) +
#   theme_bw()
# plt_muts <- wrap_plots(list(plt_muts_snv,plt_muts_indel),ncol = 1,guides = "collect") & theme(legend.position = "bottom")
# real_stats_mean <- real_stats %>% 
#   # group_by(type) %>% 
#   # summarise(mean_count=mean(n)) %>% 
#   summarise(sum_count=sum(n)) %>% 
#   mutate(purity="real")
# real_stats_mean <- rbind(real_stats_mean,real_stats_mean,real_stats_mean,real_stats_mean) %>% 
#   mutate(coverage=c(rep(50,2),rep(100,2),rep(150,2),rep(200,2)))
# 
stats <- stats %>%
  group_by(type,coverage,purity,spn) %>%
  mutate(mean_mutation_count=mean(mutation_count))

# stats_mean <- stats %>%
#   mutate(purity=as.factor(purity)) %>%
#   mutate(log_mut_counts=log10(mutation_count)) %>%
#   group_by(type,coverage,purity) %>%
#   # mutate(mean_count=mean(mutation_count)) %>%
#   mutate(sum_count=sum(mutation_count)) %>%
#   select(type,coverage,purity,sum_count) %>%
#   distinct()
# 
# 
# stats %>% 
#   # group_by(spn, type, coverage, purity) %>% 
#   # mutate(median_mutation_count = median(mutation_count)) %>% 
#   ggplot(aes(x = spn, y = log10(mutation_count), fill = as.factor(coverage))) +
#   geom_boxplot() +
#   ylab("Mutation count")+
#   xlab("SPN")+
#   labs(fill="coverage")+
#   theme(legend.title = element_text("coverage"))+
#   scale_fill_manual(values = coverage_colors)+
#   facet_wrap(~ type,scales = "free",ncol = 1)+
#   theme_minimal()

plt_muts_count <- stats %>% 
  mutate(coverage=paste0(coverage,"x")) %>% 
  ggplot(aes(x=factor(coverage, levels = c("50x", "100x", "150x", "200x")), y=log10(mutation_count), fill=as.factor(purity))) + 
  geom_boxplot(outliers = FALSE)+
  ylab("Log10 mutation count")+
  xlab("coverage")+
  labs(fill="purity")+
  #scale_fill_manual(values = purity_colors)+
  scale_fill_manual(values = c('0.3' = '#A8CCDA', '0.6' = '#759DB9', '0.9' ='#426D98' ))+
  facet_wrap(~ type,scales = "free",ncol = 1,strip.position = "right")+
  my_ggplot_theme()

# stats %>% 
#   group_by(type,coverage,purity) %>%
#   summarise(median_count_muts = median(mutation_count)) %>%
#   ggplot(aes(x=coverage,color=as.factor(purity),group = purity))+
#   geom_point(aes(y = median_count_muts)) +
#   geom_line(aes(y = median_count_muts))+
#   facet_wrap(~type,scales="free")

# 
# all_stats_mean <- rbind(real_stats_mean,stats_mean)
# plt <- all_stats_mean %>% 
#   mutate(log_mut_counts=log10(sum_count)) %>% 
#   ggplot(aes(x=as.factor(coverage),y=sum_count,color = purity)) +
#   geom_point()+
#   geom_line(aes(group = interaction(purity, type)), size = 1)+
#   scale_color_manual(values=c(purity_colors,"real"="black"))+
#   scale_y_continuous(labels = scales::scientific,
#                      breaks = scales::breaks_pretty(n = 4))+
#   facet_wrap2(ncol = 1,vars(type),scales = "free_y")+
#   ylab("Mutation Counts") +
#   xlab("coverage") +
#   theme_classic()
# 
# plt_indel <- all_stats_mean %>% 
#   filter(type=="INDEL") %>% 
#   mutate(log_mut_counts=log10(sum_count)) %>% 
#   ggplot(aes(x=as.factor(coverage),y=sum_count,color = purity)) +
#   geom_point()+
#   geom_line(aes(group = interaction(purity)), size = 1)+
#   scale_y_break(c(5e5,1.5e6),space = 2,ticklabels=c(1e6,2e6)) +
#   scale_y_break(c(1.5e6,2e6),space = 2,ticklabels=c(2e6))+
#   # scale_y_continuous(
#   #   labels = scales::scientific,
#   #   breaks = c(0, 2e5, 4e5, 6e5, 1e6, 2e6, 3e6)  # specific y-axis ticks for INDEL
#   # ) +
#   scale_color_manual(values=c(purity_colors,"real"="black"))+
#   ylab("INDEL Counts") +
#   xlab("coverage") +
#   theme_classic()+
#   theme(legend.position = "bottom")
# 
# plt_snv <- all_stats_mean %>% 
#   filter(type=="SNV") %>% 
#   mutate(log_mut_counts=log10(sum_count)) %>% 
#   ggplot(aes(x=as.factor(coverage),y=sum_count,color = purity)) +
#   # geom_boxplot(outlier.alpha = 0.3)+
#   geom_point()+
#   # scale_y_continuous(labels = scales::scientific)+
#   geom_line(aes(group = interaction(purity)), size = 1)+
#   scale_y_break(c(5e6,1.5e7),space = 2,ticklabels=c(1e7,2e7)) + 
#   scale_y_break(c(1.5e7,2e7),space = 2,ticklabels=c(2e7,3e7))+
#   scale_color_manual(values=c(purity_colors,"real"="black"))+
#   ylab("SNV Counts") +
#   xlab("coverage") +
#   theme_classic()+
#   theme(legend.position = "none")
# 
# plt_all <- (plt_snv/plt_indel)
# # real_stats_plot <- real_stats %>% 
# #   mutate(log_mut_counts=log10(n)) %>% 
# #   ggplot(aes(x=spn,y=log_mut_counts,group=type)) +
# #   geom_bar(position="dodge", stat="identity",color="black", fill="white")+
# #   # scale_fill_manual(values=c("indel"="grey","SNV"="black"))+
# #   ylab("Log10 Mutation Counts") +
# #   xlab("SPN") +
# #   coord_flip()+
# #   theme_bw()
# 
# real_stats_plot <- real_stats %>% 
#   mutate(log_mut_counts=log10(n)) %>%
#   ggplot(aes(x=spn,y=n)) +
#   geom_point()+
#   geom_line(aes(group = interaction(type)), size = 0.5,alpha = 0.2)+
#   facet_wrap2(ncol = 1,vars(type),scales = "free_y")+
#   ylab("") +
#   xlab("spn") +
#   # theme(axis. = element_text(angle = 45, vjust = 0.5, hjust=1))+
#   theme_bw()
# 
# 
# 
params_grid = expand.grid(spns, coverages,purities)
colnames(params_grid) = c("spn","coverage", "purity")

muts_counts_stats <-list()
fastq_sizes <- list()
basedir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"


params_grid = expand.grid(spns,lots_groups)
colnames(params_grid) = c("spn","lot")

for (i in 1:nrow(params_grid)){
  spn <- params_grid[i,]$spn
  pur <- 0.9
  lot <- params_grid[i,]$lot
  dir_path <- file.path(basedir,spn,paste0("sequencing/tumour/purity_",pur),"/FASTQ/")
  files <- list.files(dir_path, full.names = TRUE, recursive = FALSE,pattern = "R1|R2")

  # Get file information
  file_info <- file.info(files)

  # Create a data frame with file name and size (in KB or MB)
  tmp <- data.frame(
    File = basename(rownames(file_info)),
    Size_GB = round(file_info$size / (1024^3), 3)) %>%
    separate(col = File,into = c("lot","spn","time"),sep = "_") %>%
    separate(col = time,into = c("time","space","R","ext1","ext2"),sep = "\\.") %>%
    mutate(sample=paste0(spn,"_",time,".",space)) %>%
    select(lot,spn,sample,Size_GB)
  nsamples <- length(unique(tmp$sample))
  # fastq_sizes[[i]] <- tmp
  fastq_sizes[[i]] <- tmp %>%
    head(lot*nsamples*2) %>%
    group_by(spn) %>%
    summarize(Size_GB_total=sum(Size_GB)) %>%
    mutate(coverage=case_when(lot==10 ~ 50,
                              lot==20 ~ 100,
                              lot==30 ~ 150,
                              lot==40 ~ 200)) %>%
    mutate(spn=spn)
}
stats_space <- do.call("rbind",fastq_sizes)
stats_space <- stats_space %>% 
  mutate(Size_TB_total_all_purities=(Size_GB_total*3)/(1024))

# stats_space_SPN05 <- stats_space %>% 
#   filter(spn=="SPN01") %>% 
#   mutate(spn="SPN05")
# 
# stats_space <- rbind(stats_space,stats_space_SPN05)

stats_space_plot1 <- stats_space %>%
  mutate(n_samples=case_when(spn%in%c("SPN01","SPN05")~"3",
                    spn%in%c("SPN02","SPN04")~"2",
                    spn%in%c("SPN03")~"4",
                    spn%in%c("SPN06","SPN07")~"5")) %>% 
  filter(coverage==200) %>%
  ggplot(aes(x=spn,y=Size_TB_total_all_purities,fill=n_samples)) +
  # ggplot(aes(x=coverage,y=Size_GB_total,fill=spn)) +
  # geom_boxplot()+
  geom_bar(position="dodge", stat="identity")+
  scale_fill_manual('N samples', values  = c('2' = '#c9e4ca', '3' ='#87bba2', '4' ='#55828b', '5' ='#3b6064', '1 - Normal' = 'gray'))  +
  theme_minimal()+
  ylab("Fastq Size (TB)")+
  xlab("SPN")+
  labs(caption = "Paired-end reads for maximum sequenced coverage of 200x")+
  my_ggplot_theme()

# 
wrap_plots(list(plt_muts_count,stats_space_plot1),guides="collect",design = "AABB\nAABB") + plot_annotation(tag_levels = 'A') & theme(legend.position = "bottom")
ggsave(filename = "/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_2/edf_space.pdf",device = "pdf",width = 8,height = 4,dpi = 300)
