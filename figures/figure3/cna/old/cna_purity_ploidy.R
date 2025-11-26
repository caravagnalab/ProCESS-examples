rm(list=ls())
library(ProCESS)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
# source("compute_FGA.R")
scout_dir <-"/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN06","SPN07")
COVERAGES <- c("50","100","150")
PURITIES <- c("0.3","0.6","0.9")
SPN_colors <-c("SPN01"='steelblue', "SPN02"='seagreen', "SPN03"='goldenrod', 
               "SPN04"='coral', "SPN06"='palevioletred', "SPN07"='indianred3')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
WGD_samples <- c("SPN01_1.1","SPN01_1.3","SPN06_3.1","SPN06_3.2")

validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION"
params_grid = expand.grid(COVERAGES, PURITIES)
colnames(params_grid) = c("coverage", "purity")
df_all_SPN <- list()
for (SPN in SPNS){
  spn <- SPN
  df = lapply(1:nrow(params_grid), function(i) {
    coverage = params_grid[i,]$coverage
    purity = params_grid[i,]$purity
    comb <- paste0(coverage,"x_",purity)
    samples <- get_sample_names(spn = spn)
    validation_dir_cna <- file.path(validation_dir,spn,"cna")
    all_metrics_comb <- list()
    for (sample in samples){
      metrics_filename <- paste0(validation_dir_cna,"/",comb,"/",sample,"/metrics.rds")
      metrics_bp_filename <- paste0(validation_dir_cna,"/",comb,"/",sample,"/metrics_bp.rds")
      if (!file.exists(metrics_filename)){
        message("File not found: ", metrics_filename)
        # all_metrics_comb[[sample]] <- NA
      } else{
        metrics_df <- readRDS(metrics_filename) %>% 
          mutate(delta_purity=as.numeric(true_purity)-as.numeric(purity)) %>% 
          mutate(delta_ploidy=as.numeric(true_ploidy)-as.numeric(ploidy)) %>% 
          filter(tool!="cnvkit")
        metrics_bp_df <- readRDS(metrics_bp_filename) %>% 
          filter(chr=="genome") %>% 
          filter(tool!="cnvkit")
        all_metrics_comb[[sample]] <- inner_join(metrics_df,metrics_bp_df,by=c("tool","sample","spn","coverage","fga","fgs","true_purity"))
      }
    }
    all_combinations_SPN <- do.call("rbind",all_metrics_comb)
  })
  df_all_SPN[[spn]] <- do.call("rbind",df)  
}
df_all_combs_SPN <- do.call("rbind",df_all_SPN)


df_all_combs_SPN <- df_all_combs_SPN %>% 
  mutate(fga_class=case_when(fga>=30 ~ "High FGA",
                             TRUE ~ "Low FGA")) %>% 
  mutate(tool = case_when(
    tool == 'ascat' ~ 'ASCAT',
    tool == 'sequenza' ~ 'Sequenza',
    tool == 'battenberg' ~ 'Battenberg',
  ))


#### PLOIDY ########
grouped_data <- df_all_combs_SPN  %>% 
  group_by(sample,tool) %>%
  mutate(mean_delta_ploidy=mean(delta_ploidy)) %>% 
  select(sample,fga_class,mean_delta_ploidy,tool) %>% 
  unique()


## Proposal 1
p_ploidy_per_fga <- ggplot(df_all_combs_SPN, aes(x = tool, y = delta_ploidy, fill = fga_class, shape=)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = list("High FGA"="indianred3","Low FGA"="lightsalmon1"))+
  labs(
    y = "Δ ploidy", x = "Tool"
  ) +
  theme_minimal(base_size = 13)



ggplot(df_all_combs_SPN, aes(x = true_purity, y = purity, color=tool,shape=fga_class, fill=tool)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21) +
  geom_jitter()+
  scale_fill_manual(values=col_cna_tools)+
  scale_color_manual(values=col_cna_tools)+
  scale_shape_manual(values = c("High FGA"=24,"Low FGA"=21))+
  labs(
    y = "estimated purity", x = "simulated purity"
  )+
  theme_minimal() #+
  # facet_wrap(~fga_class)+ 
  # geom_smooth(method = 'lm')

df_all_combs_SPN %>% 
  ggplot(aes(x = true_purity, y = delta_purity)) +
  geom_boxplot(aes(fill=tool,color=tool),
    alpha = 0.7) +
  geom_jitter(aes(shape=fga_class, fill=tool, x = true_purity, y = delta_purity,color=tool), inherit.aes = F)+
  scale_fill_manual(values=col_cna_tools)+
  scale_color_manual(values=col_cna_tools)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  # scale_size_manual(values = c("50"=3,"Low FGA"=1))+
  scale_shape_manual(values = c("High FGA"=16,"Low FGA"=1))+
  labs(
    y = "simulated purity - estimated purity", x = "simulated purity"
  )+
  theme_minimal()

########### PURITY #############
df_all_combs_SPN %>%
  mutate(correct_purity_estimation_class = case_when(
    abs(delta_purity) < 0.2 ~ "correctly estimated",
    delta_purity > 0.2 ~ "underestimated purity",
    delta_purity < -0.2 ~ "overestimated purity"
  )) %>%
  mutate(true_purity = paste0('Simulated Purity =', true_purity)) %>% 
  ggplot(aes(
    x = reorder(tool, delta_purity, \(x) sd(x, na.rm = TRUE)),  # reorders by SD
    y = delta_purity,
    group = interaction(tool, true_purity)
  )) +
  geom_boxplot(
    aes(color = tool, fill = tool),
    alpha = 0.2,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(
      color = fga_class,
      fill = fga_class,
      shape = correct_purity_estimation_class
    ),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
    size = 2,
    alpha = 0.8
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey40") +
  geom_hline(yintercept = 0.2, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -0.2, linetype = "dashed", color = "grey") +
  scale_color_manual(
    name = "FGA class",
    values = c("High FGA" = "indianred3", "Low FGA" = "grey40")
  ) +
  scale_fill_manual(
    name = "Tool",
    values = col_cna_tools
  ) +
  scale_shape_manual(
    name = "Estimation class",
    values = c(
      "correctly estimated" = 16,
      "underestimated purity" = 25,
      "overestimated purity" = 24
    )
  ) +
  facet_wrap(~true_purity) +
  labs(
    y = "simulated purity - estimated purity",
    x = "simulated purity"
  ) +
  theme_minimal()+
  theme(axis.text.x = element_blank())+
  +labs(x='')+
  



#####################################
########### PLOIDY #############
df_all_combs_SPN %>%
  mutate(correct_ploidy_estimation_class = case_when(
    abs(delta_ploidy) < 1 ~ "correctly estimated",
    delta_ploidy > 1 ~ "underestimated ploidy",
    delta_ploidy < -1 ~ "overestimated ploidy"
  )) %>%
  ggplot(aes(
    x = reorder(tool, delta_ploidy, \(x) sd(x, na.rm = TRUE)),  # reorders by SD
    y = delta_ploidy,
    group = interaction(tool, true_purity)
  )) +
  geom_boxplot(
    aes(color = tool, fill = tool),
    alpha = 0.2,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(
      color = fga_class,
      fill =  fga_class,
      shape = correct_ploidy_estimation_class
    ),
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
    size = 2,
    alpha = 0.8
  ) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey40") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -1, linetype = "dashed", color = "grey") +
  scale_color_manual(
    name = "FGA class",
    values = c("High FGA" = "indianred3", "Low FGA" = "grey40")
  ) +
  scale_fill_manual(
    name = "Tool",
    values = col_cna_tools
  ) +
  scale_shape_manual(
    name = "Estimation class",
    values = c(
      "correctly estimated" = 16,
      "underestimated ploidy" = 25,
      "overestimated ploidy" = 24
    )
  ) +
  facet_wrap(~true_purity) +
  labs(
    y = "simulated purity - estimated purity",
    x = "simulated purity"
  ) +
  theme_minimal()+
  theme(axis.text.x = element_blank())



############# PLOIDY #############


df_all_combs_SPN %>% 
  mutate(correct_purity_estimation_class=case_when(abs(delta_purity)<0.2~"correct",
                                                   TRUE~"uncorrect")) %>% 
  ggplot(aes(tool, delta_purity , 
             # color = correct_purity_estimation_class,
                             # fill=correct_purity_estimation_class,
                             shape=fga_class,alpha=correct_purity_estimation_class,
                             group = interaction(tool, true_purity))) +
  # geom_boxplot(alpha = 1,outliers = F) +
  geom_jitter(position = position_jitterdodge()) +
  # scale_fill_manual(values=c("correct"="grey","uncorrect"="black"))+
  # scale_color_manual(values=c("correct"="grey","uncorrect"="black"))+
  scale_alpha_manual(values=c("correct"=0.2,"uncorrect"=1))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_shape_manual(values = c("High FGA"=16,"Low FGA"=1))+
  # scale_size_manual(values = c("High FGA"=3,"Low FGA"=1))+
  labs(
    y = "simulated purity - estimated purity", x = "simulated purity"
  )+
  theme_minimal()+
  facet_wrap(~true_purity)

df_all_combs_SPN %>% 
  ggplot(aes(x = true_purity, y = delta_purity)) +
  geom_boxplot(aes(fill = tool),
               alpha = 0.7, outlier.shape = NA) +   # hide default outliers
  geom_jitter(aes(shape = fga_class, color = tool,group=tool),
              width = 0.15, size = 2, alpha = 0.8) + # add your own outlier/point layer
  scale_fill_manual(values = col_cna_tools) +
  scale_color_manual(values = col_cna_tools) +
  scale_shape_manual(values = c("High FGA" = 16, "Low FGA" = 1)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  labs(
    y = "simulated purity - estimated purity",
    x = "simulated purity"
  ) +
  theme_minimal()

## Proposal 2
p_ploidy_per_tool <- ggplot(df_all_combs_SPN, aes(x = as.factor(true_purity), y = delta_ploidy, fill = tool,shape=fga_class)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21,outlier.alpha = 0.3) +
  geom_jitter(alpha = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_shape_manual(values = c("High FGA"=24,"Low FGA"=21))+
  scale_fill_manual(values=col_cna_tools)+
  labs(
    y = "Δ ploidy", x = "simulated purity"
  ) +
  # facet_wrap(~fga_class)+
  theme_minimal()


library(dplyr)
library(ggplot2)
library(ggrepel)


df_metrics <- df_all_combs_SPN %>%
  group_by(tool, true_purity) %>%
  summarise(
    rmse = sqrt(mean((ploidy - true_ploidy)^2, na.rm = TRUE)),
    r = cor(ploidy, true_ploidy, use = "complete.obs", method = "pearson"),
    x_pos = max(true_ploidy, na.rm = TRUE) * 0.9,
    y_pos = max(ploidy, na.rm = TRUE) * 0.9
  ) %>%
  mutate(
    label = paste0(
      "RMSE = ", round(rmse, 3), "\n",
      "R = ", round(r, 3)
    )
  ) %>%
  ungroup()

ggplot(df_all_combs_SPN, aes(x = true_ploidy, y = ploidy)) +
  geom_point(aes(color=fga_class)) +
  geom_smooth(method = "lm",alpha = 0.2,aes(fill = as.factor(true_purity),colour = as.factor(true_purity))) +
  scale_color_manual(values = purity_colors) +
  facet_wrap(~tool) +
  geom_label_repel(
    data = df_metrics,
    aes(x = x_pos, y = y_pos, label = label, fill = as.factor(true_purity)),
    inherit.aes = FALSE,
    color = "black",
    size = 3.5,
    show.legend = FALSE,
    label.padding = unit(0.04, "lines"),
    label.r = unit(0.15, "lines"),
    label.size = 0.25
  ) +
  scale_fill_manual(values = purity_colors) +
  scale_color_manual(values = list("High FGA"="indianred3","Low FGA"="grey40"))+
  theme_minimal() +
  labs(color = "True purity", fill = "True purity",x="Simulated ploidy",y="Estimated ploidy")



#### PURITY ########
grouped_data <- df_all_combs_SPN  %>% 
  group_by(sample,tool,true_purity) %>%
  mutate(mean_delta_purity=mean(delta_purity)) %>% 
  select(sample,true_purity,mean_delta_purity,tool) %>% 
  unique()



## Proposal 2
p_purity <- ggplot(df_all_combs_SPN, aes(x = true_purity, y = delta_purity, fill = tool)) +
  geom_boxplot(alpha = 0.7, outlier.shape = 21,outlier.alpha = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  # scale_fill_manual(values=purity_colors)+
  scale_fill_manual(values=col_cna_tools)+
  # scale_fill_manual(values = list("High FGA"="indianred3","Low FGA"="lightsalmon1"))+
  labs(
    y = "Δ purity", x = "simulated purity"
  ) +
  theme_minimal()+
  facet_wrap(~fga_class)

p_ploidy_per_tool/p_purity
ggsave(plot = p_ploidy_per_tool, filename = 'p_ploidy_per_tool.pdf',width = 8, height = 3, units = 'in')
ggsave(plot = p_purity, filename = 'p_purity.pdf',width = 8, height = 3, units = 'in')

wrap_plots(list(p_ploidy_per_tool,p_purity),ncol=1,guides = "collect") & theme(legend.position = "bottom")
cna_panel <- p_ploidy_per_tool + p_purity + plot_layout(nrow = 2,guides = "collect") + plot_annotation(tag_levels = 'A') & theme(legend.position = "bottom")

### PROPOSAL CORRELATION PLOIDY - PURITY
ggplot(df_all_combs_SPN, aes(x = delta_ploidy, y = delta_purity, col = spn,alpha=fga))+
  geom_point(size = 2)+
  # geom_density_2d()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = SPN_colors)+
  facet_wrap(~tool)+
  theme_bw()

ggplot(df_all_combs_SPN, aes(x = delta_ploidy, y = delta_purity, col = spn, shape = tool)) +
  geom_point(size = 2) +
  # geom_density_2d(color = "black", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = SPN_colors) +
  facet_wrap(~fga_class) +
  theme_bw()


ggplot(df_all_combs_SPN, aes(x = delta_ploidy, y = delta_purity, col = spn, shape = tool)) +
  geom_point(size = 2) +
  # geom_density_2d(color = "black", alpha = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  scale_color_manual(values = SPN_colors) +
  facet_wrap(~fga_class) +
  theme_bw()

## color by fga class
df_all_combs_SPN %>% 
  group_by(spn,tool,true_purity) %>% 
  mutate(mean_delta_purity=mean(delta_purity)) %>% 
  mutate(mean_delta_ploidy=mean(delta_ploidy)) %>% 
  ggplot(aes(x = mean_delta_ploidy, y = mean_delta_purity, col = fga_class,shape=tool))+
  geom_point(size = 3)+
  # geom_density_2d()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  # geom_rect(xmin=-0.5, xmax=0.5, ymin=-0.2,ymax=0.2,alpha=.05)+
  # scale_color_manual(values = SPN_colors)+
  scale_color_manual(values = list("High FGA"="indianred3","Low FGA"="lightsalmon1"))+
  facet_wrap(~true_purity)+
  theme_bw()
## color by spn
df_all_combs_SPN %>% 
  group_by(spn,tool,true_purity) %>% 
  mutate(mean_delta_purity=mean(delta_purity)) %>% 
  mutate(mean_delta_ploidy=mean(delta_ploidy)) %>% 
  ggplot(aes(x = mean_delta_ploidy, y = mean_delta_purity, col = spn,alpha=fga,shape=tool))+
  geom_point(size = 3)+
  # geom_density_2d()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  # geom_rect(xmin=-0.5, xmax=0.5, ymin=-0.2,ymax=0.2,alpha=.05)+
  scale_color_manual(values = SPN_colors)+
  facet_wrap(~true_purity)+
  theme_bw()

df_all_combs_SPN %>% 
  group_by(spn,tool,true_purity) %>% 
  # mutate(mean_delta_purity=mean(delta_purity)) %>% 
  # mutate(mean_delta_ploidy=mean(delta_ploidy)) %>% 
  ggplot(aes(x = delta_ploidy, y = delta_purity, col = spn,alpha=fga,shape=true_purity))+
  geom_point(size = 3)+
  # geom_density_2d()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  # geom_rect(xmin=-0.5, xmax=0.5, ymin=-0.2,ymax=0.2,alpha=.05)+
  scale_color_manual(values = SPN_colors)+
  facet_wrap(~tool)+
  theme_bw()


df_all_combs_SPN %>% 
  group_by(spn,tool,true_purity) %>% 
  # mutate(mean_delta_purity=mean(delta_purity)) %>% 
  # mutate(mean_delta_ploidy=mean(delta_ploidy)) %>% 
  ggplot(aes(x = delta_ploidy, y = delta_purity, col = fga_class,alpha=true_purity))+
  geom_point(size = 3)+
  # geom_density_2d()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  # geom_rect(xmin=-0.5, xmax=0.5, ymin=-0.2,ymax=0.2,alpha=.05)+
  facet_wrap(~tool)+
  scale_color_manual(values = list("High FGA"="indianred3","Low FGA"="lightsalmon1"))+
  theme_bw()


corr_plot <- df_all_combs_SPN %>% 
  # group_by(sample,tool,true_purity) %>% 
  # mutate(mean_delta_purity=mean(delta_purity)) %>%
  # mutate(mean_delta_ploidy=mean(delta_ploidy)) %>%
  # ggplot(aes(x = mean_delta_ploidy, y = mean_delta_purity, col = true_purity,shape=tool))+
  ggplot(aes(x = delta_ploidy, y = delta_purity, col = true_purity,shape=tool))+
  geom_point(size = 2)+
  # geom_density_2d()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  # geom_rect(xmin=-0.5, xmax=0.5, ymin=-0.2,ymax=0.2,alpha=.05)+
  facet_wrap(~fga_class,nrow=1)+
  # facet_grid(tool~fga_class)+
  scale_color_manual(values = purity_colors)+
  theme_minimal()

df_all_combs_SPN %>% 
  group_by(sample,tool,true_purity) %>% 
  mutate(mean_delta_purity=mean(delta_purity)) %>%
  mutate(mean_delta_ploidy=mean(delta_ploidy)) %>%
  ggplot(aes(x = mean_delta_purity, y = mean_delta_ploidy, col = fga_class,alpha=true_purity,shape=tool))+
  geom_point(size = 3)+
  # geom_density_2d()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  # geom_rect(xmin=-0.5, xmax=0.5, ymin=-0.2,ymax=0.2,alpha=.05)+
  scale_color_manual(values = list("High FGA"="indianred3","Low FGA"="lightsalmon1"))+
  scale_alpha_manual(values = c("0.3"=0.5,"0.6"=0.8,"0.9"=1))+
  theme_bw()


p=df_all_combs_SPN %>% 
  # group_by(spn,tool,true_purity) %>%
  # group_by(sample,tool,true_purity) %>% 
  # mutate(mean_delta_purity=mean(delta_purity)) %>%
  # mutate(mean_delta_ploidy=mean(delta_ploidy)) %>%
  ggplot(aes(x = delta_purity, y = delta_ploidy, col = true_purity,alpha=fga))+
  geom_point(size = 3)+
  # geom_density_2d()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  # geom_rect(xmin=-0.5, xmax=0.5, ymin=-0.2,ymax=0.2,alpha=.05)+
  # facet_wrap(~tool)+
  scale_color_manual(values = purity_colors)+
  theme_bw()

p2 <- ggExtra::ggMarginal(p, type="density")

p2
