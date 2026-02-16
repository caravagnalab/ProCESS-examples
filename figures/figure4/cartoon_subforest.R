library(tidyverse)
library(tidyr)
library(ggplot2)
library(ProCESS)
library(patchwork)
library(ggExtra)

my_theme = theme_light(base_size=8) +
  theme(legend.position="right",
        legend.key.height=unit(0.1, "cm"),
        legend.key.width=unit(0.1, "cm"),
        legend.key.spacing=unit(0.1, "cm"),
        panel.background=element_rect(fill="white"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title=element_text(size=8),
        legend.text=element_text(size=6, margin=margin(l=0.1, unit="cm")),
        legend.title=element_text(size=6),
        text=element_text(size=8,family="Helvetica"),
        panel.grid = element_blank())

base = '/orfeo/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/demo_simulation/'
setwd(base)
save = paste0(base, 'cartoon/')

source('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/utils_plots.R')

# 0. Create a single table process + viber ####
final_table = readRDS(paste0(save, "tables/table_viber_interpreted.rds"))

table_wide_complete = final_table %>%
  pivot_wider(values_from="VAF", names_from="sample_id")

# N = 200  # number of Subclonal rows to keep
N = 300  # number of Subclonal rows to keep

set.seed(122)

M = 1500

table_wide <- table_wide_complete %>%
  filter(is_driver_process == TRUE) %>%
  bind_rows(
    table_wide_complete %>%
      filter(is_driver_process == FALSE, cluster_id_process != "Subclonal") %>%
      slice_sample(n = M)
  ) %>%
  bind_rows(
    table_wide_complete %>%
      filter(is_driver_process == FALSE, cluster_id_process == "Subclonal") %>%
      slice_sample(n = N)
  )

sample_forest = load_sample_forest(paste0("results/sample_forest.sff"))

sample_names <- final_table$sample_id %>% unique()
sample_names = stringr::str_sort(sample_names, numeric = TRUE)

sample_names = as.character(sample_names)
cm = combn(sample_names, 2)

# 1. Create process plots (scatter + tree) ####
### Scatter ####
# As the plot functions work with is_driver, gene and cluster_id, rename the columns

table_wide_process = table_wide %>%
  dplyr::rename(gene=driver_true, is_driver=is_driver_true, cluster_id=cluster_id_process)

table_wide_process = table_wide_process %>%
  filter(!is.na(cluster_id))

color_palette_tool = RColorBrewer::brewer.pal(n = length(unique(final_table$cluster_id_tool)), name = "Dark2") %>% 
  setNames(str_sort(unique(final_table$cluster_id_tool), numeric=T))
color_palette_tool['Other'] = "gainsboro" #"#808080"

color_palette_process = RColorBrewer::brewer.pal(n = length(unique(final_table$cluster_id_tool))+1, name = "Dark2")

color_palette_process['Clone 2'] = color_palette_tool['C2']
color_palette_process['Clone 1'] = color_palette_tool['C4']
color_palette_process['Other'] = "gainsboro" #"#808080"

point_size = .1
point_alpha = .5

scatter_process <- apply(
  cm,
  2,
  function(w) plot_scatter_driver_cartoon(table_wide_process, s1 = w[1], s2 = w[2], color_palette_process, psize = point_size, palpha = point_alpha) # w is the sample name
  )[[1]]
scatter_process <- scatter_process + xlab('VAF (Time-point 1)') +  ylab('VAF (Time-point 2)')



# 2. Create tool plots (scatter + tree)  and a fake/wrong driver ####

### Scatter ####

table_wide_tool_plot = table_wide %>%
  dplyr::rename(gene=driver_true, is_driver=is_driver_true, cluster_id = cluster_id_tool) %>%
  filter(!is.na(cluster_id))

scatter_tool <- apply(
  cm,
  2,
  function(w) plot_scatter_driver_cartoon(table_wide_tool_plot, s1 = w[1], s2 = w[2], color_palette_tool, psize = point_size, palpha = point_alpha) # w is the sample name
  
)[[1]]

scatter_tool <- scatter_tool + xlab('VAF (Time-point 1)') +  ylab('VAF (Time-point 2)')

## Scatter interpreted ####

table_wide_tool_interpreted = table_wide %>%
  dplyr::rename(gene=driver_true, is_driver=is_driver_true, cluster_id = cluster_id_tool_interpreted) %>%
  filter(!is.na(cluster_id))

scatter_tool_interpreted <- apply(
  cm,
  2,
  function(w) plot_scatter_driver_cartoon(table_wide_tool_interpreted, s1 = w[1], s2 = w[2], color_palette_tool, psize = point_size, palpha = point_alpha) # w is the sample name
  )[[1]]
scatter_tool_interpreted <- scatter_tool_interpreted + xlab('VAF (Time-point 1)') +  ylab('VAF (Time-point 2)')

# Marginal univariate ####

### Process ####
marginals_table = final_table
samples = c("Sample 1","Sample 2")

process_marginals = lapply(samples, function(sample_name){
  table_sample = marginals_table %>% filter(sample_id == sample_name)

# Count the number of data points per cluster
cluster_order = table_sample %>%
  dplyr::count(cluster_id_process) %>%       
  arrange(desc(n)) %>%
  pull(cluster_id_process)

table_sample = table_sample %>%
  mutate(cluster_id_process = factor(cluster_id_process, levels = cluster_order))

# Here plot the marginal of process

if(sample_name == 'Sample 1'){
  table_sample %>% 
    filter(VAF>0.05) %>% 
    ggplot() +
    geom_histogram(aes(x=VAF, fill=cluster_id_process), position="identity", alpha=1, bins=100) +
    ylab("") +
    # ylim(0,2000)+
    xlim(0,1)+
    ggtitle(sample_name)+
    theme_minimal()+
    scale_fill_manual(values=color_palette_process, name="Cluster") +
    theme(legend.position="bottom")
}else{
  table_sample %>% 
    filter(VAF>0.05) %>% 
    ggplot() +
    geom_histogram(aes(y=VAF, fill=cluster_id_process), position="identity", alpha=1, bins=100) +
    ylab("") +
    # ylim(0,2000)+
    ylim(0,1)+
    ggtitle(sample_name)+
    theme_minimal()+
    scale_fill_manual(values=color_palette_process, name="Cluster") +
    theme(legend.position="bottom")
}


})

process_marginals

# Put marginals on scatter
p_top = process_marginals[[1]]+
  coord_cartesian(xlim = c(0,1), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )+
  ggtitle('')

p_right = process_marginals[[2]] +
  # coord_cartesian(xlim = c(0,1), expand = FALSE) +
  # coord_flip() +
  coord_cartesian(ylim = c(0,1), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )+
  ggtitle('')

zero_margins = theme(plot.margin = margin(0, 0, 0, 0))

p_scatter = scatter_process +
  theme(legend.position="bottom") + 
  zero_margins

p_top = p_top + zero_margins
p_right = p_right + zero_margins

scatter_process_w_marginals = 
  p_top + plot_spacer() + 
  p_scatter + p_right + plot_spacer() + 
  plot_layout(design = 'ABB\nCDE',
              heights = c(1, 4),
              widths = c(4, 1)
              ) 


scatter_process_w_marginals = wrap_elements(full = scatter_process_w_marginals)
scatter_process_w_marginals 

ggsave('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/cartoon/marginals_process.pdf', 
       scatter_process_w_marginals, width = 3, height=3, unit = 'in')

### Blind ####
marginals_tool_table = final_table
samples = c("Sample 1","Sample 2")

blind_marginals = lapply(samples, function(sample_name){
  table_sample = marginals_tool_table %>% filter(sample_id == sample_name)

cluster_order = table_sample %>%
  dplyr::count(cluster_id_tool) %>% 
  arrange(desc(n)) %>% 
  pull(cluster_id_tool)

table_sample = table_sample %>%
  mutate(cluster_id_tool = factor(cluster_id_tool, levels = cluster_order))

# Here plot the marginal of process
table_sample %>% 
  filter(VAF>0.05) %>% 
  ggplot() +
  geom_histogram(aes(x=VAF, fill=cluster_id_tool), position="identity", alpha=1, bins=100) +
  ylab("") +
  # ylim(0,2000)+
  xlim(0,1)+
  ggtitle(sample_name)+
  theme_minimal()+
  scale_fill_manual(values=color_palette_tool, name="Cluster") +
  theme(legend.position="bottom")

})

blind_marginals

# Put marginals on scatter
p_top = blind_marginals[[1]]+
  coord_cartesian(xlim = c(0,1), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )+
  ggtitle('')

p_right = blind_marginals[[2]] +
  coord_cartesian(xlim = c(0,1), expand = FALSE) +
  coord_flip() +
  # coord_cartesian(ylim = c(0,1), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )+
  ggtitle('')

zero_margins = theme(plot.margin = margin(0, 0, 0, 0))

p_scatter = scatter_tool +
  theme(legend.position="bottom") + 
  zero_margins

p_top = p_top + zero_margins
p_right = p_right + zero_margins

# scatter_blind_w_marginals = p_top + plot_spacer() + p_scatter + p_right +plot_spacer() + 
#   plot_layout(design = 'ABB\nCDE')
scatter_blind_w_marginals = 
  p_top +  plot_spacer() + 
  p_scatter + 
  p_right + plot_spacer() + 
  plot_layout(design = 'ABB\nCDE',
              heights = c(1, 4),
              widths = c(4, 1)
  )  &
  theme(panel.spacing = unit(0, "pt"))
scatter_blind_w_marginals

scatter_blind_w_marginals = wrap_elements(full = scatter_blind_w_marginals)
scatter_blind_w_marginals
ggsave('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/cartoon/marginals_tool_blind.pdf', 
       scatter_blind_w_marginals, width = 3, height=3, unit = 'in')

## Interpreted ####

interpreted_marginals = lapply(samples, function(sample_name){
  table_sample = marginals_tool_table %>% filter(sample_id == sample_name)
  
  # Sort clusters by size (descending order)
  cluster_order = table_sample %>%
    dplyr::count(cluster_id_tool_interpreted) %>%
    arrange(desc(n)) %>% 
    pull(cluster_id_tool_interpreted)
  
  table_sample = table_sample %>%
    mutate(cluster_id_tool_interpreted = factor(cluster_id_tool_interpreted, levels = cluster_order))
  
  # Here plot the marginal of process
  table_sample %>% 
    filter(VAF>0.05) %>% 
    ggplot() +
    geom_histogram(aes(x=VAF, fill=cluster_id_tool_interpreted), position="identity", alpha=1, bins=100) +
    ylab("") +
    # ylim(0,2000)+
    xlim(0,1)+
    ggtitle(sample_name)+
    theme_minimal()+
    scale_fill_manual(values=color_palette_tool, name="Cluster") +
    theme(legend.position="bottom")
  
})

interpreted_marginals

# Put marginals on scatter
p_top = interpreted_marginals[[1]]+
  coord_cartesian(xlim = c(0,1), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )+
  ggtitle('')

p_right = interpreted_marginals[[2]] +
  coord_cartesian(xlim = c(0,1), expand = FALSE) +
  coord_flip() +
  # coord_cartesian(ylim = c(0,1), expand = FALSE) +
  theme_void() +
  theme(
    legend.position = "none",
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0)
  )+
  ggtitle('')

zero_margins = theme(plot.margin = margin(0, 0, 0, 0))

p_scatter = scatter_tool_interpreted +
  theme(legend.position="bottom") + 
  zero_margins

p_top = p_top + zero_margins
p_right = p_right + zero_margins

# scatter_interpreted_w_marginals = p_top + plot_spacer() + p_scatter + p_right + plot_spacer() + 
#   plot_layout(design = 'ABB\nCDE')
scatter_interpreted_w_marginals = 
  p_top + plot_spacer() + 
  p_scatter + p_right + plot_spacer() + 
  plot_layout(design = 'ABB\nCDE',
              heights = c(1, 4),
              widths = c(4, 1)
  ) 

scatter_interpreted_w_marginals = wrap_elements(full = scatter_interpreted_w_marginals)
scatter_interpreted_w_marginals
ggsave('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/cartoon/marginals_tool_interpreted.pdf', 
       scatter_interpreted_w_marginals, width = 3, height=3, unit = 'in')

# Plot forest ####
source('/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/validation/Subclonal_deconvolution/demo_simulation/cartoon/utils_plots.R')

sample_forest = load_sample_forest(paste0("results/sample_forest_longitudinal.sff"))
sticks_process = get_relevant_branches(sample_forest)
sticks_process = sticks_process %>% mutate(label=ifelse(label=='Subclonal', 'Other', label)) %>% 
  mutate(label=ifelse(label=='Truncal', 'Clone 1', label))

plot_sticks_process = my_plot_sticks(forest=sample_forest, 
            labels=sticks_process, 
            cls = color_palette_process)

plot_sticks_process
ggsave(filename = paste0("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/cartoon/tree_process_longitudinal.pdf"), plot = plot_sticks_process,
       device="pdf", width=2, height=4, units="in")


plot_sticks_process + scatter_process_w_marginals + 
  scatter_blind_w_marginals +
  scatter_interpreted_w_marginals

# # OLD ####
# # Create final plot ####
# # Final plot subforest ####
# width=15
# height=15
# my_theme = theme_light(base_size=8) +
#   theme(legend.position="right",
#         legend.key.height=unit(0.1, "cm"),
#         legend.key.width=unit(0.1, "cm"),
#         legend.key.spacing=unit(0.1, "cm"),
#         panel.background=element_rect(fill="white"),
#         axis.text.x=element_text(size=6),
#         axis.text.y=element_text(size=6),
#         axis.title=element_text(size=8),
#         legend.text=element_text(size=6, margin=margin(l=0.1, unit="cm")),
#         legend.title=element_text(size=6),
#         text=element_text(size=8,family="Helvetica"),
#         panel.grid = element_blank())
# 
# patch_pr = patchwork::wrap_plots(
#   plot_sticks_process +
#     labs(title="Clonal evolution"),
#   process_marginals[[1]]+ 
#     guides(color = "none", fill = "none")+
#     labs(title="Marginal sample 1"),
#   process_marginals[[2]]+ 
#     guides(color = "none", fill = "none") +
#     labs(title="Marginal sample 2"),
#   scatter_process + 
#     # guides(color = "none", fill = "none")+ 
#     labs(title="Process clusters"),
#   scatter_tool + 
#     # guides(color = "none", fill = "none") +
#     labs(title="Blind assignment"),
#   scatter_tool_interpreted + 
#     # guides(color = "none", fill = "none") + 
#     labs(title="Interpreted assignment"),
#   design="aab\naac\nefd")&
#   my_theme &
#   theme(legend.position = 'bottom')
# patch_pr
# 
# ggsave(filename = paste0(save, "plots/patch_longitudinal.png"), plot = patch_pr,
#        device="png", width=width, height=height, units="cm", dpi=600)
# ggsave(filename = paste0(save, "plots/patch_longitudinal.pdf"), plot = patch_pr,
#        device="pdf", width=width, height=height, units="cm")
# 
# # Scatter with marginals ####
# # Final plot subforest ####
# width=15
# height=15
# my_theme = theme_light(base_size=8) +
#   theme(legend.position="right",
#         legend.key.height=unit(0.1, "cm"),
#         legend.key.width=unit(0.1, "cm"),
#         legend.key.spacing=unit(0.1, "cm"),
#         panel.background=element_rect(fill="white"),
#         axis.text.x=element_text(size=6),
#         axis.text.y=element_text(size=6),
#         axis.title=element_text(size=8),
#         legend.text=element_text(size=6, margin=margin(l=0.1, unit="cm")),
#         legend.title=element_text(size=6),
#         text=element_text(size=8,family="Helvetica"),
#         panel.grid = element_blank())
# 
# patch_pr = patchwork::wrap_plots(
#   plot_sticks_process +coord_flip()+
#     scale_y_reverse()+
#     labs(title="Clonal evolution"),
#   scatter_process_w_marginals + 
#     labs(title="Process clusters"),
#   scatter_blind_w_marginals +
#     labs(title="Blind assignment"),
#   scatter_interpreted_w_marginals +
#     labs(title="Interpreted assignment"),
#   # design="ab\nac\nad")&
#   design="aaa\nbcd")&
#   my_theme &
#   theme(legend.position = 'bottom')
# patch_pr
# 
# ggsave(filename = paste0(save, "plots/patch_longitudinal.png"), plot = patch_pr,
#        device="png", width=width, height=height, units="cm", dpi=600)
# ggsave(filename = paste0(save, "plots/patch_longitudinal.pdf"), plot = patch_pr,
#        device="pdf", width=width, height=height, units="cm")
# 
