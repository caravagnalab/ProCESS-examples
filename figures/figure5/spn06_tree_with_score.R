rm(list=ls())
library(ComplexHeatmap)
library(tidyr)
library(dplyr)
library(tibble)
library(tidygraph)
library(ggraph)
library(ggrepel)
library(grid)
library(ggforce)
setwd("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_mutect2_ascat")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/tumourevo_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure5/utils_ctree.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
# source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples//validation/Subclonal_deconvolution/utils_plots_final.R")
# source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples//validation/Subclonal_deconvolution//generate_table_main.R")
## source the function inf ctree
## ctree/R/utils_trees_transformations.R
#data <- readRDS("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_mutect2_ascat/res_int/50x_0.9p_SigProfiler_viber_df.rds")
SPN="SPN06"
cov=150
pur=0.9
vcf_caller="mutect2"
cna_caller="ascat"
sig_caller="SigProfiler"
dec_caller="viber"
data <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_",vcf_caller,"_",cna_caller,
 "/res_int/",cov,"x_",pur,"p_",sig_caller,"_",dec_caller,"_df.rds"))

sign_assignment_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/assign_signature/"
sign_type='SBS'
activity_df <- read.table(file.path(sign_assignment_dir,SPN,paste0(cov,"x_",pur,"p_",vcf_caller,"_",cna_caller),paste0(dec_caller,"_",sig_caller,"_raw"),"SBS96/Assignment_Solution/Activities/Assignment_Solution_Activities.txt"),header = T,sep="\t")
activity_df <- activity_df %>% pivot_longer(
  cols = starts_with(sign_type),
  names_to = "Signature",
  values_to = "Exposure") %>% 
  dplyr::rename(cluster=Samples)

if (dec_caller=="pyclonevi"){
  activity_df <- activity_df %>% 
    mutate(cluster=paste0("C",cluster))
}








data_sel <- data %>% 
  mutate(tier=case_when(
    is_clonal_tool                          ~ 'Tier 1',
    score_all >= .9                        ~ 'Tier 1',
    score_all >= .55 & score_all <.9           ~ 'Tier 2',
    score_all >= .2  & score_all <.55          ~ 'Tier 3',
    score_all <  .2                        ~ 'Tier 4'
  )) %>% 
  filter(spn==SPN) %>% 
  filter(coverage==cov,purity==pur) %>% 
  select(cluster,tier,score_all)

ctree <- get_tumourevo_subclonal(spn = SPN,coverage = cov,purity = pur,tool = "ctree",
                                 vcf_caller = vcf_caller,cna_caller=cna_caller,sample=SPN)
tier_colors <- c(
  "Tier 1" = "palegreen4",
  "Tier 2" = "goldenrod",
  "Tier 3" = alpha("indianred",0.6),
  "Tier 4" = alpha("grey",0.6),
  "Missing"  = alpha("gainsboro",0.6)
)

colors_cluster = c('indianred', 
                   'steelblue', 
                   'forestgreen', 
                   'goldenrod', 
                   'darkorange3', 
                   'palevioletred', 
                   'mediumpurple', 
                   'cornsilk4', 
                   'olivedrab3', 
                   'steelblue4', 
                   'indianred4',
                   'aquamarine3',
                   'saddlebrown',
                   'deeppink2',
                   'cornflowerblue',
                   'black')
names(colors_cluster) = paste0('C',0:15)
if (dec_caller=="viber"){
  ctree_obj<- readRDS(ctree$ctree_VIBER_rds)
} else{
  ctree_obj<- readRDS(ctree$ctree_pyclonevi_rds)  
}


tree = ctree_obj[[1]]
CCF_updated <- tree$CCF %>% 
  left_join(data_sel)
tree$CCF <- CCF_updated
tb_tree = tree$tb_adj_mat %>% 
  left_join(data_sel)

cex = 1

clones_orderings = igraph::topo_sort(igraph::graph_from_adjacency_matrix(DataFrameToMatrix(tree$transfer$clones)),
                                     mode = 'out')$name

nDrivers = length(clones_orderings) - 1 # avoid GL
node_palette = colorRampPalette(RColorBrewer::brewer.pal(n = 9, "Set1"))
drivers_colors = c('white', node_palette(nDrivers))
names(drivers_colors) = clones_orderings


# Add non-driver nodes, with the same colour
non_drivers = tb_tree %>%
  activate(nodes) %>%
  filter(!is.driver) %>%
  pull(cluster) # GL is not selected because is NA for is.driver

non_drivers_colors = rep("gainsboro", length(non_drivers))
names(non_drivers_colors) = non_drivers

tb_node_colors = c(drivers_colors, non_drivers_colors)

# Plot call
layout_ctree <- create_layout(tb_tree, layout = "tree")



pie_data <- layout_ctree %>%
  select(x, y, cluster,nMuts) %>%
  distinct() %>%
  left_join(activity_df, by = "cluster") %>%
  mutate(Exposure = replace_na(Exposure, 0)) %>%
  group_by(cluster) %>%
  mutate(
    frac = Exposure / sum(Exposure, na.rm = TRUE),
    ymax = cumsum(frac),
    ymin = lag(ymax, default = 0)
  ) %>%
  ungroup()
pie_data <- pie_data %>%
  mutate(
    r_node = 0.1 * cex + 0.25 * cex * (nMuts / max(nMuts, na.rm = TRUE))
  ) #%>% 
  # mutate(r_node=case_when(cluster=="GL"~max(r_node)))


ctree_plot <- ggraph(layout_ctree) +
  
  # ---------------------------
# edges
# ---------------------------
geom_edge_link(
  arrow = arrow(length = unit(2 * cex, 'mm')),
  end_cap = circle(5 * cex, 'mm'),
  start_cap = circle(5 * cex, 'mm')
) +
  
  # ---------------------------
# PIE CHART NODES
# ---------------------------
geom_arc_bar(
  data = pie_data,
  aes(
    x0 = x,
    y0 = y,
    r0 = 0,
    r = 0.25 * cex,
    start = 2 * pi * ymin,
    end   = 2 * pi * ymax,
    fill = Signature
  ),
  color = "white",
  linewidth = 0.2,
  inherit.aes = FALSE
) +
  
  # ---------------------------
# driver labels (cluster color)
# ---------------------------
geom_label_repel(
  aes(
    x = x,
    y = y,
    label = driver,
    colour = cluster
  ),
  na.rm = TRUE,
  nudge_x = .4,
  nudge_y = .4,
  size = 2.5 * cex
) +
  
  # ---------------------------
# cluster labels (cluster color)
# ---------------------------
geom_node_text(
  aes(
    label = cluster,
    colour = cluster
  ),
  vjust = 0,fontface = "bold",
  hjust = -3
) +
  
  # cluster colour scale
  scale_color_manual(values = colors_cluster) +
  
  # ---------------------------
# NEW colour scale for tier (circle border)
# ---------------------------
ggnewscale::new_scale_colour() +
  
  geom_circle(
    data = distinct(layout_ctree, x, y, cluster, tier),
    aes(
      x0 = x,
      y0 = y,
      r = 0.25 * cex,
      color = tier
    ),
    inherit.aes = FALSE,
    linewidth = 1,
    alpha = 0.6
  ) +
  
  scale_color_manual(values = tier_colors) +
  
  # ---------------------------
# layout / theme
# ---------------------------
coord_cartesian(clip = "off") +
  
  theme_void(base_size = 8 * cex) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(3 * cex, "mm")
  ) +
  
  # ---------------------------
# fill scale (pies)
# ---------------------------
scale_fill_manual(values = sbs_colors) +
  
  # ---------------------------
# guides
# ---------------------------
guides(
  fill = guide_legend("Signature"),
  colour = guide_legend("Tier")
) +
  
  # ---------------------------
# titles
# ---------------------------
labs(
  title = paste(tree$patient),
  subtitle = paste0("Coverage: ",cov,"\nPurity: ",pur,"\nSarek: ",
                    paste0(vcf_caller,", ",cna_caller),"\nTumourevo: ",
                    paste0(dec_caller,", ",sig_caller))
)
ggsave(filename = paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/spn06_phylo_ctree/spn06_ctree_",vcf_caller,"_",cna_caller,"_",sig_caller,"_",dec_caller,"_",cov,"x_",pur,"p.pdf"),
       width = 5,height = 7,plot = ctree_plot,dpi = 300)
