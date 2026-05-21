library(tidyr)
library(dplyr)
library(tibble)
library(ProCESS)
library(optparse)
setwd("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_mutect2_ascat")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/tumourevo_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/report/plotting/utils.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN07'),
                    make_option(c("--purity"), type = "double", default = 0.9),
                    make_option(c("--coverage"), type = "integer", default = 150),
                    make_option(c("--cna_caller"), type = "character", default = 'ascat'),
                    make_option(c("--vcf_caller"), type = "character", default = 'mutect2'),
                    make_option(c("--signature"), type = "character", default = 'SigProfiler'),
                    make_option(c("--tool"), type = "character", default = 'viber')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
SPN=opt$spn_id
cov=as.numeric(opt$coverage)
pur=as.numeric(opt$purity)
vcf_caller=opt$vcf_caller
cna_caller=opt$cna_caller
sig_caller=opt$signature
dec_caller=opt$tool
mut_process = get_mutations(spn=SPN, type="tumour", coverage=cov, purity=pur)

# Get interpreted table
tool=dec_caller
simulation_id = paste0(cov, "x_", pur, "p_", vcf_caller, "_", cna_caller)
basedir = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"

final_table = tryCatch(
  readRDS(file.path(basedir, "validation_subclonal_new/tables_interpreted", paste0(tool, "_", SPN, "_", simulation_id, ".rds"))),
  error = function(e) {
    message("Skipping simulation_id: ", simulation_id,
            " (", e$message, ")")
    return(NULL)
  }
)


mut_process = get_mutations(spn=SPN, type="tumour", coverage=cov, purity=pur)

# seq_results = readRDS(mut_process)
# 
# process_seq = seq_results %>%
#   mutate(mutation_id=paste0(SPN, ":", chr, ":", chr_pos,  ":", alt))

# saveRDS(object = process_seq, file = '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_tumourevo/process_seq.rds')
process_seq <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_tumourevo/process_seq.rds')

colors_cluster_process = c('indianred', 
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
names(colors_cluster_process) = paste0('C',0:15)
colors_cluster_process["Subclonal"] = "#cccccc"


sample_forest = load_sample_forest(get_sample_forest(SPN))
sample_forest_sub <- sample_forest$get_subforest_for(c('SPN07_1.1', 'SPN07_2.1'))
phylo_forest = load_phylogenetic_forest(get_phylo_forest(spn=SPN))

join_table = final_table %>%
  inner_join(process_seq) %>% 
  select(chr, chr_pos, ref, alt, mutation_id, cluster_id_tool, cluster_id_process, vaf_tool,causes)

mutations_with_cell = join_table %>% 
  rowwise() %>%
  mutate(cell_id=phylo_forest$get_first_occurrences(Mutation(
    chr, chr_pos, ref, alt
  ))[[1]]) %>%
  ungroup()

# tool
cells_labels_tool = mutations_with_cell %>% 
  # select(mutation_id, cell_id, cluster_id_tool) %>% 
  select(mutation_id, cell_id, causes) %>% 
  group_by(cell_id) %>% 
  # summarise(label_list=list(cluster_id_tool)) %>% 
  summarise(label_list=list(causes)) %>% 
  rowwise() %>% 
  mutate(label=names(sort(table(label_list[[1]]), decreasing=TRUE))[1]) %>% 
  ungroup() %>% 
  select(-label_list)

# final_labels_tool = sample_forest$get_nodes() %>% as_tibble() %>% 
#   left_join(cells_labels_tool)

final_labels_tool = sample_forest_sub$get_nodes() %>% as_tibble() %>% 
  left_join(cells_labels_tool) %>% 
  filter(label %in% c("SBS11", 'SBS25', 'SBS26'))

if (dec_caller=="pyclonevi"){
  final_labels_tool <- final_labels_tool %>% 
    mutate(label=paste0("C",label))
}






signature <- phylo_forest$get_exposures()
#palette_signature <-  get_colors_for(signature %>% pull(signature) %>% unique())

df_sign <- get_exposure_ends(phylo_forest) %>%
  mutate(time = round(time),
         end_time = round(end_time)-1)
start_time <- 0
end_time <- max(df_sign$end_time)

timeline_sbs11 <- df_sign %>% filter(signature=="SBS11")
timeline_sbs11 <- timeline_sbs11 %>%
  bind_rows(
    data.frame(
      time = 0,
      end_time   = timeline_sbs11$time[1],
      state      = "absent"
    ),
    data.frame(
      time = timeline_sbs11$end_time[1],
      end_time   = end_time,
      state      = "absent"
    )
  )

###
# > timeline
# start_time end_time  state
# 1     0.0000 289.8730 absent
# 2   289.8730 330.5606  SBS31
# 3   330.5606 446.0000 absent

tree_plot <- plot_forest(sample_forest_sub)
max_Y <- max(tree_plot$data$y, na.rm = TRUE)
phylo_spn07 <- plot_sticks(sample_forest_sub, labels=final_labels_tool,cls=sbs_colors) %>%  #, cls = colors_cluster_process) %>%
  annotate_forest(sample_forest_sub, drivers=F, samples = T, MRCAs = F) +
  theme(legend.position = 'left') + 
  ggplot2::guides(size = "none",
                  shape = "none",
                  color = ggplot2::guide_legend("Species")) +
  ggplot2::geom_rect(
    data = timeline_sbs11 %>% filter(signature=="SBS11"),
    ggplot2::aes(
      xmin = -Inf,
      xmax = Inf,
      ymin =max_Y-time,
      ymax =  max_Y-end_time,
      fill = signature
    ),alpha=0.2
  ) +
  ggplot2::scale_fill_manual(values = sbs_colors) +
  my_ggplot_theme()
