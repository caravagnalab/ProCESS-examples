library(ggplot2)
library(patchwork)
library(ProCESS)
library(dplyr)
library(tidyverse)
library(ggpubr)
library(randnet)
library(scales)
library(ggrepel)
library(RColorBrewer)
library(png)
library(grid)
setwd("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/report/")
source("plotting/signature_ProCESS.R")
source("plotting/plot_genome_wide.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
source("plotting/dynamics_ProCESS.R",local =T)
source("plotting/tables.R", local = knitr::knit_global())
metadata <- read.table(file = "SCOUT_metadata.csv",header = T,sep = "\t")

# - Cartoon with clonal evolution vs Process Muller plot -> OK
# - Longitudinal and multiregion profiling (ex: sampling plots with clone proportions) -> OK
# - Multivariate VAF distribution showing spatial and temporal heterogeneity with annotated drivers that explain the evolutionary history. -> OK
# - Marginal VAF showing the multiplicity of drivers, to identify combination of drivers and LOH.-> OK
# - Context profile to show the presence of therapy associated and hypermutamt mutation signatures in relapse samples. 
# - Compare tails profile of normal and hypermutant clones to show the mutation rate increase.

spn <- "SPN07"
sample_forest <- load_sample_forest(get_sample_forest(spn = spn))
phylo_forest <- load_phylogenetic_forest(get_phylo_forest(spn = spn))

sample_names <- get_sample_names(spn = spn) 
number_of_samples <- length(sample_names)
cna_data <- lapply(sample_names,function(s){
  readRDS(get_process_cna(spn = spn,sample = s))
}) %>% bind_rows()


info_spn <- metadata %>% filter(SPN_ID==spn)
nrow = ceiling(number_of_samples/3)
height = nrow * 3

basedir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
spndir <- file.path(basedir,spn,"process")
setwd(spndir)
sim <- recover_simulation(spn)
plot <- plot_tumour_dynamics(spn, sample_forest)

##### PANEL A #### 
# Tissue at sampling time
samples_plot = plot$plot_sampling
# Pie-chart for counts
info = sim$get_samples_info() ## requested from either the simulation recovery or as saved table
color_map_clones <- get_clone_map(sample_forest)
nodes = sample_forest$get_nodes()
clones = nodes %>% 
  dplyr::filter(!is.na(sample)) %>% 
  dplyr::group_by(sample, mutant) %>% 
  dplyr::pull(mutant) %>% 
  unique()
clones_prop = nodes %>%
  dplyr::filter(!is.na(sample)) %>% 
  dplyr::group_by(sample, mutant) %>% 
  # dplyr::mutate(mutant = gsub(" ", "_", mutant)) %>% 
  dplyr::count(mutant) %>% 
  group_by(sample) %>%
  mutate(proportion = n / sum(n))
sample_composition <- ggplot(clones_prop, aes(x = "", y = proportion, fill = mutant)) +
  geom_col(width = 1, color = "white") +
  scale_fill_manual(values=color_map_clones)+
  coord_polar(theta = "y") +
  facet_wrap(~ sample) +
  theme_void()
clones_colors = c("Clone 2"="#e46020ff",
                  "Clone 3"="#7570b0ff",
                  "Clone 4"="#f42c89ff", 
                  "Clone 5"="#59a532ff",
                  "Clone 6"="#edab30ff")
pie_charts = ggplot(clones_prop,
              aes(fill = mutant, y = proportion, x = "")) + facet_grid(~sample) +
  geom_bar(width = 1, color = 0,stat = "identity") +
  coord_polar("y", start=0) +
  #CNAqc:::my_ggplot_theme() +
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(),
        strip.text.x = element_text(size = 6),strip.text.y = element_text(size = 7))+
  scale_fill_manual(values = clones_colors)

ggsave(samples_plot, filename="/orfeo/scratch/cdslab/antonelloa/Process-figures/PanelA_samples.pdf")
ggsave(pie_charts, filename="/orfeo/scratch/cdslab/antonelloa/Process-figures/PanelA_piecharts.pdf")

##### PANEL B #### 
exposure_table <- phylo_forest$get_exposures()
df_sign <- get_exposure_ends(phylo_forest) %>%
  mutate(start_time = round(time),
         end_time = round(end_time)-1) %>%
  select(start_time, end_time, signature, exposure, type)

plot_signatures <- plot_exposure_evolution(sample_forest, phylo_forest, file.path(spndir,spn))
#plot_signatures$sign_tree
muller = plot_signatures$sign_muller
ggsave(muller, filename="/orfeo/scratch/cdslab/antonelloa/Process-figures/PanelA_muller.pdf")

### PANEL D #### 
# Multivariate VAF
purity=0.9
vcf_caller = "mutect2"
cna_caller = "ascat"
coverage = 100
tool = 'viber'
github_path = '/orfeo/cephfs/scratch/cdslab/erivar00/GitHub/ProCESS-examples/'
main_path = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
save_path = file.path(github_path, "validation/Subclonal_deconvolution/")
source(file.path(save_path, "utils_plots_final.R"))
source(file.path(save_path, "generate_table_main.R"))
simulation_id = paste0(coverage, "x_", purity, "p_", vcf_caller, "_", cna_caller)
final_table = tryCatch(
  readRDS(file.path(main_path, "validation_subclonal/tables_interpreted", paste0(tool, "_", spn, "_", simulation_id, ".rds"))),
  error = function(e) {
    message("Skipping simulation_id: ", simulation_id,
            " (", e$message, ")")
    return(NULL)
  }
)
color_palette_process = c("Other"="#cccccc", "Clonal"="#00a3dfff", clones_colors)
#plot_scatter_process
plot_scatter_process_single = function(table_wide, s1,s2, color_palette, driver){
  if(driver==F){
    ggplot()+
      geom_point(data=table_wide,aes(x=eval(parse(text = s1)),
                                     y = eval(parse(text = s2)),
                                     color = cluster_id_process),
                 size=0.5, alpha =0.2) +
      scale_color_manual(values=color_palette)+
      theme_minimal() +
      labs(
        color = "Cluster",
        x = s1,
        y = s2
      )+
      xlim(0, 1)+
      ylim(0, 1)
  }else{
    ggplot() + 
      geom_point(data = table_wide, 
                 aes(x = .data[[s1]], y =.data[[s2]], col = .data$cluster_id_process), 
                 size = .5, alpha = .2) +
      scale_color_manual(values=color_palette)+
      theme_minimal() +
      labs(
        color = "Cluster",
        x = s1,
        y = s2
      )+
      xlim(0, 1)+
      ylim(0, 1)+
      ggtitle('') +
      xlab(s1) +
      ylab(s2) + 
      geom_point(
        data = subset(table_wide, is_driver_process == TRUE),
        aes(x = .data[[s1]], y = .data[[s2]]),
        color = "black", size = 1.5, shape = 15
      ) +
      # ggrepel::geom_label_repel(
      #   data = subset(table_wide, is_driver_process == TRUE),
      #   aes(x = .data[[s1]], y = .data[[s2]], label =.data$gene, color = .data$cluster_id_process),
      #   #color = 'black',
      #   size = 2,
      #   # nudge_y = 0,
      #   # nudge_x = 0,
      #   show.legend = FALSE,
      #   max.overlaps = Inf,
      #   nudge_y = 0.5,
      #   nudge_x = 0.5
      # )  +
      # geom_label(data = subset(table_wide, is_driver_process == TRUE),
      #            aes(x = .data[[s1]], y = .data[[s2]], label =.data$gene, color = .data$cluster_id_process),
      #            size = 2,
      #            nudge_y = 0.09,
      #            nudge_x = 0.09
      #            )+
      guides(
        color = guide_legend(
          ncol = 1,
          override.aes = list(size = 3, alpha = 1)
        )
      )
  }
}
scatter_process = function(table, sample_names, color_palette_process, driver=T, vertical = F){
  # table = join_table_process
  if(spn=='SPN07'){
    table = table %>%
      distinct(patient_id, sample_id, mutation_id, .keep_all = TRUE)
  }
  table_wide = table %>%
    mutate(gene=driver_label_process) %>% 
    select(patient_id, sample_id, mutation_id, cluster_id_process, vaf_process, is_driver_process,gene) %>%
    pivot_wider(values_from="vaf_process", names_from="sample_id")
  
  table_wide <- table_wide %>%
    filter(!is.na(cluster_id_process))
  
  # table_wide[is.na(table_wide)] = 0.0
  table_wide <- table_wide %>%
    mutate(across(starts_with("Spn"), ~replace_na(., 0.0)))
  
  sample_names = as.character(sample_names)
  cm = combn(sample_names, 2)
  
  # s1 = cm[[1]]
  # s2 = cm[[2]]
  plots <- apply(
    cm,
    2,
    function(w) plot_scatter_process_single(table_wide, s1 = w[1], s2 = w[2], color_palette=color_palette_process, driver)
  )
  
  # if(vertical == F){
  #   if(cm %>% ncol() == 1){
  #     nrows = 1
  #     ncols = 1
  #     
  #   }else{
  #     num_pairs = cm %>% ncol()
  #     ncols = min(3, num_pairs) # max 3 cols
  #     nrows = ceiling(num_pairs / ncols)
  #   }
  # }else{
  #   ncols = 1
  #   nrows = length(plots)
  # }
  # 
  # plot_to_save = ggpubr::ggarrange(
  #   plotlist = plots,
  #   ncol = ncols,
  #   nrow = nrows,
  #   common.legend = T,
  #   legend = "right")
  
  # wrap_plots(plots, guides = 'collect')
  
  #plot_to_save 
  plots
}
#  MSH6 c.1082G>A   p.R361H 2-47799065-G-A
final_table = final_table %>% mutate(
  is_driver_tool = ifelse(grepl("47799065", mutation_id), TRUE, is_driver_tool),
  driver_label_process = ifelse(grepl("47799065", mutation_id), "MSH6_p.R361H", driver_label_process)
)
# scatter_process = plot_scatter_process(
#   final_table %>% mutate(cluster_id_process = ifelse(cluster_id_process=="Subclonal", "Other", cluster_id_process)),
#   paste0("SPN07_",sample_names), 
#   color_palette_process,
#   driver=T, vertical = F)
table = final_table %>% mutate(cluster_id_process = ifelse(cluster_id_process=="Subclonal", "Other", cluster_id_process))
sample_names = paste0("SPN07_",sample_names)
color_palette = color_palette_process
driver=T
vertical = F
if(spn=='SPN07'){
  table = table %>%
    distinct(patient_id, sample_id, mutation_id, .keep_all = TRUE)
}
table_wide = table %>%
  mutate(gene=driver_label_process) %>% 
  select(patient_id, sample_id, mutation_id, cluster_id_process, vaf_process, is_driver_process,gene) %>%
  pivot_wider(values_from="vaf_process", names_from="sample_id")

table_wide <- table_wide %>%
  filter(!is.na(cluster_id_process))

# table_wide[is.na(table_wide)] = 0.0
table_wide <- table_wide %>%
  mutate(across(starts_with("Spn"), ~replace_na(., 0.0)))

sample_names = as.character(sample_names)
cm = combn(sample_names, 2)

# s1 = cm[[1]]
# s2 = cm[[2]]
plots <- apply(
  cm,
  2,
  function(w) plot_scatter_process_single(table_wide, s1 = w[1], s2 = w[2], color_palette=color_palette_process, driver) + 
    theme(legend.position="none", 
          axis.title = element_blank(),
          axis.text = element_blank())
)
#scatter_process %>% ggsave(filename="/orfeo/cephfs/scratch/cdslab/antonelloa/Process-figures/multivariate.pdf", height = 10, width = 10)
plots[[1]] %>% ggsave(filename="/orfeo/cephfs/scratch/cdslab/antonelloa/Process-figures/multivariate_SPN07_1.1_vs_1.2.png",
                      height = 5, width = 5, units = "cm", dpi = 300)
plots[[2]] %>% ggsave(filename="/orfeo/cephfs/scratch/cdslab/antonelloa/Process-figures/multivariate_SPN07_1.1_vs_1.3.png",
                      height = 5, width = 5, units = "cm", dpi = 300)
plots[[3]] %>% ggsave(filename="/orfeo/cephfs/scratch/cdslab/antonelloa/Process-figures/multivariate_SPN07_1.1_vs_2.1.png",
                      height = 5, width = 5, units = "cm", dpi = 300)
plots[[4]] %>% ggsave(filename="/orfeo/cephfs/scratch/cdslab/antonelloa/Process-figures/multivariate_SPN07_1.1_vs_2.2.png",
                      height = 5, width = 5, units = "cm", dpi = 300)
plots[[5]] %>% ggsave(filename="/orfeo/cephfs/scratch/cdslab/antonelloa/Process-figures/multivariate_SPN07_1.2_vs_1.3.png",
                      height = 5, width = 5, units = "cm", dpi = 300)
plots[[6]] %>% ggsave(filename="/orfeo/cephfs/scratch/cdslab/antonelloa/Process-figures/multivariate_SPN07_1.2_vs_2.1.png",
                      height = 5, width = 5, units = "cm", dpi = 300)
plots[[7]] %>% ggsave(filename="/orfeo/cephfs/scratch/cdslab/antonelloa/Process-figures/multivariate_SPN07_1.2_vs_2.2.png",
                      height = 5, width = 5, units = "cm", dpi = 300)
plots[[8]] %>% ggsave(filename="/orfeo/cephfs/scratch/cdslab/antonelloa/Process-figures/multivariate_SPN07_1.3_vs_2.1.png",
                      height = 5, width = 5, units = "cm", dpi = 300)
plots[[9]] %>% ggsave(filename="/orfeo/cephfs/scratch/cdslab/antonelloa/Process-figures/multivariate_SPN07_1.3_vs_2.2.png",
                      height = 5, width = 5, units = "cm", dpi = 300)
plots[[10]] %>% ggsave(filename="/orfeo/cephfs/scratch/cdslab/antonelloa/Process-figures/multivariate_SPN07_2.1_vs_2.2.png",
                       height = 5, width = 5, units = "cm", dpi = 300)

## Univariate VAF
final_table = final_table %>% mutate(cluster_id_process = ifelse(cluster_id_process=="Subclonal", "Other", cluster_id_process)) 
drivers_df = final_table %>% filter(is_driver_tool==T) %>% 
  mutate(y = case_when(
    driver_label_tool == "NF1_p.R192*" & sample_id == "SPN07_SPN07_1.1" ~ 300,
    driver_label_tool == "NF1_p.R192*" & sample_id == "SPN07_SPN07_1.2" ~ 400,
    driver_label_tool == "ZFHX3_p.F1371S" & sample_id == "SPN07_SPN07_2.1" ~ 2000,
    driver_label_tool == "TP53_p.R248W"& sample_id == "SPN07_SPN07_2.2"~ 2000,
    driver_label_tool == "NF1_p.R16H" & sample_id == "SPN07_SPN07_2.2"~ 3500,
    driver_label_tool == "PIK3R1" & sample_id == "SPN07_SPN07_2.2"~ 4500,
    driver_label_tool == "ATRX_p.R907*" & sample_id == "SPN07_SPN07_1.3"~ 500,
    driver_label_tool == "ATRX_p.R907*" & sample_id == "SPN07_SPN07_1.2"~ 600,
    driver_label_tool == "ATRX_p.R907*" & sample_id == "SPN07_SPN07_2.2"~ 1300,
    driver_label_tool == "ATRX_p.R907*" & sample_id == "SPN07_SPN07_2.1"~ 4000,
    driver_label_tool == "MSH6_p.R361H" & sample_id == "SPN07_SPN07_2.1"~ 2000,
    driver_label_tool == "MSH6_p.R361H" & sample_id == "SPN07_SPN07_2.2"~ 3000
  )
  
  )
univariate_vaf = final_table %>% ggplot() +
  geom_histogram(aes(x = vaf_process, fill = cluster_id_process), bins = 100) +
  geom_segment(data = drivers_df, 
               aes(x = vaf_process, xend = vaf_process, 
                   y = 0, yend = y, color = cluster_id_process), linetype = "dashed")+
  geom_label(data = drivers_df, 
             aes(x = vaf_process, label=driver_label_tool, y = y, 
                 color = cluster_id_process), size=2.2)+
  scale_fill_manual(values = color_palette_process)+
  scale_color_manual(values = color_palette_process)+
  facet_wrap(~sample_id, scales = "free")+
  xlim(0,1)+
  CNAqc:::my_ggplot_theme()


univariate_vaf2 = final_table %>% filter(sample_id %in% c("SPN07_SPN07_2.2", "SPN07_SPN07_1.3")) %>% 
  ggplot() +
  geom_histogram(aes(x = vaf_process, fill = cluster_id_process), bins = 100) +
  geom_segment(data = drivers_df %>% filter(sample_id %in% c("SPN07_SPN07_2.2", "SPN07_SPN07_1.3")), 
               aes(x = vaf_process, xend = vaf_process, 
                   y = 0, yend = y, color = cluster_id_process), linetype = "dashed")+
  geom_label(data = drivers_df %>% filter(sample_id %in% c("SPN07_SPN07_2.2", "SPN07_SPN07_1.3")), 
             aes(x = vaf_process, label=driver_label_tool, y = y, 
                 color = cluster_id_process), size=2.2)+
  scale_fill_manual(values = color_palette_process)+
  scale_color_manual(values = color_palette_process)+
  facet_wrap(~sample_id #, #scales = "free"
             )+
  xlim(0,1)+
  CNAqc:::my_ggplot_theme()

ggsave(univariate_vaf, filename="/orfeo/scratch/cdslab/antonelloa/Process-figures/PanelB_univariata.pdf")
ggsave(univariate_vaf2, filename="/orfeo/scratch/cdslab/antonelloa/Process-figures/PanelB_univariata_ridotta.pdf")

## Context Profile + Profile of the signatures 
# BiocManager::install("MutationalPatterns")

## Exposure plot
causes_vector = final_table$causes %>% unique()
sbs_colors = c(
  "SBS1" = "#f86c70ff",
  "SBS3" = "#dd7e3bff",
  "SBS5" = "#bac4dcff",
  "SBS11" = "#c5b385ff",
  "SBS25" = "#a24765ff",
  "SBS26" = "#598777ff"
    )
exposures = final_table %>% filter(causes %in% causes_vector[grepl("SBS", causes_vector)]) %>% 
  group_by(cluster_id_process, causes) %>% summarise(n_muts = n()) %>% ungroup() %>%
  group_by(cluster_id_process) %>%
  mutate(p = n_muts / sum(n_muts)) %>% ggplot() +
  geom_col(aes(x = cluster_id_process, fill = causes, y=p)) +
  CNAqc:::my_ggplot_theme() +
  scale_fill_manual(values = sbs_colors)
ggsave(exposures, filename="/orfeo/scratch/cdslab/antonelloa/Process-figures/PanelE_clones_exposures.pdf")

## Plot signatures
cosmic_signatures <- MutationalPatterns::get_known_signatures()
subs  <- c("C>A","C>G","C>T","T>A","T>C","T>G")
bases <- c("A","C","G","T")
labels <- c()
for (s in subs) {
  for (b5 in bases) {
    for (b3 in bases) {
      labels <- c(labels, paste0(b5, "[", s, "]", b3))
    }
  }
}
sbs= cosmic_signatures[, causes_vector[grepl("SBS", causes_vector)]]
sbs = as.data.frame(sbs)
sbs$context = labels
sbs = sbs %>% reshape2::melt() %>% rowwise() %>%
  mutate(substitution = strsplit(strsplit(context, split="\\[")[[1]][2], split="\\]")[[1]][1],
         new_context = paste0(strsplit(context, split="\\[")[[1]][1], "_",
                              strsplit(context, split="\\]")[[1]][2])
         ) 
data_sbs = sbs %>%
  ggplot() +
  geom_bar(aes(value, x = new_context, fill = variable), stat = 'identity') +
  facet_grid(variable ~ substitution, scales = 'free') +
  my_ggplot_theme() +
  theme(
    axis.text.x = element_text(angle = 90, size = 4)
  ) +
  #scale_fill_manual(values = get_SBS_colors()) +
  guides(fill = 'none')  +
  labs(x = 'Context', y = "")+scale_fill_manual(values = sbs_colors)

ggsave(data_sbs, filename="/orfeo/scratch/cdslab/antonelloa/Process-figures/PanelC_sbs_profiles.pdf")

## Plot contexts
all_mutations = readRDS("/orfeo/scratch/cdslab/antonelloa/Process-figures/mutations_table.rds")
all_mutations = all_mutations #%>% filter(chr == "3") 
all_mutations = all_mutations %>% filter(!is.na(cluster_id_process)) %>% 
  mutate(cluster_id_process = ifelse(cluster_id_process=="Subclonal", "Other",cluster_id_process))%>%
  select(chr, chr_pos, ref, alt, patient_id, cluster_id_process, causes) %>% filter(nchar(ref)==1, nchar(alt)==1)
all_mutations_with_cluster = all_mutations
all_mutations <- all_mutations %>%
  mutate(
    CHROMOSOME = chr, #gsub(chr,pattern = "chr",replacement = ""),
    START = chr_pos,
    END = chr_pos,
    REFERENCE = ref,
    VARIANT = alt,
    SAMPLE = patient_id
  ) %>% filter(!is.na(REFERENCE), !is.na(VARIANT)) %>% 
  dplyr::select(SAMPLE,CHROMOSOME, START, END, REFERENCE, VARIANT)
# sbs_counts = RESOLVE::getSBSCounts(data=unique(all_mutations),
#                                    reference=BSgenome.Hsapiens.NCBI.GRCh38)

#RESOLVE::getSBSCounts
get_context = function (data, reference = NULL) {
  if (is.null(reference) | (!is(reference, "BSgenome"))) {
    stop("The reference genome provided as input needs to be a BSgenome object.")
  }
  data <- as.data.frame(data)
  colnames(data) <- c("sample", "chrom", "start", "end", "ref", 
                      "alt")
  data <- data[which(data[, "start"] == data[, "end"]), , drop = FALSE]
  data <- data[which(as.character(data[, "ref"]) %in% c("A", 
                                                        "C", "G", "T")), , drop = FALSE]
  data <- data[which(as.character(data[, "alt"]) %in% c("A", 
                                                        "C", "G", "T")), , drop = FALSE]
  data <- data[, c("sample", "chrom", "start", "ref", "alt"), 
               drop = FALSE]
  colnames(data) <- c("sample", "chrom", "pos", "ref", "alt")
  data <- unique(data)
  data <- data[order(data[, "sample"], data[, "chrom"], data[, 
                                                             "pos"]), , drop = FALSE]
  data <- GRanges(data$chrom, IRanges(start = (data$pos -1), #- 1
                                      width = 3), ref = DNAStringSet(data$ref), alt = DNAStringSet(data$alt), 
                  sample = data$sample)
  if (length(setdiff(seqnames(data), seqnames(reference))) > 
      0) {
    warning("Check chromosome names, not all match reference genome.")
  }
  data$context <- getSeq(reference, data)
  data %>% as.data.frame() %>% rowwise() %>% 
    mutate(context = paste0(strsplit(context, "")[[1]][1], "_",strsplit(context, "")[[1]][3])) #%>% as.data.frame()
}
library(BSgenome.Hsapiens.NCBI.GRCh38)
all_mutations = get_context(all_mutations, reference=BSgenome.Hsapiens.NCBI.GRCh38)
all_mutations_selected = all_mutations %>% mutate(chr=as.character(seqnames), chr_pos=start) %>% 
  select(chr, chr_pos, ref, alt, context)
all_mutations = left_join(
  all_mutations_selected, 
  all_mutations_with_cluster %>% mutate(chr_pos = chr_pos-1), 
  by=c("chr", "chr_pos","ref","alt") 
  )

all_mutations = all_mutations %>% rowwise() %>% 
  mutate(class = paste0(strsplit(context, "_")[[1]][1], "[", ref, ">", alt, "]", strsplit(context, "_")[[1]][2]))

convert_sbs = function(my_class){
  c1 = strsplit(my_class, "\\[")[[1]][1]
  ref = strsplit(strsplit(my_class, "\\[")[[1]][2], ">")[[1]][1]
  alt = strsplit(strsplit(my_class, ">")[[1]][2], "\\]")[[1]][1]
  c2 = strsplit(my_class, "\\]")[[1]][2]
  
  corr = c("A"="T", "C"="G", "G"="C", "T"="A")
  # paste0(corr[c1], "[", corr[ref], ">", corr[alt], "]", corr[c1])
  #paste0(corr[c1], "[", corr[ref], ">", corr[alt], "]", corr[c2])
  paste0(corr[c2], "[", corr[ref], ">", corr[alt], "]", corr[c1])
  #paste0(c1, "[", corr[ref], ">", corr[alt], "]", c2)
}

all_mutations = all_mutations %>% rowwise() %>% mutate(correct_class = ifelse(!(class %in% labels), convert_sbs(class), class))

all_mutations_summary = all_mutations %>% 
  rowwise() %>%
  mutate(substitution = strsplit(strsplit(correct_class, "\\[")[[1]][2], "\\]")[[1]][1],
         context_corr = paste0(strsplit(correct_class, "\\[")[[1]][1], "_", strsplit(correct_class, "\\]")[[1]][2])
  ) %>% 
  #filter(causes == "SBS1") %>%
  select(substitution, correct_class,cluster_id_process, causes) #%>%
  #group_by(substitution, correct_class,cluster_id_process) %>%
  #summarise(n=n())

context_plot= all_mutations_summary %>%
  ggplot() +
  geom_bar(aes(correct_class, #y=n, 
               fill= causes) #, stat = 'identity'
           ) +
  facet_grid(cluster_id_process ~ substitution, scales = 'free') +
  CNAqc:::my_ggplot_theme() +
  theme(
    axis.text.x = element_text(angle = 90, size = 4)
  ) +
  scale_fill_manual(values = sbs_colors) +
  #guides(fill = 'none')  +
  labs(x = 'Context', y = "")

ggsave(context_plot, filename="/orfeo/scratch/cdslab/antonelloa/Process-figures/PanelD_contexts.pdf")


# cosmic_signatures <- MutationalPatterns::get_known_signatures()
# subs  <- c("C>A","C>G","C>T","T>A","T>C","T>G")
# bases <- c("A","C","G","T")
# labels <- c()
# for (s in subs) {
#   for (b5 in bases) {
#     for (b3 in bases) {
#       labels <- c(labels, paste0(b5, "[", s, "]", b3))
#     }
#   }
# }
# sbs= cosmic_signatures[, causes_vector[grepl("SBS", causes_vector)]]
# sbs = as.data.frame(sbs)
# sbs$context = labels
# sbs %>% reshape2::melt() %>% rowwise() %>%
#   mutate(substitution = strsplit(strsplit(context, split="\\[")[[1]][2], split="\\]")[[1]][1],
#          new_context = paste0(strsplit(context, split="\\[")[[1]][1], "_", 
#                               strsplit(context, split="\\]")[[1]][2])
#          ) %>%
#   ggplot(aes(x = context, y = value)) + geom_col() + 
#   facet_grid(variable~substitution, scales = "free_y")+
#   theme()

# plt_sample <- plot_forest(sample_forest,
#                           color_map = color_map_clones) +
#   theme(legend.box = "vertical",
#         legend.box.just = "left",
#         legend.spacing = unit(0.01, "cm"))
# plt_sample_annotated <- annotate_forest(plt_sample, phylo_forest, samples = TRUE, MRCAs = TRUE, driver = F)


  
  
  
