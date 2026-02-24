# reticulate::use_python("/orfeo/scratch/area/lvaleriani/myconda/bin")
# reticulate::py_config()
# library(SigProfilerMatrixGeneratorR)
# library(SigProfilerAssignmentR)
# 
# somatic <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/spn_plot/SPN05/somatic.rds')
# muts <- somatic %>%
#   ungroup() %>% 
#   filter(DP != 0) %>% 
#   filter(NV != 0) %>% 
#   select(chr,  from, ref, alt, sample_name) %>% 
#   dplyr::mutate(Project = 'spn', Genome = 'GRCh38', mut_type = 'SNP', Type = 'SOMATIC', ID = sample_name, Sample = sample_name) %>%
#   dplyr::rename(chrom = chr, pos_start = from) %>%
#   filter(ref != alt) %>% 
#   rowwise() %>%
#   dplyr::mutate(pos_end = pos_start + abs(str_count(ref) - str_count(alt))) %>%
#   dplyr::select(Project, Sample, ID, Genome, mut_type, chrom, pos_start, pos_end, ref, alt, Type) %>%
#   distinct()
# 
# dir.create('/orfeo/cephfs/scratch/area/lvaleriani/races/spn_plot/SPN05/muts_all/', showWarnings = F, recursive = T)
# write.table(muts, file = '/orfeo/cephfs/scratch/area/lvaleriani/races/spn_plot/SPN05/muts_all/mutations.txt', quote = F, sep = '\t', row.names = F)
# SigProfilerMatrixGeneratorR(project = 'muts_all', genome = "GRCh38", matrix_path = '/orfeo/cephfs/scratch/area/lvaleriani/races/spn_plot/SPN05/muts_all/', plot=T)

y_breaks <- function(x) {
  x <- x[x != 0]                # exclude 0
  if(length(x) == 0) return(NULL)
  max_val <- max(x, na.rm = TRUE) - 0.5
  half_val <- max_val / 2
  max_val <- round(max_val, 0)
  half_val <- round(half_val, 0)
  return(c(half_val, max_val))
}


data_id <- read.table('/orfeo/cephfs/scratch/area/lvaleriani/races/spn_plot/SPN05/muts_all/output/ID/muts_all.ID83.all', header = T, sep = '\t')
rownames(data_id) <- data_id$MutationType
data_id <- data_id %>% select(!MutationType)
data_id_counts <- t(data_id) %>% as.data.frame()
rownames(data_id_counts) <- sub("^[^_]+_(.*)", "\\1", rownames(data_id_counts))
data_id_counts = data_id_counts[rowSums(data_id_counts) > 0,]

id <- data_id_counts %>%
  as.data.frame() %>%
  rownames_to_column("samples") %>%
  pivot_longer(
    cols = -samples,
    names_to = "features",
    values_to = "value"
  ) %>% 
  group_by(samples) %>%
  mutate(
    percentage = 100 * value / sum(value, na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  mutate(spn = 'SPN05')

df_plot_id <- id %>% 
  bascule:::reformat_contexts(what = "ID") %>%
  separate(col = variant, into = c('N', 'Type', 'Value')) %>% 
  mutate(name = paste(N, Type, sep = ' ')) %>% 
  #separate(samples, c("spn", "sample_id"), sep = "_") %>% 
  mutate(name = ifelse(Value == 'R' & Type == 'Del', '> 1bp Deletion', name)) %>%
  mutate(name = ifelse(Value == 'R' & Type == 'Ins', '> 1bp Insertion', name)) %>% 
  mutate(Value = ifelse(Value == 'R', N, Value)) %>% 
  mutate(name = ifelse(Value == 'M', 'Microhomology', name)) %>% 
  mutate(Value = ifelse(Value == 'M', N, Value)) %>% 
  mutate(name = case_when(
    name == '1 Del' ~ '1bp Deletion',
    name == '1 Ins' ~ '1bp Insertion',.default = name
  ))

id <- read.table('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/GRCh38/indel_signatures.txt', header = T)
rownames(id) <- id$Type
id <- id %>% select(!Type)
id <- t(id) %>% as.data.frame()
rownames(id) <- sub("^[^_]+_(.*)", "\\1", rownames(id))
id = id[rowSums(id) > 0,]

id7 <- id %>%
  as.data.frame() %>%
  rownames_to_column("samples") %>%
  pivot_longer(
    cols = -samples,
    names_to = "features",
    values_to = "value"
  ) %>% 
  group_by(samples) %>%
  mutate(
    percentage = 100 * value / sum(value, na.rm = TRUE)
  ) %>%
  ungroup()  %>% 
  bascule:::reformat_contexts(what = "ID") %>%
  separate(col = variant, into = c('N', 'Type', 'Value')) %>% 
  mutate(name = paste(N, Type, sep = ' ')) %>% 
  mutate(name = ifelse(Value == 'R' & Type == 'Del', '> 1bp Deletion', name)) %>%
  mutate(name = ifelse(Value == 'R' & Type == 'Ins', '> 1bp Insertion', name)) %>% 
  mutate(Value = ifelse(Value == 'R', N, Value)) %>% 
  mutate(name = ifelse(Value == 'M', 'Microhomology', name)) %>% 
  mutate(Value = ifelse(Value == 'M', N, Value)) %>% 
  mutate(name = case_when(
    name == '1 Del' ~ '1bp Deletion',
    name == '1 Ins' ~ '1bp Insertion',.default = name
  )) %>% 
  filter(samples == 'ID7') %>% 
  mutate(spn = '')

df_plot_id <- bind_rows(df_plot_id, id7)

df_plot_id$name <- factor(df_plot_id$name, levels = c('1bp Deletion', '1bp Insertion', '> 1bp Deletion', '> 1bp Insertion', 'Microhomology'))


plt_id <- ggplot(df_plot_id, 
                 aes(x = context, y = percentage)) +
  geom_col() +
  ggh4x::facet_nested(spn + samples~name+Value, 
                      scales="free", 
                      space="free_x") + 
  ylab("Percentage of ID") +
  my_ggplot_theme() +
  scale_y_continuous(breaks = y_breaks) +  # control y ticks
  theme(
    axis.text.x = element_text(angle = 0, hjust = 1,size=6),
    axis.text.y=element_text(size=5),
    panel.spacing.x = unit(0.05, "in"), 
    panel.spacing.y = unit(0.01, "in"),
    strip.text.x = element_text(size = 8, margin = margin(t = 0, b = 0), colour = 'gray20'),
    strip.text.y = element_text(size = 8, margin = margin(l = 0, r = 0), colour = 'gray20'),
    strip.background = element_rect(colour = "gray70", fill = "gainsboro"),
    #panel.background = element_blank(),      # removes grey background
    panel.grid.major = element_blank(),      # removes major grid lines
    panel.grid.minor = element_blank(),      # removes minor grid lines
    #panel.border = element_rect(color = "black", fill = NA)
  )
plt_id



data_sbs <- read.table('/orfeo/cephfs/scratch/area/lvaleriani/races/spn_plot/SPN05/muts_all/output/SBS/muts_all.SBS96.all', header = T, sep = '\t')
rownames(data_sbs) <- data_sbs$MutationType
data_sbs <- data_sbs %>% select(!MutationType)
data_sbs_counts <- t(data_sbs) %>% as.data.frame()
rownames(data_sbs_counts) <- sub("^[^_]+_(.*)", "\\1", rownames(data_sbs_counts))
data_sbs_counts = data_sbs_counts[rowSums(data_sbs_counts) > 0,]

sbs <- data_sbs_counts %>%
  as.data.frame() %>%
  rownames_to_column("samples") %>%
  pivot_longer(
    cols = -samples,
    names_to = "features",
    values_to = "value"
  ) %>% 
  group_by(samples) %>%
  mutate(
    percentage = 100 * value / sum(value, na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  mutate(spn = 'SPN05')

df_plot_sbs <- sbs %>% 
  bascule:::reformat_contexts(what = "SBS") %>%
  separate(variant, into = c("v1", "v2"), sep = ">", remove = FALSE) %>% 
  separate(context, into = c("c1", "c2"), sep = "_", remove = FALSE) %>% 
  mutate(context = paste0(c1, v1, c2)) 


sbs <- read.table('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/GRCh38/SBS_signatures.txt', header = T)
rownames(sbs) <- sbs$Type
sbs <- sbs %>% select(!Type)
sbs <- t(sbs) %>% as.data.frame()
rownames(sbs) <- sub("^[^_]+_(.*)", "\\1", rownames(sbs))
sbs = sbs[rowSums(sbs) > 0,]

sbs3 <- sbs %>%
  as.data.frame() %>%
  rownames_to_column("samples") %>%
  pivot_longer(
    cols = -samples,
    names_to = "features",
    values_to = "value"
  ) %>% 
  group_by(samples) %>%
  mutate(
    percentage = 100 * value / sum(value, na.rm = TRUE)
  ) %>%
  ungroup()  %>% 
  bascule:::reformat_contexts(what = "SBS") %>%
  separate(variant, into = c("v1", "v2"), sep = ">", remove = FALSE) %>% 
  separate(context, into = c("c1", "c2"), sep = "_", remove = FALSE) %>% 
  mutate(context = paste0(c1, v1, c2))  %>% 
  filter(samples == 'SBS3') %>% 
  mutate(spn = '')


df_plot_sbs <- bind_rows(df_plot_sbs, sbs3)

plt_sbs <- ggplot(df_plot_sbs, 
                  aes(x = context, y = percentage)) +
  geom_col() +
  ggh4x::facet_nested(spn + samples~variant, 
                      scales="free", 
                      space="free_x") + 
  ylab("Percentage of SBS") +
  my_ggplot_theme() +
  scale_y_continuous(breaks = y_breaks) +  # control y ticks
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, size=6),
    axis.text.y=element_text(size=5),
    panel.spacing.x = unit(0.05, "in"), 
    panel.spacing.y = unit(0.01, "in"),
    strip.text.x = element_text(size = 8, margin = margin(t = 0, b = 0), colour = 'gray20'),
    strip.text.y = element_text(size = 8, margin = margin(l = 0, r = 0), colour = 'gray20'),
    strip.background = element_rect(colour = "gray70", fill = "gainsboro"),
    #panel.background = element_blank(),      # removes grey background
    panel.grid.major = element_blank(),      # removes major grid lines
    panel.grid.minor = element_blank(),      # removes minor grid lines
    #panel.border = element_rect(color = "black", fill = NA)
  ) 
p <- plt_sbs + plt_id + plot_layout(nrow = 2)
ggsave(filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/spn_plot/SPN05/signature_plot.pdf',
       width = 8, height = 5, units = 'in')

