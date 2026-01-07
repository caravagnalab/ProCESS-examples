library(ComplexUpset)
library(ggplot2)
library(dplyr)

df_battenberg_all_proc <- df_battenberg_all %>% 
  #mutate(id = paste(chr, from, to, sample_id, sep = ':')) %>% 
  mutate(is_match = ifelse(is_match == 'match subclone', 'no match', is_match)) %>% 
  mutate(is_match = ifelse(is_match == 'switch subclone', 'no match', is_match)) %>% 
  mutate(Battenberg = ifelse(is_match == 'match', TRUE, FALSE)) %>% 
  select(chr, from, to, sample_id, type, purity, coverage, Battenberg)

df_ascat_all_proc <- df_ascat_all %>% 
  #mutate(id = paste(chr, from, to, sample_id, sep = ':')) %>% 
  mutate(ASCAT = ifelse(is_match == 'match', TRUE, FALSE))  %>% 
  select(chr, from, to, sample_id, type, purity, coverage, ASCAT)

df_sequenza_all_proc <- df_sequenza_all %>% 
  #mutate(id = paste(chr, from, to, sample_id, sep = ':')) %>% 
  mutate(Sequenza = ifelse(is_match == 'match', TRUE, FALSE))  %>% 
  select(chr, from, to, sample_id, type, purity, coverage, Sequenza)
  
df <- full_join(df_battenberg_all_proc, df_ascat_all_proc, by = join_by(chr, from, to, sample_id, type,purity, coverage)) %>% 
  full_join(df_sequenza_all_proc, by = join_by(chr, from, to, sample_id, type,purity, coverage)) %>% 
  mutate(ProCESS = TRUE) %>% 
  filter(!is.na(type)) %>% 
  filter(type!= 'subclonal')


df_abs <- relative_to_absolute_coordinates_new(df) 

df_new = df_abs %>% 
  pivot_longer(cols = c(Battenberg, ASCAT, Sequenza, ProCESS)) %>% 
  filter(value == TRUE)

chrom_lengths_abs <- df_new %>%
  mutate(chr = factor(chr, levels = chr_level)) %>%
  arrange(chr) %>%
  group_by(chr) %>%
  summarise(chr_end = max(to), .groups = "drop",
            chr_start = min(from)) %>%
  mutate(
    chr_mid   = (chr_start + chr_end)/2
  ) %>%
  filter(!is.na(chr))

chrom_lengths_abs <- chrom_lengths_abs %>%
  mutate(shade = as.factor(row_number() %% 2))



all_density <- df_new %>% 
  ggplot() + 
  geom_rect(data = chrom_lengths_abs,
            aes(xmin = chr_start, xmax = chr_end , ymin = 0, ymax = 4e-10, fill = shade), alpha = 0.1, inherit.aes = FALSE) +
  scale_fill_manual(values = c("grey50", "white"), guide = "none") +
  geom_vline(xintercept = chrom_lengths_abs$chr_start[-1], color = "grey70", linetype = "dashed") +
    geom_density(data = df_new, mapping = aes(x = from, col = name), position = 'identity', alpha = .2, bw = 1e7) + 
    geom_vline(xintercept = chrom_lengths_abs$chr_start[-1], color = "grey70", linetype = "dashed") +
    scale_color_manual('Tool', values = c('ASCAT' = 'mediumpurple3', 'Sequenza' = 'darkorange2', 'Battenberg' = 'steelblue3', 'ProCESS' = 'gray'))+ 
    scale_x_continuous(
      breaks = chrom_lengths_abs$chr_mid,
      labels = chrom_lengths_abs$chr
    )  +
    my_ggplot_theme()


# run upset
upset <- ComplexUpset::upset(df %>% select(ProCESS, ASCAT, Battenberg, Sequenza, type), c("ASCAT", "ProCESS", 'Sequenza', "Battenberg"), width_ratio = 0.1,  base_annotations=list(
  'Intersection size'=ComplexUpset::intersection_size(
    counts=F,
    mapping=aes(fill=type)
  )
), set_sizes = F)


