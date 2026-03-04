library(ProCESS)
phylo_forest <- load_phylogenetic_forest("/orfeo/scratch/cdslab/shared/SCOUT/SPN07/process/phylo_forest.sff")

# let us simulate a 30x sequencing of the four sample
seq_results <- simulate_seq(phylo_forest, coverage = 10, chromosomes="chr1")

#saveRDS(seq_results, "/orfeo/scratch/cdslab/shared/SCOUT/SPN07/process/10x")

mutations = seq_results$mutations %>% filter(classes != "germinal") %>%
  mutate(
    class = case_when(
      SPN07_1.1.occurrences > 0 & SPN07_1.2.occurrences == 0 & SPN07_1.3.occurrences == 0 ~ "SPN07_1.1.private",
      SPN07_1.2.occurrences > 0 & SPN07_1.1.occurrences == 0 & SPN07_1.3.occurrences == 0 ~ "SPN07_1.2.private",
      SPN07_1.3.occurrences > 0 & SPN07_1.1.occurrences == 0 & SPN07_1.2.occurrences == 0 ~ "SPN07_1.3.private",
      
      SPN07_2.1.occurrences > 0 & SPN07_2.2.occurrences == 0 & SPN07_1.1.occurrences==0  & SPN07_1.2.occurrences==0 & SPN07_1.3.occurrences == 0 ~ "SPN07_2.1.private",
      SPN07_2.2.occurrences > 0 & SPN07_2.1.occurrences == 0 & SPN07_1.1.occurrences==0 & SPN07_1.2.occurrences==0 & SPN07_1.3.occurrences == 0 ~ "SPN07_2.2.private",
      .default = "shared"
    )
  )

mutations %>% filter(class != "shared",
                     SPN07_1.1.VAF > 0.2 | SPN07_1.2.VAF > 0.2 |SPN07_1.3.VAF > 0.2 |
                     SPN07_2.1.VAF > 0.2 | SPN07_2.2.VAF > 0.2 ) %>% 
  reshape2::melt() %>% group_by(class, causes) %>%
  summarise(n = n()) %>% mutate(type = ifelse(grepl("ID", causes), "ID", "SNV")) %>%
  filter(type == "SNV") %>% ggplot(aes(x=class, y=n, fill=causes)) + geom_col()


