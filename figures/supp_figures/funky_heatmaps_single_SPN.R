library(funkyheatmap)
library(funkyheatmap)
library(dplyr, warn.conflicts = FALSE)
library(tibble, warn.conflicts = FALSE)
library(purrr, warn.conflicts = FALSE)
library(tidyr)
make_palette <- function(col, n = 100) {
  colorRampPalette(c("white", col))(n)
}
my_neutral <- make_palette(col = "grey80",10)
my_subclonal <- make_palette(col = "purple3",10)
my_clonal <- make_palette(col = "goldenrod3",10)
data_t <- df_all_combs_SPN_somatic %>% 
  filter(spn%in%c("SPN01")) %>%
  filter(mut_type=="SNV") %>% 
  mutate(mutation_class=case_when(CCF_bin_class=="Subclonal High CCF" ~ "subclonal",
                                  CCF_bin_class=="Subclonal Low CCF" ~ "neutral",
                                  TRUE ~ "clonal")) %>% 
  mutate(comb=paste0(coverage,"_",purity)) %>% 
  mutate(coverage=paste0(coverage,"x")) %>% 
  group_by(spn,comb,mutation_class,caller) %>% 
  mutate(mean_sensitivity_all_samples=mean(mean_sensitivity)) %>% 
  # mutate(mean_FDR=mean(FDR)) %>% 
  select(mean_sensitivity_all_samples,coverage,purity,mutation_class,caller) %>%
  unique() %>% 
  pivot_wider(
    names_from  = c(mutation_class, caller),
    values_from = mean_sensitivity_all_samples,
    names_glue  = "mean_sensitivity_{mutation_class}_{caller}"
  ) %>% 
  ungroup()
data_t1 <- df_all_combs_SPN %>% 
  # filter(spn%in%c("SPN01")) %>% 
  filter(mut_type=="SNV") %>% 
  # mutate(coverage=paste0(coverage,"x")) %>% 
  mutate(mutation_class=case_when(CCF_bin_class=="Subclonal High CCF" ~ "subclonal",
                                  CCF_bin_class=="Subclonal Low CCF" ~ "neutral",
                                  TRUE ~ "clonal")) %>% 
  mutate(comb=paste0(coverage,"_",purity)) %>% 
  mutate(coverage=paste0(coverage,"x")) %>% 
  group_by(spn,comb,caller) %>% 
  # mutate(mean_sensitivity_all_samples=mean(mean_sensitivity)) %>% 
  mutate(mean_FDR=mean(FDR)) %>%
  select(mean_FDR,coverage,purity,caller) %>% 
  unique() %>% 
  pivot_wider(
    names_from  = c(caller),
    values_from = mean_FDR,
    names_glue  = "mean_FDR_{caller}"
  ) %>% 
  ungroup()
data_final <- left_join(data_t,data_t1) 
data_funky <- data_final %>%
  mutate(
    # coverage = as.numeric(as.character(coverage)),
    purity = as.numeric(as.character(purity)),
    cell_composition = map(
      purity,
      ~ c(
        tumour_cells = .x,
        normal_cells = 1 - .x
      )
    )
  ) %>%
  # rename(id = sample) %>% 
  select(spn,cell_composition,coverage,mean_sensitivity_clonal_strelka,mean_sensitivity_neutral_strelka,mean_sensitivity_subclonal_strelka,mean_sensitivity_clonal_mutect2,
         mean_sensitivity_subclonal_mutect2,mean_sensitivity_neutral_mutect2,mean_sensitivity_clonal_freebayes,
         mean_sensitivity_subclonal_freebayes,mean_sensitivity_neutral_freebayes,mean_FDR_strelka,mean_FDR_mutect2,mean_FDR_freebayes)
# data_funky$spn <- as.numeric(gsub(x=data_funky$spn,pattern = "SPN0",""))
# SPN_colors <- c("steelblue","seagreen","goldenrod","coral",
#                 "magenta4","palevioletred","indianred3")
# names(SPN_colors) <- seq(1,7)
palettes <- list(spn=SPN_colors,mean_sensitivity_neutral = my_neutral, 
                 coverage = "Blues",
                 combination = "grey",
                 cell_composition = c("tumour_cells"="purple3","normal_cells"="pink2"),
                 mean_sensitivity_subclonal = my_subclonal, mean_sensitivity_clonal = my_clonal,
                 neutral="grey30",subclonal="purple3",clonal="goldenrod3",
                 fdr="indianred4",
                 mean_fdr=make_palette(col = "indianred4",10))

cinfo <- tibble(
  group =c("","combination","combination","neutral","neutral","neutral",rep(x = "subclonal",3),rep(x = "clonal",3),rep("fdr",3)),
  id = c("id","cell_composition","coverage","mean_sensitivity_neutral_strelka","mean_sensitivity_neutral_mutect2","mean_sensitivity_neutral_freebayes",
         "mean_sensitivity_subclonal_strelka","mean_sensitivity_subclonal_mutect2","mean_sensitivity_subclonal_freebayes",
         "mean_sensitivity_clonal_strelka","mean_sensitivity_clonal_mutect2","mean_sensitivity_clonal_freebayes",
         "mean_FDR_strelka","mean_FDR_mutect2","mean_FDR_freebayes"),
  name = c("","purity","coverage","strelka","mutect2","freebayes","strelka","mutect2","freebayes","strelka","mutect2","freebayes","strelka","mutect2","freebayes"),
  geom = c("text","pie","text",rep("circle",12)),
  palette = c(NA,"cell_composition",NA,"mean_sensitivity_neutral","mean_sensitivity_neutral","mean_sensitivity_neutral",
              "mean_sensitivity_subclonal","mean_sensitivity_subclonal","mean_sensitivity_subclonal",
              "mean_sensitivity_clonal","mean_sensitivity_clonal","mean_sensitivity_clonal",
              "mean_fdr","mean_fdr","mean_fdr"),
  width = c(rep(1,3),rep(2,12))
  # overlay = c(rep(F,1),rep(TRUE,2),rep(FALSE,9))

)

row_info <- data_funky %>%
  mutate(group=spn) %>%
  mutate(level=spn) %>%
  mutate(id=as.character(seq(1,nrow(data_funky)))) %>%
  select(group,level,id)

row_grp <- data_funky %>%
  mutate(group=spn) %>%
  mutate(level=spn) %>%
  # mutate(id=as.character(seq(1,nrow(data_funky)))) %>%
  select(group,level)

cgrp <- tibble(
  Category=c("Combination","Neutral","Subclonal","Clonal","FDR"),
  # Subcategory=c("Purity","Neutral","Clonal","Subclonal","Neutral","Clonal","Subclonal","Neutral","Clonal","Subclonal"),
  group=c("combination","neutral","subclonal","clonal","fdr"),
  palette=c(NA,"neutral","subclonal","clonal","fdr")
)
#funky_heatmap(data_funky, column_info = cinfo,palettes = palettes,scale_column = F)
funky_heatmap(data_funky, column_info = cinfo,palettes = palettes,
              scale_column = T,column_groups = cgrp,row_info = row_info,row_groups = row_grp)

min_neutral <-round(min(data_funky[,c("mean_sensitivity_neutral_strelka","mean_sensitivity_neutral_mutect2","mean_sensitivity_neutral_freebayes")]),2)
max_neutral <- round(max(data_funky[,c("mean_sensitivity_neutral_strelka","mean_sensitivity_neutral_mutect2","mean_sensitivity_neutral_freebayes")]),2)
min_subclonal <-round(min(data_funky[,c("mean_sensitivity_subclonal_strelka","mean_sensitivity_subclonal_mutect2","mean_sensitivity_subclonal_freebayes")]),2)
max_subclonal <- round(max(data_funky[,c("mean_sensitivity_subclonal_strelka","mean_sensitivity_subclonal_mutect2","mean_sensitivity_subclonal_freebayes")]),2)
min_clonal <-round(min(data_funky[,c("mean_sensitivity_clonal_strelka","mean_sensitivity_clonal_mutect2","mean_sensitivity_clonal_freebayes")]),2)
max_clonal <- round(max(data_funky[,c("mean_sensitivity_clonal_strelka","mean_sensitivity_clonal_mutect2","mean_sensitivity_clonal_freebayes")]),2)
min_fdr <-round(min(data_funky[,c("mean_FDR_strelka","mean_FDR_mutect2","mean_FDR_freebayes")]),2)
max_fdr <- round(max(data_funky[,c("mean_FDR_strelka","mean_FDR_mutect2","mean_FDR_freebayes")]),2)


legends <- list(
  list(
    palette = "mean_sensitivity_neutral",
    # geom = "bar",
    title = "neutral",
    labels = c(min_neutral,rep("",6),max_neutral)
  ),
  list(
    palette = "mean_sensitivity_subclonal",
    # geom = "bar",
    title = "subclonal",
    labels = c(min_subclonal,rep("",6),max_subclonal)
  ),
  list(
    palette = "mean_sensitivity_clonal",
    # geom = "bar",
    title = "clonal",
    labels = c(min_clonal,rep("",6),max_clonal)
  ),
  list(
    palette = "mean_fdr",
    # geom = "bar",
    title = "fdr",
    labels = c(min_fdr,rep("",6),max_fdr)
  ),
  list(
    palette = "cell_composition",
    geom = "pie",
    title = "cell_composition",
    labels = c("tumour cells","normal cells")
  )
)
pdf("funky_heatmap.pdf",height = 20,width = 16)
funky_heatmap(data_funky, column_info = cinfo,palettes = palettes,
              scale_column = T,column_groups = cgrp,row_info = row_info,row_groups = row_grp,legends = legends)
dev.off()
# 
#               