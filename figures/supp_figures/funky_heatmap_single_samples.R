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
my_pu_pl <- make_palette(col = "olivedrab4",10)

wgd_samples <- c("SPN01_1.1","SPN01_1.3")
hypermutant_samples <- c("SPN02_1.1","SPN02_1.2")
polyclonal_samples <- c("SPN01_1.1","SPN03_2.1","SPN03_4.1")

df_all_combs_SPN <- df_all_combs_SPN_somatic

spns <- "SPN03"
data_cna_1 <- df_all_combs_SPN_cna%>% 
  filter(spn%in%spns) %>% 
  mutate(comb=paste0(coverage,"_",true_purity)) %>% 
  mutate(coverage=paste0(coverage,"x")) %>% 
  dplyr::rename(estimated_purity=purity) %>% 
  dplyr::rename(purity=true_purity) %>% 
  mutate(delta_purity=abs(delta_purity)) %>% 
  select(sample,coverage,purity,tool,delta_purity) %>%
  # mutate(fga=round(fga,0)) %>% 
  # mutate(fgs=round(fgs,0)) %>% 
  unique() %>% 
  pivot_wider(
    names_from  = c(tool),
    values_from = delta_purity,
    names_glue  = "delta_purity_{tool}"
  )

data_cna_2 <- df_all_combs_SPN_cna%>% 
  filter(spn%in%spns) %>% 
  mutate(comb=paste0(coverage,"_",true_purity)) %>% 
  mutate(coverage=paste0(coverage,"x")) %>% 
  dplyr::rename(estimated_purity=purity) %>% 
  dplyr::rename(purity=true_purity) %>% 
  mutate(delta_ploidy=abs(delta_ploidy)) %>% 
  select(sample,coverage,purity,tool,delta_ploidy) %>%
  unique() %>% 
  pivot_wider(
    names_from  = c(tool),
    values_from = delta_ploidy,
    names_glue  = "delta_ploidy_{tool}"
  )

data_cna_3 <- df_all_combs_SPN_cna%>% 
  filter(spn%in%spns) %>% 
  mutate(comb=paste0(coverage,"_",true_purity)) %>% 
  mutate(coverage=paste0(coverage,"x")) %>% 
  dplyr::rename(estimated_purity=purity) %>% 
  dplyr::rename(purity=true_purity) %>% 
  mutate(delta_ploidy=abs(delta_ploidy)) %>% 
  select(sample,coverage,purity,tool,recall) %>%
  unique() %>% 
  pivot_wider(
    names_from  = c(tool),
    values_from = recall,
    names_glue  = "recall_{tool}"
  )


data_cna_4 <- df_all_combs_SPN_cna%>% 
  filter(spn%in%spns) %>% 
  mutate(comb=paste0(coverage,"_",true_purity)) %>% 
  mutate(coverage=paste0(coverage,"x")) %>% 
  dplyr::rename(estimated_purity=purity) %>% 
  dplyr::rename(purity=true_purity) %>% 
  mutate(delta_ploidy=abs(delta_ploidy)) %>% 
  select(sample,coverage,purity,tool,precision) %>%
  unique() %>% 
  pivot_wider(
    names_from  = c(tool),
    values_from = precision,
    names_glue  = "precision_{tool}"
  )

data_cna <- inner_join(x = data_cna_1,y = data_cna_2) %>% 
  inner_join(data_cna_3) %>% 
  inner_join(data_cna_4)
data_t <- df_all_combs_SPN %>% 
  filter(spn%in%spns) %>% 
  filter(mut_type=="SNV") %>% 
  mutate(mutation_class=case_when(CCF_bin_class=="Subclonal High CCF" ~ "subclonal",
                                  CCF_bin_class=="Subclonal Low CCF" ~ "neutral",
                                  TRUE ~ "clonal")) %>% 
  mutate(comb=paste0(coverage,"_",purity)) %>% 
  mutate(coverage=paste0(coverage,"x")) %>% 
  # group_by(spn,comb,mutation_class,caller) %>% 
  # mutate(mean_sensitivity_all_samples=mean(mean_sensitivity)) %>% 
  # mutate(mean_FDR=mean(FDR)) %>% 
  select(sample,mean_sensitivity,coverage,purity,mutation_class,caller) %>%
  unique() %>% 
  pivot_wider(
    names_from  = c(mutation_class, caller),
    values_from = mean_sensitivity,
    names_glue  = "mean_sensitivity_{mutation_class}_{caller}"
  )
data_t1 <- df_all_combs_SPN %>% 
  filter(spn%in%spns) %>% 
  filter(mut_type=="SNV") %>% 
  # mutate(coverage=paste0(coverage,"x")) %>% 
  mutate(mutation_class=case_when(CCF_bin_class=="Subclonal High CCF" ~ "subclonal",
                                  CCF_bin_class=="Subclonal Low CCF" ~ "neutral",
                                  TRUE ~ "clonal")) %>% 
  mutate(comb=paste0(coverage,"_",purity)) %>% 
  mutate(coverage=paste0(coverage,"x")) %>% 
  # group_by(spn,comb,caller) %>% 
  # mutate(mean_sensitivity_all_samples=mean(mean_sensitivity)) %>% 
  # mutate(mean_FDR=mean(FDR)) %>%
  select(sample,FDR,coverage,purity,caller) %>% 
  unique() %>% 
  pivot_wider(
    names_from  = c(caller),
    values_from = FDR,
    names_glue  = "FDR_{caller}"
  )
  # ungroup()
data_final <- left_join(data_t,data_t1) %>% 
  inner_join(data_cna)
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
  select(sample,cell_composition,coverage,mean_sensitivity_clonal_strelka,mean_sensitivity_neutral_strelka,mean_sensitivity_subclonal_strelka,mean_sensitivity_clonal_mutect2,
         mean_sensitivity_subclonal_mutect2,mean_sensitivity_neutral_mutect2,mean_sensitivity_clonal_freebayes,
         mean_sensitivity_subclonal_freebayes,mean_sensitivity_neutral_freebayes,FDR_strelka,FDR_mutect2,FDR_freebayes,
         delta_purity_ascat,delta_ploidy_ascat,
         delta_purity_sequenza,delta_ploidy_sequenza,
         delta_purity_battenberg,delta_ploidy_battenberg,
         recall_ascat,recall_sequenza,recall_battenberg,
         precision_ascat,precision_sequenza,precision_battenberg) %>% 
  arrange((sample)) %>% 
  mutate(WGD_status=case_when(sample%in%wgd_samples ~ "x",
                              TRUE ~ NA)) %>% 
  mutate(sample_type=case_when(sample%in%polyclonal_samples ~ "poly",
                              TRUE ~ "mono")) %>% 
  mutate(hypermutant_status=case_when(sample%in%hypermutant_samples ~ "high TMB",
                               TRUE ~ "low TMB"))
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
                 purity=my_pu_pl,ploidy=my_pu_pl,
                 recall = make_palette(col = "royalblue3",10),
                 precision = make_palette(col = "royalblue3",10),
                 fdr="indianred4",
                 purity_ploidy="olivedrab4",
                 recall_precision="royalblue3",
                 mean_fdr=make_palette(col = "indianred4",10),
                 sample_info="forestgreen")

cinfo <- tibble(
  group =c("","sample_info","sample_info","sample_info","combination","combination","neutral","neutral","neutral",rep(x = "subclonal",3),rep(x = "clonal",3),rep("fdr",3),rep("ploidy",3),rep("purity",3),rep("recall",3),rep("precision",3)),
  # group =c("","sample_info","sample_info","combination","combination",rep("somatic",12),rep("cna",12)),
  id = c("id","sample_type","hypermutant_status","WGD_status","cell_composition","coverage","mean_sensitivity_neutral_strelka","mean_sensitivity_neutral_mutect2","mean_sensitivity_neutral_freebayes",
         "mean_sensitivity_subclonal_strelka","mean_sensitivity_subclonal_mutect2","mean_sensitivity_subclonal_freebayes",
         "mean_sensitivity_clonal_strelka","mean_sensitivity_clonal_mutect2","mean_sensitivity_clonal_freebayes",
         "FDR_strelka","FDR_mutect2","FDR_freebayes","delta_ploidy_ascat","delta_ploidy_sequenza","delta_ploidy_battenberg",
         "delta_purity_ascat","delta_purity_sequenza","delta_purity_battenberg",
         "recall_ascat","recall_sequenza","recall_battenberg",
         "precision_ascat","precision_sequenza","precision_battenberg"),
  name = c("","sample_type","hypermutant_status","wgd_status","purity","coverage","strelka","mutect2","freebayes","strelka","mutect2","freebayes","strelka","mutect2","freebayes","strelka","mutect2","freebayes",
           "ascat","sequenza","battenberg","ascat","sequenza","battenberg",
           "ascat","sequenza","battenberg","ascat","sequenza","battenberg"),
  geom = c("text","text","text","text","pie","text",rep("circle",12),rep("bar",12)),
  palette = c(NA,NA,NA,NA,"cell_composition",NA,"mean_sensitivity_neutral","mean_sensitivity_neutral","mean_sensitivity_neutral",
              "mean_sensitivity_subclonal","mean_sensitivity_subclonal","mean_sensitivity_subclonal",
              "mean_sensitivity_clonal","mean_sensitivity_clonal","mean_sensitivity_clonal",
              "mean_fdr","mean_fdr","mean_fdr","ploidy","ploidy","ploidy",
              "purity","purity","purity",
              rep("recall",3),
              rep("precision",3)),
  width = c(rep(3,1),rep(3,5),rep(3,24)),
  size = c(rep(1,1),rep(3,3),rep(NA,1),rep(3,1),rep(NA,24))
  # options = c(rep(x = NA,29))
  # overlay = c(rep(F,1),rep(TRUE,2),rep(FALSE,9))
  
)
# cinfo <- cinfo %>% 
#   add_row(id = "fga", group = "sample_info", name = "", geom = "text", size=3,options = lst(lst(overlay = TRUE)), palette = "black",.before = 4) %>%
#   add_row(id = "fgs", group = "sample_info", name = "", geom = "text", size=3,options = lst(lst(overlay = TRUE)), palette = "black",.before = 5)
# palettes$black <- c(rep("black", 2))
row_info <- data_funky %>%
  mutate(group=sample) %>%
  mutate(level=sample) %>%
  mutate(id=as.character(seq(1,nrow(data_funky)))) %>%
  select(group,level,id)

row_grp <- data_funky %>%
  mutate(group=sample) %>%
  mutate(level=sample) %>%
  # mutate(id=as.character(seq(1,nrow(data_funky)))) %>%
  select(group,level)

cgrp <- tibble(
  # Category=c("Sample","Seq",rep("Somatic mutations",12),rep("CNA",12)),
  c("Sample","Seq","Neutral","Subclonal","Clonal","FDR","Ploidy","Purity","Recall","Precision"),
  # Subcategory=c("Sample","Seq","Neutral","Subclonal","Clonal","FDR","Ploidy","Purity","Recall","Precision"),
  Subcategory=c("\n","\n","CCF < 0.10","0.10 < CCF < 0.95","CCF > 0.95","fdr","Delta ploidy","Delta purity","BP recall","BP precision"),
  group=c("sample_info","combination","neutral","subclonal","clonal","fdr","ploidy","purity","recall","precision"),
  # group=c("sample_info","combination","somatic","somatic","somatic","somatic","cna","cna","cna","cna"),
  # palette=c("sample_info","combination","neutral","neutral","neutral","neutral","ploidy","ploidy","ploidy","ploidy")
  palette=c("sample_info","combination","neutral","subclonal","clonal","fdr","purity_ploidy","purity_ploidy","recall_precision","recall_precision")
)
#funky_heatmap(data_funky, column_info = cinfo,palettes = palettes,scale_column = F)
# funky_heatmap(data_funky, column_info = cinfo,palettes = palettes,
#               scale_column = T,column_groups = cgrp,row_info = row_info,row_groups = row_grp)

min_neutral <-round(min(data_funky[,c("mean_sensitivity_neutral_strelka","mean_sensitivity_neutral_mutect2","mean_sensitivity_neutral_freebayes")]),2)
max_neutral <- round(max(data_funky[,c("mean_sensitivity_neutral_strelka","mean_sensitivity_neutral_mutect2","mean_sensitivity_neutral_freebayes")]),2)
min_subclonal <-round(min(data_funky[,c("mean_sensitivity_subclonal_strelka","mean_sensitivity_subclonal_mutect2","mean_sensitivity_subclonal_freebayes")]),2)
max_subclonal <- round(max(data_funky[,c("mean_sensitivity_subclonal_strelka","mean_sensitivity_subclonal_mutect2","mean_sensitivity_subclonal_freebayes")]),2)
min_clonal <-round(min(data_funky[,c("mean_sensitivity_clonal_strelka","mean_sensitivity_clonal_mutect2","mean_sensitivity_clonal_freebayes")]),2)
max_clonal <- round(max(data_funky[,c("mean_sensitivity_clonal_strelka","mean_sensitivity_clonal_mutect2","mean_sensitivity_clonal_freebayes")]),2)
min_fdr <-round(min(data_funky[,c("FDR_strelka","FDR_mutect2","FDR_freebayes")]),2)
max_fdr <- round(max(data_funky[,c("FDR_strelka","FDR_mutect2","FDR_freebayes")]),2)
max_delta_purity <- round(max(data_funky[,c("delta_purity_ascat","delta_purity_sequenza","delta_purity_battenberg")]),2)
min_delta_purity <- round(min(data_funky[,c("delta_purity_ascat","delta_purity_sequenza","delta_purity_battenberg")]),2)
max_delta_ploidy <- round(max(data_funky[,c("delta_ploidy_ascat","delta_ploidy_sequenza","delta_ploidy_battenberg")]),2)
min_delta_ploidy <- round(min(data_funky[,c("delta_ploidy_ascat","delta_ploidy_sequenza","delta_ploidy_battenberg")]),2)

max_recall <- round(max(data_funky[,c("recall_ascat","recall_sequenza","recall_battenberg")]),2)
min_recall <- round(min(data_funky[,c("recall_ascat","recall_sequenza","recall_battenberg")]),2)
max_precision <- round(max(data_funky[,c("precision_ascat","precision_sequenza","precision_battenberg")]),2)
min_precision <- round(min(data_funky[,c("precision_ascat","precision_sequenza","precision_battenberg")]),2)

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
  ),
  list(
    palette = "purity",
    # geom = "bar",
    title = "delta purity",
    labels = c(min_delta_purity,rep("",6),max_delta_purity)
  ),
  list(
    palette = "ploidy",
    # geom = "bar",
    title = "delta ploidy",
    labels = c(min_delta_ploidy,rep("",6),max_delta_ploidy)
  ),
  list(
    palette = "recall",
    # geom = "bar",
    title = "recall",
    labels = c(min_recall,rep("",6),max_recall)
  ),
  list(
    palette = "precision",
    # geom = "bar",
    title = "precision",
    labels = c(min_precision,rep("",6),max_precision)
  )
)
funky_heatmap(data_funky, column_info = cinfo,palettes = palettes,
              scale_column = T,column_groups = cgrp,row_info = row_info,
              row_groups = row_grp,legends = legends)

