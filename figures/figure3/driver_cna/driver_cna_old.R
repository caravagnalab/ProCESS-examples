rm(list=ls())
library(ProCESS)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
library(tidyverse)
library(ggnewscale)
library(ggrepel)
library(broom)
library(tidyr)
library(tibble)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/sarek_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/nf-validation/bin/cna/utils.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/report/plotting/tables.R")

driver_info_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/drivers/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN06","SPN07")

simulated_drivers <- lapply(SPNS, function(x){
  t = readRDS(file.path(driver_info_dir,x,"process_drivers.rds"))
}) %>% bind_rows()
  


COVERAGES <- c("50","100","150")
PURITIES <- c("0.3","0.6","0.9")
SPN_colors <-c("SPN01"='steelblue', "SPN02"='seagreen', "SPN03"='goldenrod', 
               "SPN04"='coral', "SPN06"='palevioletred', "SPN07"='indianred3')
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/utils.R")
WGD_samples <- c("SPN01_1.1","SPN01_1.3","SPN06_3.1","SPN06_3.2")

validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION"
params_grid = expand.grid(COVERAGES, PURITIES)
colnames(params_grid) = c("coverage", "purity")


df_all_SPN <- list()
for (SPN in SPNS){
  spn <- SPN
  sample_forest <- load_sample_forest(get_sample_forest(spn = spn))
  phylo_forest <- load_phylogenetic_forest(get_phylo_forest(spn = spn))
  clone_ccf_table <- process_ccf(sample_forest = sample_forest,phylo_forest = phylo_forest)
  
  df = lapply(1:nrow(params_grid), function(i) {
    coverage = params_grid[i,]$coverage
    purity = params_grid[i,]$purity
    comb <- paste0(coverage,"x_",purity,"p")
    samples <- get_sample_names(spn = spn)
    validation_dir_somatic <- file.path(validation_dir,spn,"somatic")
    simulated_drivers <- readRDS(file.path(driver_info_dir,spn,"process_drivers.rds")) %>% 
      filter(type=="SID") %>% 
      dplyr::rename(from=start)
    selected_chroms <- simulated_drivers %>% pull(chr)
    all_metrics_comb <- list()
    samples <- get_sample_names(spn = spn)
    for (sample in samples){
      all_metrics_comb[[sample]] <- lapply(selected_chroms, function(chr){
        muts_process <- readRDS(file.path(validation_dir_somatic,comb,"process",sample,paste0("chr",chr,".rds"))) %>% 
          bind_rows()
        muts_mutect2 <- readRDS(file.path(validation_dir_somatic,comb,"mutect2",sample,paste0("chr",chr,".rds"))) %>% 
          bind_rows() %>% 
          mutate(from=case_when(!ref%in%c("A","C","G","T") ~ from+1,
                 TRUE ~ from)) %>% 
          mutate(to=case_when(!ref%in%c("A","C","G","T") ~ to+1,
                                TRUE ~ to))
        
        simulated_drivers_with_vaf_process <- inner_join(simulated_drivers,muts_process)
        simulated_drivers_with_vaf_mutect2 <- simulated_drivers %>% 
          mutate(chr=paste0("chr",chr)) %>% 
          dplyr::rename(to="end") %>% 
          inner_join(muts_mutect2,by=c("chr","from"))

        cna_process <- readRDS(get_process_cna(spn = spn,sample = sample))
        cna_ascat <- read_ASCAT(spn_id = spn,sample_id = sample,coverage = coverage,purity = purity)
        ascat_purity <- as.numeric(cna_ascat$purity_ploidy$AberrantCellFraction)
        mut_with_cn_ascat <- map_driver_to_cna_ascat(simulated_drivers_with_vaf_mutect2,cna_ascat$CNA,coverage,purity) %>% 
          mutate(sample=sample) %>% 
          mutate(estimated_purity=ascat_purity) %>% 
          inner_join(clone_ccf_table)
        
        mut_with_cn_process <- map_driver_to_cna_process(simulated_drivers_with_vaf_process,cna_process,coverage,purity) %>% 
          mutate(sample=sample) %>% 
          # mutate(estimated_purity=ascat_purity) %>% 
          inner_join(clone_ccf_table) %>% 
          mutate(multiplicity=major)
        expected_VAF_driver_process <- expected_vaf(m = mut_with_cn_process$multiplicity,
                                            ccf = mut_with_cn_process$ccf,
                                            purity = as.numeric(as.character(purity)),
                                            karyotype = mut_with_cn_process$karyotype)
                                          
        expected_VAF_driver <- expected_vaf(m = 1,
                                            ccf = mut_with_cn_ascat$ccf,
                                            purity = cna_ascat$purity_ploidy$AberrantCellFraction,
                                            karyotype = mut_with_cn_ascat$karyotype)
        mut_with_cn_ascat <- mut_with_cn_ascat %>% 
          # mutate(expected_CCF=expected_CCF_driver) %>% 
          mutate(expected_VAF=expected_VAF_driver)
        
        mut_with_cn_process <- mut_with_cn_process %>% 
          # mutate(expected_CCF=expected_CCF_driver) %>% 
          mutate(expected_VAF=expected_VAF_driver_process)


        
        
        #mut_with_cn <- map_driver_to_cna_process(simulated_drivers_with_vaf_process,cna_process,coverage,purity)
      }) %>% bind_rows()
    }
    all_combinations_SPN <- do.call("rbind",all_metrics_comb)
  })
  df_all_SPN[[spn]] <- do.call("rbind",df)  
}

df_all_combs_SPN <- do.call("rbind",df_all_SPN)


df_metrics <- df_all_combs_SPN %>%
  group_by(karyotype, coverage, purity) %>%
  do({
    dat <- .
    fit <- lm(expected_VAF ~ VAF, data = dat)
    x_label <- quantile(dat$VAF, 0.90, na.rm = TRUE)

    y_pred <- predict(fit, newdata = data.frame(VAF = x_label))
    
    tibble(
      rmse = sqrt(mean((dat$VAF - dat$expected_VAF)^2, na.rm = TRUE)),
      r = cor(dat$VAF, dat$expected_VAF, use = "complete.obs"),
      x_lab = x_label,
      y_lab = y_pred + 0.015,   # slight offset above line
      x_arrow = x_label,
      y_arrow = y_pred
    )
  }) %>%
  ungroup() %>%
  mutate(
    label = paste0("RMSE = ", round(rmse, 3), "\nR = ", round(r, 3))
  )


plot_driver <- df_all_combs_SPN %>% 
  # mutate(coverage_corr=paste0(coverage,"x")) %>% 
  ggplot(aes(x = VAF, y = expected_VAF)) +
  geom_point(aes(color = karyotype), alpha = 0.8, size = 2) +
  scale_color_manual(values = CNAqc:::get_karyotypes_colors(karyotypes = c("1:1","2:1","1:0","2:2")), 
                     name = "Karyotype") +
  ggnewscale::new_scale_color() +
  
  geom_smooth(
    aes(fill = as.factor(purity), color = as.factor(purity)),
    method = "lm",
    alpha = 0.2,
    linewidth = 1
  ) +
  scale_fill_manual(values = purity_colors, name = "Simulated purity") +
  scale_color_manual(values = purity_colors, name = "Simulated purity") +
  
  facet_wrap(~ coverage) +
  theme_bw() +
  
  geom_label(
    data = df_metrics,
    aes(x = x_lab, y = y_lab, label = label, fill = as.factor(purity)),
    inherit.aes = FALSE,
    color = "black",
    size = 2,
    label.padding = unit(0.04, "lines"),
    label.r = unit(0.05, "lines"),
    label.size = 0.1,
    show.legend = FALSE,
    alpha=0.4
  )+
  labs(x="Real VAF",y="Expected VAF")

plot_driver

df_all_combs_SPN_clonal <- df_all_combs_SPN %>% 
  filter(type=="clonal") %>% 
  # filter(karyotype%in%c("1:1","1:0")) %>% 
  group_by(purity,karyotype) %>% 
  summarise(mean_VAF=mean(VAF)) %>% 
  mutate(multipliciy=case_when(karyotype%in%c("1:1","1:0")~1,
                               TRUE~2))
vaf_palette <- circlize::colorRamp2(
  c(0,0.25,0.33,0.5,0.66,0.8, 1),          # low, mid, high VAF
  c("#4d88a3ff","#c3e79aff","#e8f69dff","#fefebcff","#fee18bff",
    "#fa865cff","#ce4a58ff"))




purity_vals <- seq(0, 1, by = 0.1)
multiplicity_vals <- c(1,2)

karyotypes_m1 <- c("1:0","1:1","2:1","2:2","2:0")
karyotypes_m2 <- c("2:1","2:2","2:0")

df <- data.frame()

for (p in purity_vals) {
  for (m in multiplicity_vals) {
    if (m == 1) {
      for (k in karyotypes_m1) {
        vaf <- expected_vaf(m = m, ccf = 1, purity = p, karyotype = k)
        df <- rbind(df, data.frame(
          purity = p,
          multiplicity = m,
          karyotype = k,
          vaf = vaf,
          stringsAsFactors = FALSE
        ))
      }
    } else if (m == 2) {
      for (k in karyotypes_m2) {
        vaf <- expected_vaf(m = m, ccf = 1, purity = p, karyotype = k)
        df <- rbind(df, data.frame(
          purity = p,
          multiplicity = m,
          karyotype = k,
          vaf = vaf,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

######### M=2
driver_m2_plt <-df %>% 
  filter(multiplicity==2) %>% 
  ggplot(aes(
    x = as.factor(purity),
    y = karyotype,
    fill = vaf
  )) +
  geom_tile(color = "grey70") +
  geom_text(
    data = df_all_combs_SPN_clonal %>% filter(multipliciy==2),
    aes(
      x = as.factor(purity),
      y = karyotype,
      label = sprintf("%.2f", mean_VAF)
    ),
    size = 3,
    color = "black",
    vjust = 0.4
  ) +
  xlab("")+
  ggtitle(label = "m=2")+
  scale_fill_gradientn(
    colours = c("#4d88a3ff","#c3e79aff","#fefebcff" ,"#fa865cff","#ce4a58ff" ),
    name = "Expected VAF"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank()
  )



######### M=1
driver_m1_plt <- df %>% 
  filter(multiplicity==1) %>% 
  ggplot(aes(
    x = as.factor(purity),
    y = karyotype,
    fill = vaf
  )) +
  geom_tile(color = "grey70") +
  geom_text(
    data = df_all_combs_SPN_clonal %>% filter(multipliciy==1),
    aes(
      x = as.factor(purity),
      y = karyotype,
      label = sprintf("%.2f", mean_VAF)
    ),
    size = 3,
    color = "black",
    vjust = 0.4
  ) +
  xlab("")+
  ggtitle(label = "m=1")+
  scale_fill_gradientn(
    colours = c("#4d88a3ff","#c3e79aff","#fefebcff" ,"#fa865cff","#ce4a58ff" ),
    name = "Expected VAF"
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_blank()
  )
plt_driver_cna = driver_m1_plt+driver_m2_plt+plot_layout(guides="collect",nrow=2) & theme(legend.position = "bottom")
ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/driver_cna.pdf",
       plot = plt,device = "pdf",width = 10,height = 5,dpi = 300)
