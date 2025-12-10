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
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/report/plotting/tables.R")

driver_info_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/drivers/"
SPNS <- c("SPN01","SPN02","SPN03","SPN04","SPN06","SPN07")
#SPNS <- c("SPN01","SPN03","SPN04")
simulated_drivers <- lapply(SPNS, function(x){
  t = readRDS(file.path(driver_info_dir,x,"process_drivers.rds"))
}) %>% bind_rows()



COVERAGES <- c("50","100","150")
PURITIES <- c("0.3","0.6","0.9")



WGD_samples <- c("SPN01_1.1","SPN01_1.3","SPN06_3.1","SPN06_3.2")

validation_dir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION"
params_grid = expand.grid(COVERAGES, PURITIES)
colnames(params_grid) = c("coverage", "purity")


df_all_SPN_driver <- list()
for (SPN in SPNS){
  spn <- SPN

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
          # select(mutationID,SPN,VAF,code) %>% 
          # dplyr::rename(VAF_process=VAF)
        simulated_drivers_with_vaf_mutect2 <- simulated_drivers %>% 
          mutate(chr=paste0("chr",chr)) %>% 
          dplyr::rename(to="end") %>% 
          inner_join(muts_mutect2,by=c("chr","from"),suffix = c(".process",".mutect2")) %>% 
          mutate(mutationID=paste(chr,from,ref.process,alt.process,sep = ":"))
          # select(mutationID,SPN,VAF,code) %>% 
          # dplyr::rename(VAF_mutect2=VAF)
        
        cna_process <- readRDS(get_process_cna(spn = spn,sample = sample))
        cna_ascat <- read_ASCAT(spn_id = spn,sample_id = sample,coverage = coverage,purity = purity)
        #ascat_purity <- as.numeric(cna_ascat$purity_ploidy$AberrantCellFraction)
        mut_with_cn_ascat <- map_driver_to_cna_ascat(simulated_drivers_with_vaf_mutect2,cna_ascat$CNA,coverage,purity) %>% 
          mutate(sample=sample)
        
        mut_with_cn_process <- map_driver_to_cna_process(simulated_drivers_with_vaf_process,cna_process,coverage,purity) %>% 
          mutate(sample=sample)
        
        if (nrow(mut_with_cn_process)>1){
          max_ratio_cn <- mut_with_cn_process %>% 
            pull(ratio) %>% max()
          mut_with_cn_process <- mut_with_cn_process %>% 
            filter(ratio==max_ratio_cn)
        }
        
        compare_vaf <- inner_join(mut_with_cn_ascat,mut_with_cn_process,by="mutationID",suffix = c(".ascat",".process")) %>% 
            mutate(CN.ascat=major.ascat+minor.ascat) %>% 
            mutate(CN.process=major.process+minor.process) %>% 
            select(VAF.ascat,VAF.process,CN.ascat,CN.process,code.process,sample_name,SPN.process,purity.process,coverage.process,type.process)

        
        
        
        
        #mut_with_cn <- map_driver_to_cna_process(simulated_drivers_with_vaf_process,cna_process,coverage,purity)
      }) %>% bind_rows()
    }
    all_combinations_SPN <- do.call("rbind",all_metrics_comb)
  })
  df_all_SPN_driver[[spn]] <- do.call("rbind",df)  
}

df_all_combs_SPN_driver <- do.call("rbind",df_all_SPN_driver)

fga_status <- df_all_combs_SPN_cna %>% 
  select(sample,spn,fga_class) %>% 
  unique()

df_all_combs_SPN_driver <- df_all_combs_SPN_driver %>% 
  mutate(delta_vaf=VAF.ascat-VAF.process) %>% 
  mutate(delta_cn=CN.ascat-CN.process) %>% 
  dplyr::rename("sample"=sample_name) %>% 
  left_join(fga_status) %>% 
  separate(col = code.process,into = c("gene","sost"),sep = " ",remove = F) %>% 
  mutate(estimation_class=case_when(abs(delta_cn) > 0.5 ~"incorrect",
                                    TRUE ~ "correct")) %>% 
  mutate(WGD=case_when(sample%in%WGD_samples ~ "WGD",
                       TRUE ~ "not WGD")) 

gene_types <- tibble::tibble(
  gene = c("APC","TP53","PIK3CA","KRAS","PTEN","NF1","ATRX"),
  type = c("TS", "TS", "Oncogene", "Oncogene", "TS", "TS", "TS")
)
interesting_genes_tbl <- df_all_combs_SPN_driver %>% 
  filter(estimation_class == "incorrect") %>%
  distinct(code.process,.keep_all = TRUE) %>%
  select(gene, delta_vaf, delta_cn, fga_class, WGD,CN.ascat) %>% 
  mutate(gene_label=paste0(gene," (cn=",CN.ascat,")")) %>% 
  filter(!is.na(gene))
  #left_join(gene_types, by = "gene") %>%
  #mutate(type = ifelse(is.na(type), "Unknown", type))
  



driver_cna_vaf <- df_all_combs_SPN_driver %>% 
  mutate(fga_class = ifelse(is.na(fga_class), "High FGA", fga_class)) %>% 
  ggplot(aes(x=delta_vaf,y = delta_cn))+
  #geom_density2d(aes(color=purity.process),alpha=0.5)+
  scale_color_manual(values = purity_colors,name="Purity")+
  ggnewscale::new_scale_color()+
  geom_point(aes(color=fga_class,shape=WGD))+
  scale_color_manual(name = "FGA class", values = c("High FGA"="indianred3","Low FGA"="dodgerblue3"),na.value = "indianred3") +
  scale_shape_manual(values = c("not WGD"=16,"WGD"=1)) +
  scale_x_continuous(breaks = scales::pretty_breaks(n=3),limits = c(-0.05,0.05))+
  ggnewscale::new_scale_color()+
  # geom_label_repel(
  #   data = interesting_genes_tbl,
  #   aes(label = gene_label),
  #   size = 3,
  #   max.overlaps = 100
  # )+
  # scale_color_manual(name = "Gene Role", values = c("TS"="#3a5875ff","Oncogene"="#b83d51ff")) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 30, vjust = 1),
        axis.ticks.x =element_text(angle = 90, vjust = 0.5, hjust = 1))+
  labs(shape="CNA type",x="VAF deviation",y="Total CN deviation",)+
  my_ggplot_theme()

