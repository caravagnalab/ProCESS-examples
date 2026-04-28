library(dplyr)
library(ggplot2)
library(vcfR)
library(tidyverse)
library(reshape2)
library(ggalluvial)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/process_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/getters/tumourevo_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/nf-validation//bin/signatures/utils_getters.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/nf-validation/bin/signatures/utils_validation.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/nf-validation/bin/signatures/utils.R")
spn="SPN07"
coverage=150
purity=0.9
samples <- c("SPN07_1.1","SPN07_2.1")
# combinations <- expand.grid(coverage = coverages, purity = purities)


vcf_caller="mutect2"
cna_caller="ascat"
gender <- get_process_gender(spn = spn)
if (gender=="XY"){
  chromosomes <- c(paste0("chr",seq_along(1:22)),"chrX","chrY")
} else {
  chromosomes <- c(paste0("chr",seq_along(1:22)),"chrX")
}


process_snv <- lapply(chromosomes, function(c){
  #samples <- get_sample_names(spn = spn)
  process_data <- list()
  for (sample in samples){
    process_data_ss <-readRDS(paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION/",
                                     spn,"/somatic/",coverage,"x_",purity,"p/process/",sample,"/",c,".rds"))
    process_data[[sample]] <- process_data_ss$SNV
  }
  process_data_all_samples <- do.call("rbind",process_data)
}) %>% bind_rows()


input_signatures <- readRDS(file=get_tumourevo_qc(spn = spn,coverage = coverage,purity = purity,
                                                  tool = "join_cnaqc",
                                                  vcf_caller = vcf_caller,cna_caller = "ascat")[["multi_cnaqc_ALL_rds"]])
#samples <- get_sample_names(spn = spn)
caller_snv <- lapply(samples, function(sample){
  sample_correct <- paste0(spn,"_",sample)
  ss <- input_signatures[["original_cnaqc_objc"]][[sample_correct]]$mutations %>% filter(VAF!=0) %>% 
    mutate(sample_name=sample)
})  %>% bind_rows()

passed <- caller_snv %>% 
  filter(type=="SNV") %>% 
  mutate(mutation_id=paste(chr,from,ref,alt,sep=":"))
process <- process_snv %>% 
  mutate(chr=paste0("chr",chr)) %>% 
  dplyr::rename(mutation_id=mutationID)
joined <- full_join(passed,process,by=c("mutation_id","sample_name"),suffix = c("_mutect2","_process"))
joined <- joined %>% 
  mutate(status=case_when(is.na(chr_mutect2)~"FN",
                          is.na(chr_process)~"FP",
                          TRUE~"TP"))
signatures_centroid <- joined %>% 
  group_by(causes,sample_name) %>% 
  summarise(vaf_centroid_causes=mean(VAF_process))
performance <- joined %>% 
  group_by(causes,status,sample_name) %>%
  summarise(n=n()) %>% 
  group_by(causes,sample_name) %>% 
  mutate(tot_causes=sum(n)) %>% 
  mutate(percentage_status=n/tot_causes)
all_data <- left_join(performance,signatures_centroid)
subclonal_arch <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure5/signatures_cohort/sample_classification_clonal_heterogeneity.rds")
n_tot_counts <- joined %>% 
  filter(status=="TP") %>% 
  filter(str_detect(pattern = "SBS",causes)) %>% 
  group_by(sample_name,causes) %>% 
  summarise(n_counts=n()) %>% 
  mutate(label_c=paste0(causes," mutations = ",n_counts))
pp1 <- joined %>% 
  filter(sample_name%in%samples) %>% 
  filter(status == "TP") %>% 
  ggplot(aes(x = VAF_mutect2, fill = causes)) +
  geom_histogram(binwidth = 0.01) +
  gghighlight::gghighlight(
    causes%in%c("SBS11"),
    unhighlighted_colour = "gainsboro",
    use_group_by = FALSE,
    calculate_per_facet = TRUE
  ) +
  # geom_text(
  #   data = n_tot_counts %>% filter(causes== "SBS5"),
  #   mapping = aes(label = n_counts),
  #   x = Inf,
  #   y = Inf,
  #   hjust = 1.1,
  #   vjust = 1.5,
  #   inherit.aes = FALSE
  # ) +
  facet_wrap(~sample_name, nrow = 1) +
  scale_fill_manual(values = sbs_colors) +
  coord_cartesian(clip = "off") +
  my_ggplot_theme()+
  theme(strip.text.x.top = element_blank())
samples <- c("SPN07_1.1","SPN07_2.1")
pp2 <- subclonal_arch %>% 
  filter(spn=="SPN07") %>% 
  filter(sample%in%samples) %>% 
  filter(coverage==150,
         purity==0.9,
         context=="SBS96") %>% 
  ggplot(aes(x=exposure,y=cluster_id_process,fill=causes)) +
  geom_col()+
  scale_fill_manual(values=sbs_colors)+
  facet_wrap(~sample,nrow = 1)+
  ylab("")+
  xlab("Exposure")+
  my_ggplot_theme()



message("Reading inferred signatures")                                    
all_exposures <- list()
contexts_all <- c("SBS96")
for (context in contexts_all){
  context_classes <- gsub('[[:digit:]]+', '', context)
  
  process <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/signature_exposures/",
                            spn,"/signature_exposure_",coverage,"x_",purity,"p_",context_classes,".rds"))
  samples <- get_sample_names(spn)
  colnames(process) <- samples
  process_gt <- as.data.frame(t(process))
  result <- tryCatch({
    
    # Get SigProfiler paths
    sigprofiler <- get_tumourevo_signatures(
      spn = spn,
      coverage = coverage,
      purity = purity,
      vcf_caller = vcf_caller,
      cna_caller = cna_caller,
      tool = "sigprofiler",
      context = context
    )
    
    bascule <- get_tumourevo_signatures(
      spn = spn,
      coverage = coverage,
      purity = purity,
      vcf_caller = vcf_caller,
      cna_caller = cna_caller,
      tool = "BASCULE",
      context = context_classes
    )
    # Combine and load paths
    paths <- c(
      sigprofiler$COSMIC_exposure,
      sigprofiler$COSMIC_signatures,
      sigprofiler$denovo_exposure,
      sigprofiler$denovo_signatures,
      bascule$refined_fit
    )
    
    data <- load_signature_data(paths)
    
    names(data) <- c(
      "SigProfiler_COSMIC_exposure",
      "SigProfiler_COSMIC_signatures",
      "SigProfiler_denovo_exposure",
      "SigProfiler_denovo_signatures",
      "BASCULE_refined_fit"
    )
    data
    
  })
  
  if (!is.null(result)) {
    tumourevo_signature_res <- result
  }
  
  
  sigprof_aligned <- align_callers(tumourevo_signature_res = tumourevo_signature_res,tool = "SigProfiler",spn = spn) %>% as.data.frame()
  bascule_aligned <- align_callers(tumourevo_signature_res = tumourevo_signature_res,tool = "BASCULE",spn = spn) %>% as.data.frame()
  
  all_exposures[[context]] <- bind_rows(
    process_gt %>% mutate(sample=row.names(.), Method="ProCESS"),
    bascule_aligned %>% mutate(sample=row.names(.), Method="BASCULE"),
    sigprof_aligned %>% mutate(sample=row.names(.), Method="SigProfiler")
  ) %>%
    pivot_longer(cols = starts_with(context_classes),
                 names_to = "Signature",
                 values_to = "Exposure") %>% 
    mutate(context=context)
}
signatures_found <- do.call("rbind",all_exposures)
samples <- c("SPN07_1.1","SPN07_2.1")


p_sankey <- signatures_found %>% 
  filter(sample%in%samples) %>% 
  mutate(Exposure=case_when(Exposure==0 ~NA, 
                            TRUE~Exposure)) %>% 
  ggplot(aes(x = Method, y = Exposure,
             stratum = Signature, alluvium = Signature,
             fill = Signature, label = Signature)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", color = "black") +
  geom_stratum(color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~sample,nrow=1) +
  labs(y = "Exposure", x = "Method") +
  theme_minimal() +
  scale_fill_manual(values=c(sbs_colors,id_colors))+
  my_ggplot_theme()+
  theme(legend.position = "bottom")

pp3 <- all_data %>% 
  filter(str_detect(causes, "SBS")) %>%
  filter(!str_detect(causes, "errors")) %>%
  filter(status=="TP") %>% 
  filter(sample_name%in%samples) %>% 
  ggplot(aes(x=causes,y=n,fill=causes)) +
  geom_col()+
  facet_wrap(~sample_name,nrow = 1)+
  scale_fill_manual(values=sbs_colors)+
  # geom_hline(yintercept = 250,linetype = "dashed")+
  my_ggplot_theme()+
  theme(strip.text.x.top = element_blank())
spn07 <- wrap_plots(list(pp2,pp1,pp3,p_sankey+
                           theme(strip.text.x = element_blank(),
                                 strip.text.y = element_blank())),
                    design = "AAAA\nBBBB\nCCCC\nDDDD", guides="collect") & theme(legend.position = "bottom")
spn07
