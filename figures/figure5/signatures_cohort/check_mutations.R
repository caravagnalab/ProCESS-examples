library(dplyr)
library(ggplot2)
library(vcfR)
library(tidyverse)
library(naniar)
library(MutationalPatterns)
library(reshape2)
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
  samples <- get_sample_names(spn = spn)
  process_data <- list()
  for (sample in samples){
    process_data_ss <-readRDS(paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION/",
                                     spn,"/somatic/",coverage,"x_",purity,"p/process/",sample,"/",c,".rds"))
    process_data[[sample]] <- process_data_ss$SNV
  }
  process_data_all_samples <- do.call("rbind",process_data)
}) %>% bind_rows()


input_signatures <- readRDS(file=get_tumourevo_qc(spn = spn,coverage = coverage,purity = purity,
                                                  tool = "join_cnaqc",vcf_caller = vcf_caller,cna_caller = "ascat")[["multi_cnaqc_ALL_rds"]])
samples <- get_sample_names(spn = spn)
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
expected_sig <- all_data %>% 
  filter(str_detect(causes, "SBS")) %>%
  filter(!str_detect(causes, "errors"))  %>% pull(causes) %>% unique()
signatures_bench_out <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION/",spn,"/signature/",coverage,"x_",purity,"/mutect2_ascat/metrics_SBS96_sample.rds"))



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


p_sankey <- ggplot(signatures_found,
                   aes(x = Method, y = Exposure,
                       stratum = Signature, alluvium = Signature,
                       fill = Signature, label = Signature)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", color = "black") +
  geom_stratum(color = "black") +
  scale_y_continuous(expand = c(0,0)) +
  facet_grid(context~sample) +
  labs(y = "Exposure", x = "Method") +
  theme_minimal() +
  scale_fill_manual(values=c(sbs_colors,id_colors))+
  my_ggplot_theme()+
  theme(legend.position = "bottom")

# signatures_bench_out_sankey <- readRDS(paste0("/orfeo/cephfs/scratch/cdslab/shared/SCOUT/VALIDATION/",spn,"/signature/",coverage,"x_",purity,"/mutect2_ascat/sankey_plot_SBS96.rds"))

recall_plot <- signatures_bench_out %>% 
  distinct() %>% 
  ggplot(aes(x=sample,y=recall,fill=caller,group=caller))+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c(
    'SigProfiler' = 'sienna1',
    'BASCULE'     = 'dodgerblue4'
  ))+
  my_ggplot_theme()

precision_plot <- signatures_bench_out %>% 
  distinct() %>% 
  ggplot(aes(x=sample,y=precision,fill=caller,group=caller))+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c(
    'SigProfiler' = 'sienna1',
    'BASCULE'     = 'dodgerblue4'
  ))+
  my_ggplot_theme()

p1 <- all_data %>% 
  filter(str_detect(causes, "SBS")) %>%
  filter(!str_detect(causes, "errors")) %>%
  ggplot(aes(x=causes,y=percentage_status,fill=status))+
  geom_col()+
  facet_wrap(~sample_name,nrow = 1,scales="free")+
  scale_fill_manual(values=c("FN"="firebrick3","TP"="forestgreen"))+
  my_ggplot_theme()+
  ggtitle(label = spn,subtitle = paste0("Coverage: ",coverage, "\nPurity: ",purity))


p2 <- all_data %>% 
  filter(str_detect(causes, "SBS")) %>%
  filter(!str_detect(causes, "errors")) %>%
  ggplot(aes(x=sample_name,y=causes,size=vaf_centroid_causes))+
  geom_point()+
  my_ggplot_theme()


p3 <- all_data %>% 
  filter(str_detect(causes, "SBS")) %>%
  filter(!str_detect(causes, "errors")) %>%
  filter(status=="TP") %>% 
  ggplot(aes(x=causes,y=n,fill=causes)) +
  geom_col()+
  facet_wrap(~sample_name,nrow = 1,scales="free")+
  scale_fill_manual(values=sbs_colors)+
  geom_hline(yintercept = 250,linetype = "dashed")+
  my_ggplot_theme()

cosmic_path <- '/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples-my/nf-validation/bin/data/COSMIC_v3.4_SBS_GRCh38.txt'
mat = read.table(file = cosmic_path,header = T,row.names = 1)
# Normalize columns
norm_mat <- apply(mat, 2, function(x) x / sqrt(sum(x^2)))

# Cosine similarity = dot product of normalized vectors
cos_sim_matrix <- t(norm_mat) %*% norm_mat
df <- melt(cos_sim_matrix[expected_sig,
                          expected_sig])

cos_sim_sig <- ggplot(df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  my_ggplot_theme()+
  xlab("")+
  ylab("")



wrap_plots(list(p1,p2,p3,recall_plot,precision_plot,cos_sim_sig,p_sankey),design = "AAAA\nBBBB\nCCCC\nDEFF\nGGGG", guides="collect") & theme(legend.position = "bottom")

ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure5/signatures_cohort/test.pdf",width = 10,height = 12)

