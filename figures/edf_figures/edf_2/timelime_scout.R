library(dplyr)
library(lubridate)
library(ggplot2)
library(hms)
library(ggh4x)
library(scales)
spns <- paste0('SPN0', 1:7)
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_2/get_sarek_total_time.R")
### Process Time
time_ProCESS <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/benchmark/time_ProCESS.rds")


time_ProCESS_spn <- time_ProCESS %>% 
  filter(sample %in% spns) %>% 
  select(time_sample_forest,time_phylo_forest, sample) %>% 
  pivot_longer(
               cols = c(time_sample_forest, time_phylo_forest),
               names_to = "step",
               values_to = "time") %>% 
  mutate(step = gsub("^time_", "", step)) %>%
  # mutate(step="ProCESS Simulation") %>%
  group_by(step, sample) %>% 
  summarise(time=sum(time)) %>% 
  mutate(pipeline_step ="ProCESS Simulation") %>% 
  mutate(substep=step) %>% 
  mutate(time=hms::as_hms(time))
time_ProCESS_spn <- time_ProCESS_spn[rep(seq_len(nrow(time_ProCESS_spn)), 3), ]
time_ProCESS_spn <- time_ProCESS_spn %>% 
  mutate(coverage = c(rep(x = 50,length(spns)),rep(x = 100,length(spns)),rep(x = 150,length(spns)))) 

### Sequencing Time
time_Build_cohort <- readRDS("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/edf_figures/edf_2/time_sequencing.rds")

time_Build_cohort_spn <- time_Build_cohort %>% 
  filter(SPN %in% spns) %>% 
  filter(step!="merge_rds") %>% 
  group_by(step, SPN) %>% 
  summarize(cpu_time_secs_mean=mean(cpu_time_secs)) %>% 
  mutate(time=hms::as_hms(cpu_time_secs_mean)) %>% 
  mutate(pipeline_step ="ProCESS Build Cohort") %>% 
  mutate(substep=step) %>% 
  select(SPN, time,step,substep,pipeline_step)

time_Build_cohort_spn <- time_Build_cohort_spn[rep(seq_len(nrow(time_Build_cohort_spn)), 3), ]
time_Build_cohort_spn <- time_Build_cohort_spn %>% 
  mutate(coverage = c(rep(x = 50,7),rep(x = 100,7),rep(x = 150,7))) 


time_mergin_spn <- time_Build_cohort %>% 
  filter(SPN %in% spns) %>% 
  filter(step=="merge_rds") %>% 
  arrange(cpu_time_secs, SPN) %>% 
  mutate(coverage = c(rep(x = 50,7*3),rep(x = 100,7*3),rep(x = 150,7*3),rep(x = 200,7*3))) %>% 
  group_by(coverage, SPN) %>% 
  mutate(time=mean(cpu_time_secs)) %>% 
  mutate(time=hms::as_hms(time)) %>% 
  mutate(step=paste0("merging rds ",coverage)) %>% 
  mutate(pipeline_step ="ProCESS Build Cohort") %>% 
  mutate(substep=step) %>% 
  select(time,step,substep,pipeline_step,coverage) %>% 
  distinct()







# Convert duration strings to seconds
desired_order_substep <- rev(c(
  "sample_forest",
  "phylo_forest",
  "sequencing",
  "samtools_merge",
  "samtools_split",
  "fastq",
  "merging rds 50",
  "merging rds 100",
  "merging rds 150",
  "merging rds 200",
  "preprocess",
  "cna calling",
  "snv/indel calling",
  "variant/driver annotation",
  "quality control",
  "signature deconvolution",
  "subclonal deconvolution"
))

desired_order_step <- rev(c(
  "sample_forest",
  "phylo_forest",
  "sequencing",
  "samtools_merge",
  "samtools_split",
  "fastq",
  "merging rds 50",
  "merging rds 100",
  "merging rds 150",
  "merging rds 200",
  "preprocess 50","preprocess 100","preprocess 150",
  "variant_calling 50","variant_calling 100","variant_calling 150",
  "tumourevo 50","tumourevo 100","tumourevo 150"
))



time_sarek <- df_nexflow %>%
  select(step,substep,coverage,mean_duration, SPN) %>% 
  mutate(time=hms::as_hms(mean_duration)) %>% 
  mutate(pipeline_step=case_when(step=="tumourevo" ~ "Tumourevo",
                                 TRUE~"Sarek")) %>% 
  mutate(step=paste(step,coverage,sep=" "))


time_all <- time_Build_cohort_spn %>%
  full_join(time_ProCESS_spn %>% ungroup %>% dplyr::rename(SPN = sample), by = c("SPN", "step","substep","time","pipeline_step","coverage")) %>% 
  full_join(time_mergin_spn %>% ungroup, by = c("SPN", "step","substep","time","pipeline_step","coverage")) %>% 
  full_join(time_sarek, by = c("SPN", "step","substep","time","pipeline_step","coverage")) %>% 
  mutate(step = factor(step, levels = desired_order_step),
         substep = factor(substep, levels = desired_order_substep)) %>%
  mutate(pipeline_step = factor(pipeline_step, levels = c("ProCESS Simulation","ProCESS Build Cohort","Sarek","Tumourevo"))) %>% 
  arrange(desc(substep)) %>%
  arrange(desc(step)) %>%
  group_by(SPN, coverage) %>% 
  mutate(duration = time,
         duration_secs = as.numeric(duration),
         start = hms::as_hms(cumsum(lag(duration_secs, default = 0))),
         end   = hms::as_hms(start + duration_secs)) %>% 
  mutate(substep = case_when(
    str_detect(step, "preprocess")  ~ "preprocessing",
    TRUE ~ substep
  ))


time_all <- time_all %>% 
  mutate(new_substep = case_when(
    substep == "sample_forest" ~ 'ProCESS simulation',
    substep == "phylo_forest" ~ 'ProCESS simulation',
    substep == "sequencing" ~ 'ProCESS sequencing',
    substep == "samtools_merge" ~ 'ProCESS sequencing',
    substep == "samtools_split" ~ 'ProCESS sequencing',
    substep == "samtools_split" ~ 'ProCESS sequencing',
    substep == "fastq" ~ 'Generate FASTQ',
    substep == "merging rds 50" ~ 'Merge RDS',
    substep == "merging rds 100" ~'Merge RDS',
    substep == "merging rds 150" ~ 'Merge RDS',
    substep == "merging rds 200" ~ 'Merge RDS',
    substep == "preprocessing" ~ 'Sarek Preprocessing',
    substep == "cna calling" ~ 'Sarek CN calling',
    substep == "snv/indel calling" ~ 'Sarek SNV/INDEL calling',
    substep =="variant/driver annotation" ~ 'Tumourevo QC/Annotation',
    substep == "quality control" ~ 'Tumourevo QC/Annotation',
    substep == "signature deconvolution" ~ 'Tumourevo Signature Deconvolution',
    substep == "subclonal deconvolution" ~ 'Tumourevo Subclonal Deconvolution'
  ))

color_palette_substep <-c(
    "sample_forest" = "deepskyblue3",
    "phylo_forest"= "deepskyblue3",
    "sequencing" = "goldenrod1",
    "samtools_merge"= "goldenrod2",
    "samtools_split"= "goldenrod3",
    "fastq"= "khaki2",
    "merging rds 50" ="khaki1",
    "merging rds 100"="khaki1",
    "merging rds 150" ="khaki1",
    "merging rds 200"="khaki1",
    "preprocessing" = "darkseagreen1",
    "cna calling" = "darkseagreen3",
    "snv/indel calling" = "darkseagreen4",
    "variant/driver annotation" ="coral",
    "quality control"="coral2",
    "signature deconvolution"="coral3",
    "subclonal deconvolution"="coral4"
)

unique(time_all$new_substep)
color_palette_newsubstep <-c(
  "ProCESS simulation" = "rosybrown2",
  "ProCESS sequencing"  = "palevioletred4",
  "Generate FASTQ"  = "lightgoldenrod",
  "Merge RDS" = "goldenrod2",
  "Sarek Preprocessing" = "paleturquoise3",
  "Sarek CN calling" = "cadetblue",
  "Sarek SNV/INDEL calling"   = "skyblue4",
  "Tumourevo QC/Annotation" ="lightsalmon1",
  "Tumourevo Signature Deconvolution"="salmon2",
  "Tumourevo Subclonal Deconvolution"="indianred"
)

color_palette_pipeline_step <-c(
  "ProCESS Simulation" = "darkseagreen3",
  "ProCESS Build Cohort"= "deepskyblue4",
  "Sarek" = "lemonchiffon1"
)

# breaks_process_simulation <- time_all %>% filter(pipeline_step == "ProCESS Simulation") %>% pull(start) %>% min()%>% as_hms()
# breaks_process_simulation<-c(breaks_process_simulation,time_all %>% filter(pipeline_step == "ProCESS Simulation") %>% pull(end)%>% as_hms())
# 
# breaks_process_build_cohort <- time_all %>% filter(pipeline_step == "ProCESS Build Cohort") %>% pull(start) %>% min()%>% as_hms()
# breaks_process_build_cohort<-c(breaks_process_build_cohort,time_all %>% filter(pipeline_step == "ProCESS Build Cohort") %>% pull(end)%>% as_hms())
# 
# breaks_sarek <- time_all %>% filter(pipeline_step == "Sarek") %>% pull(start) %>% min() %>% as_hms()
# breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="preprocess 50") %>% tail(1) %>% pull(end)%>% as_hms())
# breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="preprocess 100") %>% tail(1) %>% pull(end)%>% as_hms())
# breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="preprocess 150") %>% tail(1) %>% pull(end)%>% as_hms())
# breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="variant_calling 50",substep=="cna calling") %>% tail(1) %>% pull(end)%>% as_hms())
# breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="variant_calling 50",substep=="snv/indel calling") %>% tail(1) %>% pull(end)%>% as_hms())
# breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="variant_calling 100",substep=="cna calling") %>% tail(1) %>% pull(end)%>% as_hms())
# breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="variant_calling 100",substep=="snv/indel calling") %>% tail(1) %>% pull(end)%>% as_hms())
# breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="variant_calling 150",substep=="cna calling") %>% tail(1) %>% pull(end)%>% as_hms())
# breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="variant_calling 150",substep=="snv/indel calling") %>% tail(1) %>% pull(end)%>% as_hms())
# 
# start_sarek <- time_all %>% filter(pipeline_step == "Sarek") %>% pull(start) %>% min() %>% as_hms()
# end_sarek <- time_all %>% filter(pipeline_step == "Sarek") %>% pull(end) %>% max()%>% as_hms()
# breaks_sarek<-c(breaks_sarek,end_sarek)
# 
# breaks_tumourevo <- time_all %>% filter(pipeline_step == "Tumourevo") %>% pull(start) %>% min() %>% as_hms()
# breaks_tumourevo<-c(breaks_tumourevo,time_all %>% filter(pipeline_step == "Tumourevo",step=="tumourevo 50") %>% tail(1) %>% pull(end)%>% as_hms())
# breaks_tumourevo<-c(breaks_tumourevo,time_all %>% filter(pipeline_step == "Tumourevo",step=="tumourevo 100") %>% tail(1) %>% pull(end)%>% as_hms())
# breaks_tumourevo<-c(breaks_tumourevo,time_all %>% filter(pipeline_step == "Tumourevo",step=="tumourevo 150") %>% tail(1) %>% pull(end)%>% as_hms())
# 
# start_tumourevo <- time_all %>% filter(pipeline_step == "Tumourevo") %>% pull(start) %>% min() %>% as_hms()
# end_tumourevo <- time_all %>% filter(pipeline_step == "Tumourevo") %>% pull(end) %>% max()%>% as_hms()
# breaks_tumourevo<-c(breaks_tumourevo,end_tumourevo)
# 
# bk = c(
#   round_hms(breaks_sarek, 1),
#   round_hms(breaks_tumourevo, 1)
# )
# lb = bk

        
### New plot
time_all$end <- hms::trunc_hms(time_all$end, 1)
plt_gannt <- time_all %>% 
  filter(coverage!=200) %>%
  ggplot(aes(x = start, xend = end, y = as.factor(coverage), yend = as.factor(coverage), color = as.factor(new_substep))) +
    geom_segment(size = 6) +
    scale_color_manual('Steps', values =color_palette_newsubstep )+
    labs(
         y = "Coverage",
         x = "Time (H:M:S)") +
    theme_minimal(base_size = 10) +
    # theme(
    #   panel.grid.major.y = element_blank(),
    #   panel.grid.minor.y = element_blank()
    # )+
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.minor.ticks=element_blank(),
          axis.ticks = element_blank()) +
  facet_grid(SPN~.) +
  scale_x_continuous(labels = function(x) format(hms::as_hms(x), "%H:%M:%S"))
plt_gannt
# plt <- ggplot(time_all, aes(x = start, xend = end, y = step, yend = step, color = as.factor(substep))) +
#   geom_segment(size = 4) +
#   scale_color_manual(values =color_palette_substep )+
#   labs(title = "SCOUT Timeline: SPN01",
#        y = "Pipeline Step",
#        x = "Time") +
#   theme_bw(base_size = 10) +
#   theme(
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank()
#   )+
#   # geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=pipeline_step),
#   #           inherit.aes = FALSE)+
#   scale_fill_manual(values = color_palette_pipeline_step)+
#   facet_grid(~pipeline_step, scales = "free_x") +
#   theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 10))
# 
# 
# plt <- ggplot(time_all, aes(x = start, xend = end, y = step, yend = step, color = as.factor(substep))) +
#   geom_segment(size = 4) +
#   scale_color_manual(values =color_palette_substep )+
#   labs(title = paste0("SCOUT Timeline: )",spn),
#        y = "Pipeline Step",
#        x = "Time") +
#   theme_bw(base_size = 10) +
#   theme(
#     panel.grid.major.y = element_blank(),
#     panel.grid.minor.y = element_blank()
#   )+
#   # geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=pipeline_step),
#   #           inherit.aes = FALSE)+
#   scale_fill_manual(values = color_palette_pipeline_step)+
#   facet_grid(coverage~pipeline_step, scales = "free") +
#   theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 10))
# 
# 
# plt +
#   facetted_pos_scales(
#     x = list(
#       "ProCESS Simulation"    = scale_x_time(breaks =breaks_process_simulation , labels = time_format("%H:%M:%S")),
#       "ProCESS Build Cohort"  = scale_x_time(breaks = breaks_process_build_cohort, labels = time_format("%H:%M:%S")),
#       "Sarek"  = scale_x_time(breaks = breaks_sarek %>% round(),limits=(c(start_sarek,end_sarek))),
#       "Tumourevo"  = scale_x_time(breaks = breaks_tumourevo %>% round(),limits=(c(start_tumourevo,end_tumourevo)))
#     )
#   )
#ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/benchmark/plot_all_scout_timeline_spn01.pdf",width = 10,height =5,dpi = 300)
