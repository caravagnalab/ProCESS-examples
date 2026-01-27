library(dplyr)
library(lubridate)
library(ggplot2)
library(hms)
library(ggh4x)
library(scales)
spn="SPN01"

### Process Time
time_ProCESS <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/benchmark/time_ProCESS.rds")


time_ProCESS_spn <-time_ProCESS %>% 
  filter(sample==spn) %>% 
  select(time_sample_forest,time_phylo_forest) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "step",
    values_to = "time") %>% 
  mutate(step = gsub("^time_", "", step)) %>%
  # mutate(step="ProCESS Simulation") %>%
  group_by(step) %>% 
  summarise(time=sum(time)) %>% 
  mutate(pipeline_step ="ProCESS Simulation") %>% 
  mutate(substep=step) %>% 
  mutate(time=hms::as_hms(time)) 

### Sequencing Time
time_Build_cohort <- readRDS("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/benchmark/time_sequencing_partial.rds")
time_Build_cohort_spn <- time_Build_cohort %>% 
  filter(SPN==spn) %>% 
  filter(step!="merge_rds") %>% 
  group_by(step) %>% 
  summarize(cpu_time_secs_mean=mean(cpu_time_secs)) %>% 
  mutate(time=hms::as_hms(cpu_time_secs_mean)) %>% 
  mutate(pipeline_step ="ProCESS Build Cohort") %>% 
  mutate(substep=step) %>% 
  select(time,step,substep,pipeline_step)


time_mergin_spn <- time_Build_cohort %>% 
  filter(SPN==spn) %>% 
  filter(step=="merge_rds") %>% 
  arrange(cpu_time_secs) %>% 
  mutate(coverage = c(rep(x = 50,3),rep(x = 100,3),rep(x = 150,3),rep(x = 200,3))) %>% 
  group_by(coverage) %>% 
  summarise(time=mean(cpu_time_secs)) %>% 
  mutate(time=hms::as_hms(time)) %>% 
  mutate(step=paste0("merging rds ",coverage)) %>% 
  mutate(pipeline_step ="ProCESS Build Cohort") %>% 
  mutate(substep=step) %>% 
  select(time,step,substep,pipeline_step)







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
  select(step,substep,coverage,mean_duration) %>% 
  mutate(time=hms::as_hms(mean_duration)) %>% 
  mutate(pipeline_step=case_when(step=="tumourevo" ~ "Tumourevo",
                                 TRUE~"Sarek")) %>% 
  mutate(step=paste(step,coverage,sep=" "))


time_all <- time_Build_cohort_spn %>%
  full_join(time_ProCESS_spn, by = c("step","substep","time","pipeline_step")) %>% 
  full_join(time_mergin_spn, by = c("step","substep","time","pipeline_step")) %>% 
  full_join(time_sarek, by = c("step","substep","time","pipeline_step")) %>% 
  mutate(step = factor(step, levels = desired_order_step),
         substep = factor(substep, levels = desired_order_substep)) %>%
  mutate(pipeline_step = factor(pipeline_step, levels = c("ProCESS Simulation","ProCESS Build Cohort","Sarek","Tumourevo"))) %>% 
  arrange(desc(substep)) %>%
  arrange(desc(step)) %>%
  mutate(duration = time,
         duration_secs = as.numeric(duration),
         start = hms::as_hms(cumsum(lag(duration_secs, default = 0))),
         end   = hms::as_hms(start + duration_secs))


color_palette_substep <-c(
  "sample_forest" = "deepskyblue3",
  "phylo_forest"= "deepskyblue3",
  "sequencing" = "goldenrod",
  "samtools_merge"= "goldenrod",
  "samtools_split"= "goldenrod",
  "fastq"= "goldenrod",
  "merging rds 50" ="goldenrod",
  "merging rds 100"="goldenrod",
  "merging rds 150" ="goldenrod",
  "merging rds 200"="goldenrod",
  "preprocess" = "darkseagreen1",
  "cna calling" = "darkseagreen3",
  "snv/indel calling" = "darkseagreen4",
  "variant/driver annotation" ="coral",
  "quality control"="coral2",
  "signature deconvolution"="coral3",
  "subclonal deconvolution"="coral4"
)

color_palette_pipeline_step <-c(
  "ProCESS Simulation" = "darkseagreen3",
  "ProCESS Build Cohort"= "deepskyblue4",
  "Sarek" = "lemonchiffon1"
)

breaks_process_simulation <- time_all %>% filter(pipeline_step == "ProCESS Simulation") %>% pull(start) %>% min()%>% as_hms()
breaks_process_simulation<-c(breaks_process_simulation,time_all %>% filter(pipeline_step == "ProCESS Simulation") %>% pull(end)%>% as_hms())

breaks_process_build_cohort <- time_all %>% filter(pipeline_step == "ProCESS Build Cohort") %>% pull(start) %>% min()%>% as_hms()
breaks_process_build_cohort<-c(breaks_process_build_cohort,time_all %>% filter(pipeline_step == "ProCESS Build Cohort") %>% pull(end)%>% as_hms())

breaks_sarek <- time_all %>% filter(pipeline_step == "Sarek") %>% pull(start) %>% min() %>% as_hms()
breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="preprocess 50") %>% tail(1) %>% pull(end)%>% as_hms())
breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="preprocess 100") %>% tail(1) %>% pull(end)%>% as_hms())
breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="preprocess 150") %>% tail(1) %>% pull(end)%>% as_hms())
breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="variant_calling 50",substep=="cna calling") %>% tail(1) %>% pull(end)%>% as_hms())
breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="variant_calling 50",substep=="snv/indel calling") %>% tail(1) %>% pull(end)%>% as_hms())
breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="variant_calling 100",substep=="cna calling") %>% tail(1) %>% pull(end)%>% as_hms())
breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="variant_calling 100",substep=="snv/indel calling") %>% tail(1) %>% pull(end)%>% as_hms())
breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="variant_calling 150",substep=="cna calling") %>% tail(1) %>% pull(end)%>% as_hms())
breaks_sarek<-c(breaks_sarek,time_all %>% filter(pipeline_step == "Sarek",step=="variant_calling 150",substep=="snv/indel calling") %>% tail(1) %>% pull(end)%>% as_hms())

start_sarek <- time_all %>% filter(pipeline_step == "Sarek") %>% pull(start) %>% min() %>% as_hms()
end_sarek <- time_all %>% filter(pipeline_step == "Sarek") %>% pull(end) %>% max()%>% as_hms()
breaks_sarek<-c(breaks_sarek,end_sarek)

breaks_tumourevo <- time_all %>% filter(pipeline_step == "Tumourevo") %>% pull(start) %>% min() %>% as_hms()
breaks_tumourevo<-c(breaks_tumourevo,time_all %>% filter(pipeline_step == "Tumourevo",step=="tumourevo 50") %>% tail(1) %>% pull(end)%>% as_hms())
breaks_tumourevo<-c(breaks_tumourevo,time_all %>% filter(pipeline_step == "Tumourevo",step=="tumourevo 100") %>% tail(1) %>% pull(end)%>% as_hms())
breaks_tumourevo<-c(breaks_tumourevo,time_all %>% filter(pipeline_step == "Tumourevo",step=="tumourevo 150") %>% tail(1) %>% pull(end)%>% as_hms())

start_tumourevo <- time_all %>% filter(pipeline_step == "Tumourevo") %>% pull(start) %>% min() %>% as_hms()
end_tumourevo <- time_all %>% filter(pipeline_step == "Tumourevo") %>% pull(end) %>% max()%>% as_hms()
breaks_tumourevo<-c(breaks_tumourevo,end_tumourevo)


plt <- ggplot(time_all, aes(x = start, xend = end, y = step, yend = step, color = as.factor(substep))) +
  geom_segment(size = 4) +
  scale_color_manual(values =color_palette_substep )+
  labs(title = "SCOUT Timeline: SPN01",
       y = "Pipeline Step",
       x = "Time") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank()
  )+
  # geom_rect(aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=pipeline_step),
  #           inherit.aes = FALSE)+
  scale_fill_manual(values = color_palette_pipeline_step)+
  facet_grid(~pipeline_step, scales = "free_x") +
  theme(legend.position = "bottom",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 10))

plt +
  facetted_pos_scales(
    x = list(
      "ProCESS Simulation"    = scale_x_time(breaks =breaks_process_simulation , labels = time_format("%H:%M:%S")),
      "ProCESS Build Cohort"  = scale_x_time(breaks = breaks_process_build_cohort, labels = time_format("%H:%M:%S")),
      "Sarek"  = scale_x_time(breaks = breaks_sarek %>% round(),limits=(c(start_sarek,end_sarek))),
      "Tumourevo"  = scale_x_time(breaks = breaks_tumourevo %>% round(),limits=(c(start_tumourevo,end_tumourevo)))
    )
  )
ggsave(filename = "/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/benchmark/plot_all_scout_timeline.pdf",width = 12,height =4,dpi = 300)
