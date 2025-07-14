options(bitmapType='cairo')
library(dplyr)
library(ProCESS)
library(optparse)
library(tidyr)
library(ggplot2)
library(future.apply)
library(progressr)
library(patchwork)
source("../../getters/sarek_getters.R")
source("../../getters/process_getters.R")
source("utils.R")

ascat_df <- absolute_to_relative_coordinates(readRDS('data/ascat.rds') %>% dplyr::rename(chr=V1))
cnvkit_df <- absolute_to_relative_coordinates(readRDS('data/cnvkit.rds'), centromere = T) 

############ Parse command-line arguments
option_list <- list(make_option(c("--spn_id"), type = "character", default = 'SPN03'),
                    make_option(c("--purity"), type = "character", default = '0.3'),
                    make_option(c("--coverage"), type = "character", default = '50')
		                )

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
data_dir = '/orfeo/scratch/cdslab/shared/SCOUT/'

spn_id = opt$spn_id
coverage = opt$coverage
purity = opt$purity

samples = get_sample_names(spn_id)
chromosomes = c(paste0('chr',1:22))

process_rds <- readRDS(get_mutations(spn = spn_id, coverage = coverage, purity = purity, type = 'tumour'))
process_normal <- readRDS(get_mutations(spn = spn_id, coverage = coverage, purity = purity, type = 'normal'))
process_normal_long <- ProCESS::seq_to_long(process_normal)
process_normal_long <- process_normal_long %>% ungroup() %>% select(chr, from, to, NV, DP, VAF)
rm(process_normal)
process_rds_long <- ProCESS::seq_to_long(process_rds)
rm(process_rds)
process_rds_long <- process_rds_long %>% ungroup() %>% filter(classes == 'germinal') %>% select(chr, from, to, ref, alt, NV, DP, VAF, sample_name)
cov_T = as.numeric(coverage)
cov_N = 30
ratio = cov_N/cov_T

process_rds_long <- process_rds_long %>% 
  left_join(process_normal_long, by = join_by(chr, from, to), suffix = c('_T', '_N')) %>% 
  filter(ref %in% c('A', 'T', 'C', 'G')) %>% 
  filter(alt %in% c('A', 'T', 'C', 'G')) %>%
  select(chr, from, to, NV_T, DP_T, VAF_T, sample_name, DP_N, VAF_N) %>% 
  mutate(DR = (DP_T/DP_N)*ratio) %>% 
  filter(VAF_N > 0.3, VAF_N < 0.7) 


for (sample_id in samples){
  message(paste0("Reading ProCESS data for sample ", sample_id))
  
  process_rds <- process_rds_long %>% filter(sample_name == sample_id) %>% filter(DP_T > 10) %>% sample_n(1e4)
  process_rds <- absolute_to_relative_coordinates_muts(process_rds %>% mutate(chr = paste0('chr',chr)))
  plt_snp <- CNAqc:::blank_genome(chromosomes = chromosomes) + 
    geom_point(data = process_rds, aes(x = from, y = VAF_T), size = .1) + 
    ylab('BAF') + 
    CNAqc:::blank_genome(chromosomes = chromosomes) + 
    geom_point(data = process_rds, aes(x = from, y = DR), size = .1) + 
    ylab('DR') +
    ylim(-1,3) + 
    plot_layout(nrow=2)
  
  ProCESS_output = read_ProCESS(spn_id,sample_id,coverage,purity)
  CNA_ProCESS = ProCESS_output[["CNA"]] 
  
  CNA_ProCESS = CNA_ProCESS %>% mutate(ratio = ifelse(ratio > 0.09 & ratio < 0.1, 0.1, round(ratio,1))) %>% 
           filter(ratio != 0) 
  
  # CNA_ProCESS <- lapply(chromosomes, FUN = function(c){
  #   df <- CNA_ProCESS %>% filter(chr == c) 
  #   
  #   if (nrow(df %>% filter(ratio !=1))>0){
  #   
  #     df = df %>% 
  #       mutate(ratio = ifelse(ratio > 0.09 & ratio < 0.1, 0.1, round(ratio,1))) %>% 
  #       filter(ratio != 0) %>%
  #       group_by(chr, major, minor, ratio) %>% 
  #       summarise(from = min(from), to = max(to)) %>% 
  #       ungroup()
  #   
  #     breakpoint <- sort(unique(c(df$from, df$to, 1)))
  #     
  #     segments <- tibble(
  #       chr = c,
  #       seg_from = head(breakpoint, -1),
  #       seg_to = tail(breakpoint, -1),
  #       seg_id = paste(chr, seg_from, seg_to, sep = ':')
  #     )
  #     
  #     df <- df %>% 
  #       left_join(segments, by = join_by(chr), relationship = "many-to-many") %>% 
  #       filter((from >= seg_from & to <= seg_to & ratio < 1) | (ratio == 1 & seg_from == from | seg_to == to) | (ratio < 1 & seg_from == from | seg_to == to))   
  #     sub_clonal_seg <- df %>% filter(ratio != 1) %>% pull(seg_id)
  #     
  #     df <- df %>% filter(ratio == 1 & !(seg_id %in%  sub_clonal_seg) | ratio !=1) %>% 
  #       select(-from, -to, -seg_id) %>% 
  #       dplyr::rename(from = seg_from, to = seg_to) 
  #     
  #     return(df)
  #   } else{
  #     return(df)
  #   }
  #   }) %>% bind_rows()

  # Compute real ploidy
  ploidy = compute_true_ploidy(CNA_ProCESS)
  
  # Compute FGA
  tot_genome = CNA_ProCESS  %>% filter(!(chr %in% c('chrX', 'chrY')))  %>% select(chr, from, to) %>% distinct() %>%  mutate(len=to-from) %>% pull(len) %>% unique() %>% sum()
  altered = CNA_ProCESS %>% 
    filter(!(chr %in% c('chrX', 'chrY'))) %>% 
    mutate(len = to-from, CN = paste(major, minor, sep=':')) %>% 
    filter(ratio < 1 | CN !='1:1') %>% 
    select(-ratio, -CN, -major, -minor) %>% 
    distinct() %>% 
    pull(len) %>%
    unique() %>% 
    sum()
  fga = (altered/tot_genome)*100
  
  
  subclonal = CNA_ProCESS %>% 
    filter(!(chr %in% c('chrX', 'chrY'))) %>% 
    mutate(len = to-from, CN = paste(major, minor, sep=':')) %>% 
    filter(ratio < 1) %>% 
    select(-ratio, -CN, -major, -minor) %>% 
    distinct() %>% 
    pull(len) %>% 
    unique() %>% 
    sum()
  
  fgs = (subclonal/tot_genome)*100
  
  CNA_ProCESS  = CNA_ProCESS %>% mutate(type = ifelse(ratio != 1, 'subclonal', 'clonal'))
  
  #### ASCAT data
  message("Reading ASCAT data")
  
  ASCAT_output = read_ASCAT(spn_id,sample_id,coverage,purity)
  CNA_ascat = ASCAT_output[["CNA"]] 
  purity_ploidy_ascat = ASCAT_output[["purity_ploidy"]]
  
  #### Sequenza data
  message("Reading Sequenza data")
  
  Sequenza_output = read_Sequenza(spn_id,sample_id,coverage,purity)
  CNA_sequenza = Sequenza_output[["CNA"]] 
  purity_ploidy_sequenza = Sequenza_output[["purity_ploidy"]] 
  
  sequenza_centromere <- CNA_sequenza %>% 
    mutate(centromere = from - c(0, CNA_sequenza$to[1:nrow(CNA_sequenza)-1])) %>% 
    filter(centromere > 1e5) %>% select(chr, centromere)
  
  sequenza_df <- CNA_sequenza %>% 
    group_by(chr) %>% 
    summarize(start = min(from), end = max(to)) %>% 
    left_join(sequenza_centromere)
  sequenza_df <- absolute_to_relative_coordinates(sequenza_df, centromere = T)
  
  
  #### CNVkit data
  message("Reading CNVkit data")
  
  CNVkit_output = read_CNVkit(spn_id,sample_id,coverage,purity)
  CNA_cnvkit = CNVkit_output[["CNA"]]
  
  CNA_ProCESS_original <- CNA_ProCESS
  CNA_ProCESS <- CNA_ProCESS %>% filter(type == 'clonal')
  
  ############ Process data
  message("Create joint table ProCESS and ASCAT calls") 
  joint_segmentation_ascat = create_joint_segmentation(CNA_ProCESS, CNA_target=CNA_ascat, caller='ascat', chromosomes)
  joint_segmentation_ascat_long = joint_segmentation_ascat[['joint_segmentation_long']] %>% 
    filter(to - from > 1) %>% 
    left_join(ascat_df %>% dplyr::rename(chromosome = chr)) %>% 
    filter(from >= start) %>%  
    filter(to <= end) %>% 
    select(-start, -end)
  
  joint_segmentation_ascat[['joint_segmentation']] = joint_segmentation_ascat[['joint_segmentation']] %>% 
    filter(to - from > 1) %>% 
    left_join(ascat_df %>% dplyr::rename(chromosome = chr)) %>% 
    filter(from >= start) %>%  
    filter(to <= end) %>% 
    select(-start, -end)
  
  
  message("Create joint table ProCESS and Sequenza calls") 
  joint_segmentation_sequenza = create_joint_segmentation(CNA_ProCESS, CNA_target=CNA_sequenza, caller='sequenza', chromosomes)
  joint_segmentation_sequenza_long = joint_segmentation_sequenza[['joint_segmentation_long']] %>% 
    filter(to - from > 1) %>% 
    left_join(sequenza_df %>% dplyr::rename(chromosome = chr)) %>% 
    filter(from >= start) %>%  
    filter(to <= end) %>% 
    select(-start, -end) %>% 
    filter(!(to %in% sequenza_df$centromere)) %>% 
    filter(!is.na(INFERRED_CN))
  
  joint_segmentation_sequenza[['joint_segmentation']] = joint_segmentation_sequenza[['joint_segmentation']] %>% 
    filter(to - from > 1) %>% 
    left_join(sequenza_df %>% dplyr::rename(chromosome = chr)) %>% 
    filter(from >= start) %>%  
    filter(to <= end) %>% 
    select(-start, -end) %>% 
    filter(!(to %in% sequenza_df$centromere)) %>% 
    filter(!is.na(INFERRED_CN))
  
  
  message("Create joint table ProCESS and CNVkit calls") 
  joint_segmentation_cnvkit = create_joint_segmentation(CNA_ProCESS, CNA_target=CNA_cnvkit, caller='cnvkit', chromosomes)
  joint_segmentation_cnvkit_long = joint_segmentation_cnvkit[['joint_segmentation_long']] %>% 
    filter(to - from > 1) %>% 
    left_join(cnvkit_df %>% dplyr::rename(chromosome = chr)) %>% 
    filter(from >= start) %>%  
    filter(to <= end) %>% 
    select(-start, -end) %>% 
    filter(!(to %in% cnvkit_df$centromere))
  
  joint_segmentation_cnvkit[['joint_segmentation']] = joint_segmentation_cnvkit[['joint_segmentation']] %>% 
    filter(to - from > 1) %>% 
    left_join(cnvkit_df %>% dplyr::rename(chromosome = chr)) %>% 
    filter(from >= start) %>%  
    filter(to <= end) %>% 
    select(-start, -end) %>% 
    filter(!(to %in% cnvkit_df$centromere))
  
  
  # #### Calculate metrics
  # message("Compute metrics")
  # 
  # process <- CNA_ProCESS %>%  filter(ratio ==1) %>% mutate(cn = major+minor) %>% dplyr::rename(start = from, end = to, chrom = chr) %>%  select(chrom, start,end, cn)
  # ascat <- CNA_ascat %>%  mutate(cn = major+minor) %>% dplyr::rename(start = from, end = to, chrom = chr) %>%  select(chrom, start,end, cn)
  # metric_ascat <- calculate_metrics(gt = process, df_caller = ascat)
  # 
  # sequenza <- CNA_sequenza %>%  mutate(cn = major+minor) %>% dplyr::rename(start = from, end = to, chrom = chr) %>%  select(chrom, start,end, cn)
  # metric_sequenza <- calculate_metrics(gt = process, df_caller = sequenza)
  # 
  # cnvkit <- CNA_cnvkit %>% dplyr::rename(start = from, end = to, chrom = chr, cn = CN) %>%  select(chrom, start,end, cn)
  # metric_cnvkit <- calculate_metrics(gt = process, df_caller = cnvkit)
  # 
  ############ Compute CNA correctness
  message("Compute CNA correctness")
  
  ascat_correctness = compute_correctness(joint_segmentation_ascat_long, caller = 'ascat') 
  sequenza_correctness = compute_correctness(joint_segmentation_sequenza_long, caller = 'sequenza') 
  cnvkit_correctness =  compute_correctness(joint_segmentation_cnvkit_long, caller = 'cnvkit') 
  
  purity_correctness_ascat = purity_ploidy_ascat$AberrantCellFraction - as.double(purity)
  purity_correctness_sequenza = mean(purity_ploidy_sequenza$cellularity) - as.double(purity)
  
  ## Segmentation correctness
  # find out what proportion of each chromosome is not covered by the segmentation
  # breakpoints:
  #     - find matching BP and compute distance 
  #     - find missed/addes BP
  
  segmentation_ascat = segmentation_analysis(CNA_ProCESS, CNA_ascat, chromosomes, th=2e7)
  segmentation_sequenza = segmentation_analysis(CNA_ProCESS, CNA_sequenza, chromosomes, th=2e7)
  segmentation_cnvkit = segmentation_analysis(CNA_ProCESS, CNA_cnvkit, chromosomes, th=2e7)
  
  breakpoints_ascat = segmentation_ascat[["breakpoints"]]
  seg_summary_ascat = segmentation_ascat[["summary"]]
  breakpoints_sequenza = segmentation_sequenza[["breakpoints"]]
  seg_summary_sequenza = segmentation_sequenza[["summary"]]
  breakpoints_cnvkit = segmentation_cnvkit[["breakpoints"]]
  seg_summary_cnvkit = segmentation_cnvkit[["summary"]]
  
  ########### Plots
  message("Generate plots")
  
  CNA_ProCESS_rel <- absolute_to_relative_coordinates(CNA_ProCESS_original %>% dplyr::rename(start = from, end = to))%>% 
    mutate(type = ifelse(ratio != 1, 'sub-clonal', 'clonal')) %>% 
    mutate(c = ifelse(ratio < 1 & ratio > 0.5, '1', '2')) %>% 
    mutate(c = ifelse(ratio ==1, 1, c)) %>% 
    select(-ratio) %>% 
    pivot_longer(cols = c(major, minor)) %>% 
    mutate(name = ifelse(name == 'major', 'Major', 'minor')) %>% 
    mutate(name = paste0(name, c))  %>% 
    filter(chr %in% chromosomes)
  
  ProCESS_plt <-  CNAqc:::blank_genome(chromosomes = chromosomes) +
    geom_rect(data = CNA_ProCESS_rel, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill=type), alpha = .1) +
    geom_segment(data = CNA_ProCESS_rel %>% 
                   mutate(value = case_when(
                     name == 'Major1' ~ value + .04,
                     name == 'minor1' ~ value - .04,
                     name == 'Major2' ~ value + .10,
                     name == 'minor2' ~ value - .10,
                   )), aes(x=start, xend=end, y=value, col = name), size=1) +
    scale_color_manual('', values = color_process) + 
    scale_fill_manual('', values = color_type) + 
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.text=element_text(size=10)) +
    ggtitle(paste0('ProCESS simulation\n', spn_id, '-', sample_id),
            subtitle = element_text(paste0(
              'Coverage: ', coverage,
              '\nTrue purity: ',purity,
              '\nTrue ploidy: ',round(ploidy,2),
              '\nNumber of subclonal segments: ', CNA_ProCESS_original %>% filter(ratio != 1) %>% select(-ratio, -major, -minor) %>% unique() %>% nrow(),
              '\nFraction of genome altered: ', round(fga,2), '%',
              '\nFraction of genome subclonal: ', round(fgs,2), '%')))
  
    
  ascat_plt <- CNAqc:::blank_genome(chromosomes = chromosomes) + 
    geom_rect(data = joint_segmentation_ascat_long,
    aes(xmin=from, xmax=to, ymin=-Inf, ymax=Inf, fill=is_match)) + 
    geom_segment(data=joint_segmentation_ascat_long %>% filter(Type %in% c('INFERRED_Major1', 'INFERRED_minor1')) %>% 
                   mutate(Value = case_when(
                     Type == 'INFERRED_Major1' ~ Value + .06,
                     Type == 'INFERRED_minor1' ~ Value - .06,
                   )),
                 aes(x=from, xend=to, y=Value, color=Type), size=1)+
    scale_color_manual(values=color_by_state)+
    scale_fill_manual(values=fill_by_match) + 
    labs(fill = "", color = "")+
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.text=element_text(size=10)) +
    ggtitle('ASCAT',
              subtitle = element_text(paste0(
                'Proportion of genome inferred correctly - clonal: ', round(ascat_correctness$clonal,2)*100,'%',
                #'\nProportion of genome inferred correctly - all: ', round(ascat_correctness$all,2)*100,'%',
                '\nInferred purity: ', purity_ploidy_ascat$AberrantCellFraction,
                '\nInferred ploidy: ', round(purity_ploidy_ascat$Ploidy,2),
                '\nAverage breakpoint distance: ',round(seg_summary_ascat$av_distance))))
  
  
  
  sequenza_plt <- CNAqc:::blank_genome(chromosomes = chromosomes) + 
    geom_rect(data = joint_segmentation_sequenza_long,
              aes(xmin=from, xmax=to, ymin=-Inf, ymax=Inf, fill=is_match)) + 
    geom_segment(data=joint_segmentation_sequenza_long %>% filter(Type %in% c('INFERRED_Major1', 'INFERRED_minor1')) %>% 
                   mutate(Value = case_when(
                     Type == 'INFERRED_Major1' ~ Value + .06,
                     Type == 'INFERRED_minor1' ~ Value - .06,
                   )),
                 aes(x=from, xend=to, y=Value, color=Type), size=1)+
    scale_color_manual(values=color_by_state)+
    scale_fill_manual(values=fill_by_match) + 
    labs(fill = "", color = "")+
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.text=element_text(size=10))  +
    ggtitle('sequenza',
            subtitle = element_text(paste0(
              'Proportion of genome inferred correctly - clonal: ', round(sequenza_correctness$clonal,2)*100,'%',
              #'\nProportion of genome inferred correctly - all: ', round(sequenza_correctness$all,2)*100,'%',
              '\nInferred purity: ', purity_ploidy_sequenza$cellularity[[1]],
              '\nInferred ploidy: ', round(purity_ploidy_sequenza$ploidy.estimate[[1]],2),
              '\nAverage breakpoint distance: ',round(seg_summary_sequenza$av_distance))))
  
  
  
  cnvkit_plt = CNAqc:::blank_genome(chromosomes = chromosomes) + 
    geom_rect(data = joint_segmentation_cnvkit_long,
              aes(xmin=from, xmax=to, ymin=-Inf, ymax=Inf, fill=is_match)) +
    geom_segment(data=joint_segmentation_cnvkit_long %>% filter(Type == 'INFERRED_CN'), 
                 aes(x=from, xend=to, y=Value, color=Type), size=1)+
    scale_color_manual(values=color_by_state)+
    scale_fill_manual(values=fill_by_match)+
    labs(fill = "", color = "")+
    guides(fill = guide_legend(override.aes = list(alpha = 1))) +
    theme(legend.text=element_text(size=10))+
    ylab('Total CN') + 
    ggtitle('CNVkit',
            subtitle = element_text(paste0(
              'Proportion of genome inferred correctly - clonal: ', round(cnvkit_correctness$clonal,2)*100,'%',
              '\nAverage breakpoint distance: ',round(seg_summary_cnvkit$av_distance))))
  
  
  plt <- plt_snp + ProCESS_plt + ascat_plt + sequenza_plt + cnvkit_plt + patchwork::plot_layout(nrow = 6)
  
  ### Save reports 
  outdir <- paste0(data_dir,spn_id,"/validation/cna/",spn_id,"/",coverage,"x_",purity,'p/',sample_id,'/')
  dir.create(outdir, recursive = T, showWarnings = F)
  
  reportdir <- paste0(data_dir, spn_id, "/validation/cna/report/")
  dir.create(reportdir, recursive = T, showWarnings = F)
  
  filename <- paste(spn_id, coverage, purity, sample_id, sep='_')
  file_path <- file.path(reportdir, filename)
  
  ggsave(plt, file = paste0(outdir,'report.png'), height = 14, width = 8)
  ggsave(plt, file = paste0(file_path,'.png'), height = 14, width = 8)
  
  
  table_metric <- tibble('true_purity' = purity,
                         'true_ploidy' = ploidy,
                         "purity" = c(purity_ploidy_ascat$AberrantCellFraction,NA,purity_ploidy_sequenza$cellularity[[1]]),
                         "ploidy" = c(purity_ploidy_ascat$Ploidy,NA,purity_ploidy_sequenza$ploidy.estimate[[1]]),
                         "correctness_all" = c(ascat_correctness$all, cnvkit_correctness$all, sequenza_correctness$all),
                         "correctness_clonal" = c(ascat_correctness$clonal, cnvkit_correctness$clonal, sequenza_correctness$clonal),
                         "bp_distance" = c(seg_summary_ascat$av_distance,seg_summary_cnvkit$av_distance, seg_summary_sequenza$av_distance),
                         "tool" = c('ascat', 'cnvkit', 'sequenza'),
                         "sample" = sample_id,
                         "spn" = spn_id,
                         "coverage" = coverage,
                         "fga" = fga,
                         "fgs" = fgs
  )
    
  
  saveRDS(table_metric, file = paste0(outdir, 'metrics.rds'))
  saveRDS(list(ascat = ASCAT_output, ascat_long = joint_segmentation_ascat_long,
               sequenza = Sequenza_output, sequenza_long = joint_segmentation_sequenza_long,
               cnvkit = CNVkit_output, cnvkit_long = joint_segmentation_cnvkit_long), file = paste0(outdir, 'data.rds'))
  saveRDS(sequenza_df, file = paste0(outdir, 'sequenza_bp.rds'))
  
  message("Report saved for combination: purity=", purity, ", cov=", coverage, ', sample=', sample_id)
}

