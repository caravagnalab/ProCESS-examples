process CNA_REPORT_COMBINATION {
    tag "${meta.spn}-generate_report"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v1' :
        'docker.io/lvaleriani/process_validation:v1' }"

    input:
    tuple val(meta), val(path)

    output:
    tuple val(meta), path("metrics.rds"),	    emit:rds_metric
    tuple val(meta), path("data.rds"),	        emit:rds_data
    tuple val(meta), path("metrics_bp.rds"),	emit:rds_metric_bp

    tuple val(meta), path("*png"),	emit:report
    
    publishDir "${params.outdir}/${meta.spn}/cna/${meta.coverage}x_${meta.purity}/${meta.sample}", mode: 'copy'

    script:

    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(ProCESS)
    library(optparse)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
    
    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/getters/sarek_getters.R")
    source("${projectDir}/bin/cna/utils.R")

    ascat_df_rel = readRDS('${projectDir}/bin/data/ascat.rds') %>% dplyr::rename(chr=V1)
    ascat_df_abs = absolute_to_relative_coordinates(ascat_df_rel)

    cnvkit_df_rel = readRDS('${projectDir}/bin/data/cnvkit.rds')
    cnvkit_df_abs = absolute_to_relative_coordinates(cnvkit_df_rel, centromere = T)     

    data_dir = "${path}"
    spn_id = "${meta.spn}"
    coverage = "${meta.coverage}"
    purity = "${meta.purity}"
    sample_id = "${meta.sample}"

    chromosomes = c(paste0('chr',1:22))
    process_rds = readRDS(get_mutations(spn = spn_id, coverage = coverage, purity = purity, type = 'tumour'))
    process_rds_long = ProCESS::seq_to_long(process_rds)
    rm(process_rds)

    process_normal = readRDS(get_mutations(spn = spn_id, coverage = coverage, purity = purity, type = 'normal'))
    process_normal_long = ProCESS::seq_to_long(process_normal)
    rm(process_normal)
    process_normal_long = process_normal_long %>% ungroup() %>% select(chr, from, to, NV, DP, VAF)

    process_rds_long = process_rds_long %>% 
        ungroup() %>% 
        filter(classes == 'germinal') %>% 
        select(chr, from, to, ref, alt, NV, DP, VAF, sample_name)

    cov_T = as.numeric(coverage)
    cov_N = 30
    ratio = cov_N/cov_T

    process_rds_long = process_rds_long %>% 
    left_join(process_normal_long, by = join_by(chr, from, to), suffix = c('_T', '_N')) %>% 
    filter(ref %in% c('A', 'T', 'C', 'G')) %>% 
    filter(alt %in% c('A', 'T', 'C', 'G')) %>%
    select(chr, from, to, NV_T, DP_T, VAF_T, sample_name, DP_N, VAF_N) %>% 
    mutate(DR = (DP_T/DP_N)*ratio) %>% 
    filter(VAF_N > 0.3, VAF_N < 0.7) 
    
    process_rds = process_rds_long %>% filter(sample_name == sample_id) %>% filter(DP_T > 10) %>% sample_n(1e4)
    process_rds = absolute_to_relative_coordinates_muts(process_rds %>% mutate(chr = paste0('chr',chr)))
    
    ProCESS_output = read_ProCESS(spn_id,sample_id,coverage,purity)
    CNA_ProCESS = ProCESS_output[["CNA"]] 
    
    CNA_ProCESS = CNA_ProCESS %>% 
        filter(ratio > 0.09) %>% 
        group_by(seg_id) %>% 
        mutate(ratio = ratio/sum(ratio)) %>% 
        #mutate(ratio = round(ratio,1)) %>% 
        filter(ratio != 0)  %>% 
        ungroup()
    
    final_CNA_ProCESS = lapply(chromosomes, FUN = function(c){
        df = CNA_ProCESS %>% 
        filter(chr == c) 
        
        segments = unique(df\$seg_id)
        new_df = df %>% filter(seg_id == segments[1]) %>% arrange(ratio)
        
        if (length(segments)>=2){
        for (s in 2:length(segments)){
            curr_seg = segments[s]
            curr_df = df %>% filter(seg_id == curr_seg)  %>% arrange(ratio)
        
            prev_seg = segments[s-1]
            prev_df = new_df %>% filter(seg_id == prev_seg) %>% arrange(ratio)
            
            if (identical(prev_df %>% select(major, minor, ratio), curr_df %>% select(major, minor, ratio))){
            tmp_df = curr_df %>% mutate(from = prev_df\$from)
            new_df = new_df  %>% filter(seg_id != prev_seg)
            new_df = new_df %>% bind_rows(tmp_df)
            } else {
            new_df = new_df %>% bind_rows(curr_df)
            }
        }
        }
        return(new_df)
    }) %>% bind_rows() %>% select(-seg_id)
    
    # Compute real ploidy
    ploidy = compute_true_ploidy(CNA_ProCESS)
    
    # Compute FGA 
    tot_genome = final_CNA_ProCESS  %>% 
        filter(!(chr %in% c('chrX', 'chrY')))  %>% 
        select(chr, from, to) %>% 
        distinct() %>%  
        mutate(len=to-from) %>% 
        pull(len) %>% 
        unique() %>% 
        sum()

    altered = final_CNA_ProCESS %>% 
        filter(!(chr %in% c('chrX', 'chrY'))) %>% 
        mutate(len = to-from, CN = paste(major, minor, sep=':')) %>% 
        filter(ratio < 1 | CN !='1:1') %>% 
        select(-ratio, -CN, -major, -minor) %>% 
        distinct() %>% 
        pull(len) %>%
        unique() %>% 
        sum()
    fga = (altered/tot_genome)*100
    
    subclonal = final_CNA_ProCESS %>% 
        filter(!(chr %in% c('chrX', 'chrY'))) %>% 
        mutate(len = to-from, CN = paste(major, minor, sep=':')) %>% 
        filter(ratio < 1) %>% 
        select(-ratio, -CN, -major, -minor) %>% 
        distinct() %>% 
        pull(len) %>% 
        unique() %>% 
        sum()
    
    fgs = (subclonal/tot_genome)*100
    CNA_ProCESS  = final_CNA_ProCESS %>% mutate(type = ifelse(ratio != 1, 'subclonal', 'clonal'))
    
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
    
    sequenza_centromere = CNA_sequenza %>% 
        mutate(centromere = from - c(0, CNA_sequenza\$to[1:nrow(CNA_sequenza)-1])) %>% 
        filter(centromere > 1e5) %>% select(chr, centromere)
    
    sequenza_df = CNA_sequenza %>% 
        group_by(chr) %>% 
        summarize(start = min(from), end = max(to)) %>% 
        left_join(sequenza_centromere)
    sequenza_df = absolute_to_relative_coordinates(sequenza_df, centromere = T)
    
    
    #### CNVkit data
    message("Reading CNVkit data")
    
    CNVkit_output = read_CNVkit(spn_id,sample_id,coverage,purity)
    CNA_cnvkit = CNVkit_output[["CNA"]]
    
    
    #### Battenberg data
    Battenberg_output = read_Battenberg(spn_id,sample_id,coverage,purity)
    CNA_battenberg = Battenberg_output[["CNA"]] 
    purity_ploidy_battenberg = Battenberg_output[["purity_ploidy"]] 
    
    battenberg_centromere = CNA_battenberg %>% 
      mutate(centromere = from - c(0, CNA_battenberg\$to[1:nrow(CNA_battenberg)-1])) %>% 
      filter(centromere > 1e5) %>% select(chr, centromere)
    
    battenberg_df = CNA_battenberg %>% 
      group_by(chr) %>% 
      summarize(start = min(from), end = max(to)) %>% 
      left_join(battenberg_centromere)
    battenberg_df = absolute_to_relative_coordinates(battenberg_df, centromere = T)

    
    CNA_ProCESS_original = CNA_ProCESS
    
    ############ Process data
    message("Create joint table ProCESS and ASCAT calls") 
    joint_segmentation_ascat = create_joint_segmentation(final_CNA_ProCESS, CNA_target=CNA_ascat, caller='ascat', chromosomes)
    joint_segmentation_ascat_long = joint_segmentation_ascat[['joint_segmentation_long']] %>% 
        filter(to - from > 1) %>% 
        left_join(ascat_df_abs %>% dplyr::rename(chromosome = chr)) %>% 
        filter(from >= start) %>%  
        filter(to <= end) %>% 
        select(-start, -end)
    
    joint_segmentation_ascat[['joint_segmentation']] = joint_segmentation_ascat[['joint_segmentation']] %>% 
        filter(to - from > 1) %>% 
        left_join(ascat_df_abs %>% dplyr::rename(chromosome = chr)) %>% 
        filter(from >= start) %>%  
        filter(to <= end) %>% 
        select(-start, -end)
    
    
    message("Create joint table ProCESS and Sequenza calls") 
    joint_segmentation_sequenza = create_joint_segmentation(final_CNA_ProCESS, CNA_target=CNA_sequenza, caller='sequenza', chromosomes)
    joint_segmentation_sequenza_long = joint_segmentation_sequenza[['joint_segmentation_long']] %>% 
        filter(to - from > 1) %>% 
        left_join(sequenza_df %>% dplyr::rename(chromosome = chr)) %>% 
        filter(from >= start) %>%  
        filter(to <= end) %>% 
        select(-start, -end) %>% 
        filter(!(to %in% sequenza_df\$centromere)) %>% 
        filter(!is.na(INFERRED_CN))
    
    joint_segmentation_sequenza[['joint_segmentation']] = joint_segmentation_sequenza[['joint_segmentation']] %>% 
        filter(to - from > 1) %>% 
        left_join(sequenza_df %>% dplyr::rename(chromosome = chr)) %>% 
        filter(from >= start) %>%  
        filter(to <= end) %>% 
        select(-start, -end) %>% 
        filter(!(to %in% sequenza_df\$centromere)) %>% 
        filter(!is.na(INFERRED_CN))
    
    
    message("Create joint table ProCESS and CNVkit calls") 
    joint_segmentation_cnvkit = create_joint_segmentation(final_CNA_ProCESS, CNA_target=CNA_cnvkit, caller='cnvkit', chromosomes)
    joint_segmentation_cnvkit_long = joint_segmentation_cnvkit[['joint_segmentation_long']] %>% 
        filter(to - from > 1) %>% 
        left_join(cnvkit_df_abs %>% dplyr::rename(chromosome = chr)) %>% 
        filter(from >= start) %>%  
        filter(to <= end) %>% 
        select(-start, -end) %>% 
        filter(!(to %in% cnvkit_df_abs\$centromere))
    
    joint_segmentation_cnvkit[['joint_segmentation']] = joint_segmentation_cnvkit[['joint_segmentation']] %>% 
        filter(to - from > 1) %>% 
        left_join(cnvkit_df_abs %>% dplyr::rename(chromosome = chr)) %>% 
        filter(from >= start) %>%  
        filter(to <= end) %>% 
        select(-start, -end) %>% 
        filter(!(to %in% cnvkit_df_abs\$centromere))
        
    message("Create joint table ProCESS and Battenberg calls") 
    joint_segmentation_batt = create_joint_segmentation(CNA_ProCESS = final_CNA_ProCESS, CNA_target=CNA_battenberg, caller='battenberg', chromosomes = chromosomes)
    joint_segmentation_batt_long = joint_segmentation_batt[['joint_segmentation_long']] %>% 
      filter(to - from > 1) %>% 
      left_join(battenberg_df %>% dplyr::rename(chromosome = chr)) %>% 
      filter(from >= start) %>%  
      filter(to <= end) %>% 
      select(-start, -end) %>% 
      filter(!(to %in% battenberg_df\$centromere)) %>% 
      filter(!is.na(INFERRED_CN))
    
    joint_segmentation_batt[['joint_segmentation']] = joint_segmentation_batt[['joint_segmentation']] %>% 
      filter(to - from > 1) %>% 
      left_join(battenberg_df %>% dplyr::rename(chromosome = chr)) %>% 
      filter(from >= start) %>%  
      filter(to <= end) %>% 
      select(-start, -end) %>% 
      filter(!(to %in% battenberg_df\$centromere)) %>% 
      filter(!is.na(INFERRED_CN))
    

    # compute bp metrics
    # Split segemnts that span centromeres
    message("Compute BP correctness")
    all_bp_metrics = tibble()
    for (tool in c("ascat", "sequenza", "cnvkit", "battenberg")) {
        print(tool)
        
        if (tool == "ascat") {
        callable_regions_df = ascat_df_rel
        CNA_pred = lapply(chromosomes, function(c) {
            filter_genomic_signal(
            signal_df = CNA_ascat %>% dplyr::filter(chr == c) %>% dplyr::mutate(ratio = 1, type = 1), 
            centromere_df = convert_centromeres_to_relative(CNAqc::chr_coordinates_GRCh38) %>% dplyr::filter(chr == c), 
            callable_regions_df = callable_regions_df
            )
        }) %>% do.call("bind_rows", .)  
        } else if (tool == "sequenza") {
        callable_regions_df = relative_to_absolute_coordinates(sequenza_df %>% select(-centromere))
        CNA_pred = lapply(chromosomes, function(c) {
            df = filter_genomic_signal(
            signal_df = CNA_sequenza %>% dplyr::filter(chr == c) %>% dplyr::mutate(ratio = 1, type = 1), 
            centromere_df = convert_centromeres_to_relative(CNAqc::chr_coordinates_GRCh38) %>% dplyr::filter(chr == c), 
            callable_regions_df = callable_regions_df
            )
            if (nrow(df) > 1){
            return(df)
            }
        }) %>% do.call("bind_rows", .)  
        } else if (tool == "cnvkit") {
        callable_regions_df = cnvkit_df_rel %>% select(-centromere)
        CNA_pred = lapply(chromosomes, function(c) {
            df = filter_genomic_signal(
            signal_df = CNA_cnvkit %>% dplyr::filter(chr == c) %>% dplyr::mutate(ratio = 1, type = 1, minor = "0") %>% dplyr::rename(major = CN), 
            centromere_df = convert_centromeres_to_relative(CNAqc::chr_coordinates_GRCh38) %>% dplyr::filter(chr == c), 
            callable_regions_df = callable_regions_df %>% dplyr::filter(chr == c)
            )
            if (nrow(df) > 1){
            return(df)
            }
        }) %>% do.call("bind_rows", .)  
        } else if (tool == 'battenberg'){
          callable_regions_df = relative_to_absolute_coordinates(battenberg_df %>% select(-centromere))
          CNA_pred = lapply(chromosomes, function(c) {
            df = filter_genomic_signal(
              signal_df = CNA_battenberg %>% dplyr::filter(chr == c) %>% dplyr::rename(major = Major1, minor = minor1) %>% dplyr::mutate(ratio = 1, type = 1) %>% select(-Major2, -minor2, -frac2_A, -frac1_A), 
              centromere_df = convert_centromeres_to_relative(CNAqc::chr_coordinates_GRCh38) %>% dplyr::filter(chr == c), 
              callable_regions_df = callable_regions_df
            )
            if (nrow(df) > 1){
              return(df)
            }
          }) %>% do.call("bind_rows", .) 
        } else {
          stop()
        }
        
        CNA_ProCESS_curr = lapply(chromosomes, function(c) {
        remove_centromeric_signal_gr(
            signal_df = CNA_ProCESS %>% dplyr::filter(chr == c), 
            centromere_df = convert_centromeres_to_relative(CNAqc::chr_coordinates_GRCh38) %>% dplyr::filter(chr == c), 
            callable_regions_df = callable_regions_df
        )
        }) %>% do.call("bind_rows", .)  
        
        
        CNA_confusion_df = lapply(chromosomes, function(c) {
        
        gt_bp = unique(
            c(
            CNA_ProCESS_curr %>% dplyr::filter(chr == c) %>% dplyr::pull(from),
            CNA_ProCESS_curr %>% dplyr::filter(chr == c) %>% dplyr::pull(to) - 1
            )
        )
        
        pred_bp = unique(
            c(
            CNA_pred %>% dplyr::filter(chr == c) %>% dplyr::pull(from),
            CNA_pred %>% dplyr::filter(chr == c) %>% dplyr::pull(to) - 1
            )
        )
        
        compute_confusion_matrix(pred_bp, truth = gt_bp, delta = 1e6, delta_merge = 0) %>% dplyr::mutate(chr = c)
        }) %>% do.call("bind_rows", .)
        
        
        CNA_metrics_df = compute_CNA_metrics(CNA_confusion_df) %>% mutate(tool = tool)
        all_bp_metrics = bind_rows(all_bp_metrics,CNA_metrics_df)
    }
    
    metric_plot = all_bp_metrics %>% 
        filter(chr == 'genome') %>% 
        pivot_longer(cols = c(precision, recall, f1))  %>% 
        ggplot() + 
        geom_col(aes(x = name, y = value, fill = tool), position=position_dodge()) + 
        scale_fill_manual(values = color_caller) + 
        xlab('metric') + 
        theme_bw() + theme(legend.text=element_text(size=10), legend.position = 'bottom')
    
    ############ Compute CNA correctness
    message("Compute CNA correctness")
    
    ascat_correctness = compute_correctness(df = joint_segmentation_ascat_long, caller = 'ascat') 
    sequenza_correctness = compute_correctness(joint_segmentation_sequenza_long, caller = 'sequenza') 
    cnvkit_correctness =  compute_correctness(joint_segmentation_cnvkit_long, caller = 'cnvkit') 
    battenberg_correctness = compute_correctness(df = joint_segmentation_batt_long, caller = 'battenberg') 

    ########### Plots
    message("Generate plots")
    
    CNA_ProCESS_rel = absolute_to_relative_coordinates(final_CNA_ProCESS %>% dplyr::rename(start = from, end = to))%>% 
        mutate(type = ifelse(ratio != 1, 'sub-clonal', 'clonal')) %>% 
        mutate(c = ifelse(ratio < 1 & ratio > 0.5, '1', '2')) %>% 
        mutate(c = ifelse(ratio == 1, 1, c)) %>% 
        pivot_longer(cols = c(major, minor)) %>% 
        mutate(name = ifelse(name == 'major', 'Major', 'minor')) %>% 
        mutate(name = paste0(name, c))  %>% 
        filter(chr %in% chromosomes)
    
    plt_snp = CNAqc:::blank_genome(chromosomes = chromosomes) + 
        geom_point(data = process_rds, aes(x = from, y = VAF_T), size = .1) + 
        ylab('BAF') + 
        theme(legend.text=element_text(size=10)) + 
        ggtitle(paste0('ProCESS simulation\n', spn_id, '-', sample_id),
                subtitle = element_text(paste0(
                'Coverage: ', coverage,
                '\nTrue purity: ',purity,
                '\nTrue ploidy: ',round(ploidy,2),
                '\nNumber of clonal segments: ', CNA_ProCESS %>% filter(ratio == 1) %>% select(-ratio, -major, -minor) %>% unique() %>% nrow(),
                '\nNumber of subclonal segments: ', CNA_ProCESS %>% filter(ratio != 1) %>% select(-ratio, -major, -minor) %>% unique() %>% nrow(),
                '\nFraction of genome altered: ', round(fga,2), '%',
                '\nFraction of genome subclonal: ', round(fgs,2), '%'))) + 
        CNAqc:::blank_genome(chromosomes = chromosomes) + 
        geom_point(data = process_rds, aes(x = from, y = DR), size = .1) + 
        ylab('DR') +
        ylim(-1,3) + 
        plot_layout(nrow=2) +
        theme(legend.text=element_text(size=10))
    
    ProCESS_annotations <- CNA_ProCESS_rel %>%
      filter(type == "sub-clonal") %>% 
      mutate(
        ratio_label = paste0(round(ratio, 1)),
        xpos = (start + end) / 2,
        ypos = .5  # Adjust as needed based on your y-axis scale
      ) %>% 
      group_by(chr, start, end, xpos, ypos)
    
    if (nrow(ProCESS_annotations) > 1){
      ProCESS_annotations = ProCESS_annotations %>% summarise(ratio_label = max(ratio_label))
    }
    
    ProCESS_plt =  CNAqc:::blank_genome(chromosomes = chromosomes) +
        geom_rect(data = CNA_ProCESS_rel, aes(xmin=start, xmax=end, ymin=-Inf, ymax=Inf, fill=type), alpha = .1) +
        geom_segment(data = CNA_ProCESS_rel %>% 
                    mutate(value = case_when(
                        name == 'Major1' ~ value + .04,
                        name == 'minor1' ~ value - .04,
                        name == 'Major2' ~ value + .10,
                        name == 'minor2' ~ value - .10,
                    )), aes(x=start, xend=end, y=value, col = name), size=1) +
        geom_text(data = ProCESS_annotations, aes(x = xpos, y = ypos, label = ratio_label), size = 3) +
        scale_color_manual('', values = color_process) + 
        scale_fill_manual('', values = color_type) + 
        guides(fill = guide_legend(override.aes = list(alpha = 1))) +
        theme(legend.text=element_text(size=10))
    
    joint_segmentation_ascat_long = joint_segmentation_ascat_long %>% mutate(is_match = ifelse(type == 'subclonal', 'subclonal', is_match))
    ascat_plt = CNAqc:::blank_genome(chromosomes = chromosomes) + 
        geom_rect(data = joint_segmentation_ascat_long %>% select(from, to, chromosome, is_match) %>% distinct(),
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
                    'Percentage of clonal genome inferred correctly: ', round(ascat_correctness\$clonal,2)*100,'%',
                    '\nInferred purity: ', purity_ploidy_ascat\$AberrantCellFraction,
                    '\nInferred ploidy: ', round(purity_ploidy_ascat\$Ploidy,2))))
    
    
    joint_segmentation_sequenza_long = joint_segmentation_sequenza_long %>% mutate(is_match = ifelse(type == 'subclonal', 'subclonal', is_match))
    sequenza_plt = CNAqc:::blank_genome(chromosomes = chromosomes) + 
        geom_rect(data = joint_segmentation_sequenza_long %>% select(from, to, chromosome, is_match) %>% distinct(),
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
                'Percentage of clonal genome inferred correctly: ', round(sequenza_correctness\$clonal,2)*100,'%',
                '\nInferred purity: ', purity_ploidy_sequenza\$cellularity[[1]],
                '\nInferred ploidy: ', round(purity_ploidy_sequenza\$ploidy.estimate[[1]],2))))
    
    
    joint_segmentation_cnvkit_long = joint_segmentation_cnvkit_long %>% mutate(is_match = ifelse(type == 'subclonal', 'subclonal', is_match))
    cnvkit_plt = CNAqc:::blank_genome(chromosomes = chromosomes) + 
        geom_rect(data = joint_segmentation_cnvkit_long %>% select(from, to, chromosome, is_match) %>% distinct(),
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
                'Percentage of clonal genome inferred correctly: ', round(cnvkit_correctness\$clonal,2)*100,'%')))
    
    batt_annotations <- joint_segmentation_batt_long %>%
        filter(type == "subclonal") %>%
        filter(!is.na(INFERRED_ratio)) %>%
        filter(INFERRED_ratio < 1) %>% 
        mutate(
          ratio_label = round(INFERRED_ratio, 1),
          xpos = (from + to) / 2,
          ypos = .5
        ) %>% select(chromosome, from, to,ratio_label,xpos, ypos) %>% distinct()
    
    plt_battenberg <- CNAqc:::blank_genome(chromosomes = chromosomes) +
        geom_rect(data = joint_segmentation_batt_long %>% select(from, to, chromosome, is_match) %>% distinct(), 
        aes(xmin=from, xmax=to, ymin=-Inf, ymax=Inf, fill=is_match)) +
        geom_segment(data = joint_segmentation_batt_long %>% 
                      filter(Type %in% c('INFERRED_Major1', 'INFERRED_minor1', 'INFERRED_Major2','INFERRED_minor2')) %>% 
                       mutate(Value = case_when(
                         Type == 'INFERRED_Major1' ~ Value + .04,
                         Type == 'INFERRED_minor1' ~ Value - .04,
                         Type == 'INFERRED_Major2' ~ Value + .10,
                         Type == 'INFERRED_minor2' ~ Value - .10,
                       )), aes(x=from, xend=to, y=Value, col = Type), size=1) +
        geom_text(data = batt_annotations, aes(x = xpos, y = ypos, label = ratio_label), size = 3) +
        scale_color_manual('', values = color_by_state) + 
        scale_fill_manual('', values = fill_by_match) + 
        guides(fill = guide_legend(override.aes = list(alpha = 1))) +
        theme(legend.text=element_text(size=10)) +  
        ggtitle('Battenberg',
               subtitle = element_text(paste0(
                 'Percentage of clonal genome inferred correctly: ', round(battenberg_correctness\$clonal,2)*100,'%',
                 '\nPercentage of subclonal genome inferred correctly: ', round(battenberg_correctness\$subclonal,2)*100,'%',
                 '\nInferred purity: ', round(purity_ploidy_battenberg\$purity,2),
                 '\nInferred ploidy: ', round(purity_ploidy_battenberg\$ploidy,2))))
    
    
    plt = plt_snp + ProCESS_plt + ascat_plt + sequenza_plt + cnvkit_plt + plt_battenberg + metric_plot + patchwork::plot_layout(nrow = 8)
    
    ### Save reports     
    ggsave(plt, file = paste0('report.png'), height = 16, width = 8)
    table_metric = tibble('true_purity' = purity,
                            'true_ploidy' = ploidy,
                            "purity" = c(purity_ploidy_ascat\$AberrantCellFraction,NA,purity_ploidy_sequenza\$cellularity[[1]], purity_ploidy_battenberg\$purity),
                            "ploidy" = c(purity_ploidy_ascat\$Ploidy,NA,purity_ploidy_sequenza\$ploidy.estimate[[1]], purity_ploidy_battenberg\$ploidy),
                            "correctness_clonal" = c(ascat_correctness\$clonal, cnvkit_correctness\$clonal, sequenza_correctness\$clonal, battenberg_correctness\$clonal),
                            "correctness_subclonal" = c(NA, NA, NA, battenberg_correctness\$subclonal),
                            "tool" = c('ascat', 'cnvkit', 'sequenza', 'battenberg'),
                            "sample" = sample_id,
                            "spn" = spn_id,
                            "coverage" = coverage,
                            "fga" = fga,
                            "fgs" = fgs
    )
    
    saveRDS(table_metric, file = paste0('metrics.rds'))
    saveRDS(all_bp_metrics %>% mutate("fga" = fga, 
                                        "fgs" = fgs,
                                        "sample" = sample_id,
                                        "spn" = spn_id,
                                        "coverage" = coverage,
                                        'true_purity' = purity), file = paste0('metrics_bp.rds'))
    saveRDS(list(ascat = ASCAT_output, ascat_long = joint_segmentation_ascat_long,
                sequenza = Sequenza_output, sequenza_long = joint_segmentation_sequenza_long,
                cnvkit = CNVkit_output, cnvkit_long = joint_segmentation_cnvkit_long,
                batt = Battenberg_output, Batt = joint_segmentation_batt_long), file = paste0('data.rds'))
    
    message("Report saved for combination: purity=", purity, ", cov=", coverage, ', sample=', sample_id)
    """
}