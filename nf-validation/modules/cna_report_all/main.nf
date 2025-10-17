process CNA_REPORT_ALL {
    tag "${meta.spn}-generate_report"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v2' :
        'docker.io/lvaleriani/process_validation:v2' }"

    input:
    tuple val(meta), val(rds_list), val(samples), val(combinations)

    output:
    tuple val(meta), path("*rds"),	emit:rds
    tuple val(meta), path("*pdf"),	emit:plot
    
    publishDir "${params.outdir}/${meta.spn}/cna/all", mode: 'copy'


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

    data_dir = '/orfeo/scratch/cdslab/shared/SCOUT/'
    spn_id = "${meta.spn}"

    combinations = substr("$combinations", 2, nchar("$combinations")-1)
    combinations = strsplit(combinations, ", ")[[1]]

    samples = strsplit("$samples", '], ') %>% 
            unlist() %>% 
            stringr::str_replace_all('\\\\[', '') %>% 
            stringr::str_replace_all('\\\\]', '') %>% 
            strsplit(', ') 
    names(samples) = combinations

    rds = strsplit("$rds_list", '], ') %>% 
            unlist() %>% 
            stringr::str_replace_all('\\\\[', '') %>% 
            stringr::str_replace_all('\\\\]', '') %>% 
            strsplit(', ') 
    names(rds) = combinations

    for (c in combinations){
        names(rds[[c]]) = samples[[c]]
    }

    df_metric = lapply(combinations, function(comb) {
        tmp = strsplit(as.character(comb), '_') %>% unlist()
        coverage = tmp[[1]]
        purity = tmp[[2]]
        
        data = rds[[comb]]
        
        tmp_df = lapply(names(data), FUN = function(sample){
            metrics = readRDS(data[[sample]]) %>% mutate(true_purity = as.numeric(true_purity), coverage = as.numeric(coverage))
        }) %>% bind_rows()
    }) %>% bind_rows()
    saveRDS(df_metric, paste0(spn_id,'_final_df.rds'))
    
    
    plt_data <- df_metric %>% 
      select(sample, fga,fgs) %>% 
      distinct() %>% 
      pivot_longer(cols = c(fga,fgs)) %>% 
      ggplot() + 
      geom_col(aes(x = sample, y = value,fill=name ),position = position_dodge()) +
      scale_fill_manual('', values = c('indianred', 'orange'))  +
      ylab('% of genome') +
      theme_minimal() 
    
    plt_ploidy <- df_metric %>% 
      select(sample, true_ploidy) %>% 
      distinct() %>% 
      mutate(true_ploidy = round(true_ploidy,1)) %>% 
      ggplot(aes(x = sample, y = true_ploidy )) + 
      geom_col() +
      theme_minimal()
    
    #plt_bp <- df_bp %>% 
    #filter(chr == 'genome') %>% 
    #pivot_longer(cols = c(precision, recall, f1)) %>% 
    #ggplot() +
    #geom_boxplot(aes(x = as.factor(true_purity), y = value, col = tool)) +
    #scale_color_manual(values = color_caller) +
    #facet_grid(name~coverage)  +
    #xlab('purity') + 
    #theme_bw() 
    
    plt <- plt_ploidy + 
      plt_data + 
      
      df_metric %>%
      mutate(delta_purity = as.numeric(true_purity) - as.numeric(purity)) %>%
      ggplot(aes(x = sample, y = delta_purity, fill = tool)) +
      geom_col(position = position_dodge(width = 0.35), width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = color_caller) +
      ylab("true_purity - inferred_purity") +
      facet_grid(as.numeric(coverage) ~ as.numeric(true_purity)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      
      df_metric %>%
      mutate(delta_ploidy = as.numeric(true_ploidy) - as.numeric(ploidy)) %>%
      ggplot(aes(x = sample, y = delta_ploidy, fill = tool)) +
      geom_col(position = position_dodge(width = 0.35), width = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = color_caller) +
      ylab("true_ploidy - inferred_ploidy") +
      facet_grid(as.numeric(coverage) ~ as.numeric(true_purity)) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      
      df_metric %>% 
      ggplot() + 
      geom_col(aes(x = sample, y = correctness_clonal*100, fill = tool), position = position_dodge()) + 
      scale_fill_manual(values = color_caller) +
      ylab('% CN correctness') +
      theme_bw() +
      facet_grid(as.numeric(coverage)  ~ as.numeric(true_purity))  +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      #plt_bp + 
    
    plot_layout(design = 'AB\nCC\nDD\nEE\nFF') +
    plot_annotation(title = spn_id) + plot_layout(guides = 'collect')
    
    ggsave(filename = paste0(spn_id, '.pdf'), 
        plot = plt, 
        dpi = 400, 
        height = 14, width = 12, units = 'in')

    """
}