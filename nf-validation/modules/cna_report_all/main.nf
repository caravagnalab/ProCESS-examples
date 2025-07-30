process CNA_REPORT_ALL {
    tag "${meta.spn}-generate_report"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v1' :
        'docker.io/lvaleriani/process_validation:v1' }"

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

    plt <- df_metric %>% 
        ggplot() + 
        geom_point(aes(x = sample, y = as.numeric(true_purity) - as.numeric(purity), col = tool)) +
        scale_color_manual(values = color_caller) +
        geom_hline(aes(yintercept = 0)) +
        ylab('true_purity - inferred_purity') +
        theme_bw() +
        facet_grid(as.numeric(coverage) ~ as.numeric(true_purity)) +

        df_metric %>% 
        ggplot() + 
        geom_point(aes(x = sample, y = as.numeric(true_ploidy) - as.numeric(ploidy), col = tool)) + 
        scale_color_manual(values = color_caller) +
        geom_hline(aes(yintercept = 0)) +
        ylab('true_ploidy - inferred_ploidy') +
        theme_bw() +
        facet_grid(as.numeric(coverage)  ~ as.numeric(true_purity))  +
        
        df_metric %>% 
        ggplot() + 
        geom_point(aes(x = sample, y = correctness_clonal, col = tool)) + 
        scale_color_manual(values = color_caller) +
        ylab('% correctness') +
        theme_bw() +
        facet_grid(as.numeric(coverage)  ~ as.numeric(true_purity))  +
        
        plot_layout(nrow = 3) +
        plot_annotation(title = spn_id) + plot_layout(guides = 'collect')

    ggsave(filename = file.path(paste0(spn_id, '.pdf')), plot = plt, dpi = 400, height = 12, width = 12, units = 'in')


    p1 = df_metric %>%
        filter(tool != 'cnvkit') %>%
        group_by(sample) %>%
        dplyr::mutate(sample_fga = mean(fga, na.rm = TRUE)) %>%
        ungroup() %>%
        dplyr::mutate(
            purity_err =  true_purity - purity,
            sample_fga_rank = rank(sample_fga)
        ) %>%
        ggplot(mapping = aes(x = true_purity, y = purity_err, col = sample_fga_rank, group = sample, shape = sample)) +
        geom_point(size = 3) +
        geom_line() +
        facet_grid(coverage~tool) +
        scale_fill_viridis_c(name = "FGA Rank") +
        scale_color_viridis_c(name = "FGA Rank") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(x = "True Purity", y = "Delta Purity (True - Predicted)", fill = "FGA", title = 'FGA')

    p2 = df_metric %>%
        filter(tool != 'cnvkit') %>%
        group_by(sample) %>%
        dplyr::mutate(sample_fgs = mean(fgs, na.rm = TRUE)) %>%
        ungroup() %>%
        dplyr::mutate(
            purity_err =  true_purity - purity,
            sample_fgs_rank = rank(sample_fgs)
        ) %>%
        ggplot(mapping = aes(x = true_purity, y = purity_err, col = sample_fgs_rank, group = sample, shape = sample)) +
        geom_point(size = 3) +
        geom_line() +
        facet_grid(coverage~tool) +
        scale_fill_viridis_c(name = "FGS Rank") +
        scale_color_viridis_c(name = "FGS Rank") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(x = "True Purity", y = "Delta Purity (True - Predicted)", fill = "FGS", title ='FGS')

    p3 = df_metric %>%
        filter(tool != 'cnvkit') %>%
        group_by(sample) %>%
        dplyr::mutate(sample_fga = mean(fga, na.rm = TRUE)) %>%
        ungroup() %>%
        dplyr::mutate(
            ploidy_err =  true_ploidy - ploidy,
            sample_fga_rank = rank(sample_fga)
        ) %>%
        ggplot(mapping = aes(x = true_purity, y = ploidy_err, col = sample_fga_rank, group = sample, shape = sample)) +
        geom_point(size = 3) +
        geom_line() +
        facet_grid(coverage~tool) +
        scale_fill_viridis_c(name = "FGA Rank") +
        scale_color_viridis_c(name = "FGA Rank") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(x = "True Purity", y = "Delta Ploidy (True - Predicted)", fill = "FGA")

    p4 = df_metric %>%
        filter(tool != 'cnvkit') %>%
        group_by(sample) %>%
        dplyr::mutate(sample_fgs = mean(fgs, na.rm = TRUE)) %>%
        ungroup() %>%
        dplyr::mutate(
            ploidy_err =  true_ploidy - ploidy,
            sample_fgs_rank = rank(sample_fgs)
        ) %>%
        ggplot(mapping = aes(x = true_purity, y = ploidy_err, col = sample_fgs_rank, group = sample, shape = sample)) +
        geom_point(size = 3) +
        geom_line() +
        facet_grid(coverage~tool) +
        scale_fill_viridis_c(name = "FGS Rank") +
        scale_color_viridis_c(name = "FGS Rank") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        geom_hline(yintercept = 0, linetype = "dashed") +
        labs(x = "True Purity", y = "Delta Ploidy (True - Predicted)", fill = "FGS")

    p5 = df_metric %>%
        group_by(sample) %>%
        dplyr::mutate(sample_fga = mean(fga, na.rm = TRUE)) %>%
        ungroup() %>%
        dplyr::mutate(
            ploidy_err =  true_ploidy - ploidy,
            sample_fga_rank = rank(sample_fga)
        ) %>%
        ggplot(mapping = aes(x = true_purity, y = correctness_clonal, col = sample_fga_rank, group = sample, shape = sample)) +
        geom_point(size = 3) +
        geom_line() +
        facet_grid(coverage~tool) +
        scale_fill_viridis_c(name = "FGA Rank") +
        scale_color_viridis_c(name = "FGA Rank") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "True Purity", y = "% of correct CN", fill = "FGA")

    p6 = df_metric %>%
        group_by(sample) %>%
        dplyr::mutate(sample_fgs = mean(fgs, na.rm = TRUE)) %>%
        ungroup() %>%
        dplyr::mutate(
            ploidy_err =  true_ploidy - ploidy,
            sample_fgs_rank = rank(sample_fgs)
        ) %>%
        ggplot(mapping = aes(x = true_purity, y = correctness_clonal, col = sample_fgs_rank, group = sample, shape = sample)) +
        geom_point(size = 3) +
        geom_line() +
        facet_grid(coverage~tool) +
        scale_fill_viridis_c(name = "FGS Rank") +
        scale_color_viridis_c(name = "FGS Rank") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x = "True Purity", y = "% of correct CN", fill = "FGS")


    pp <- p1 + p2 + p3 + p4 + p5 +theme(legend.position = 'none') + p6 +theme(legend.position = 'none') + plot_layout(ncol=2, guides = 'collect')
    ggsave(filename = file.path(paste0(spn_id, '_v2.pdf')), plot = pp, dpi = 400, height = 12, width = 12, units = 'in')
    """
}