process GENERATE_SOMATIC_REPORT_ALL_COMBINATION {
    tag "${meta.spn}_somatic"

    input:
    tuple val(meta), val(rds_list), val(combinations)

    output:
    tuple val(meta), path("*.pdf"),	emit:report
    
    publishDir "${params.outdir}/${meta.spn}/somatic/report/", mode: 'copy'

    script:
    """
    #!/usr/bin/env Rscript
    
    rm(list = ls())
    options(bitmapType='cairo')
    require(tidyverse)
    library(optparse)
    
    source("${projectDir}/bin/getters/sarek_getters.R")
    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/somatic/utils/utils.R")
    source("${projectDir}/bin/somatic/utils/plot_utils.R")

    spn_id = "${meta.spn}"
    gender = "${meta.sex}"

    combination = substr("$combinations", 2, nchar("$combinations")-1) 
    combination = strsplit(combination, ", ")[[1]]

    data =  strsplit("${rds_list}", '], ') %>% unlist() %>% 
            stringr::str_replace_all('\\\\[', '') %>% 
            stringr::str_replace_all('\\\\]', '') %>% 
            strsplit(', ')
    names(data) = combination

    MUT_TYPES = c("INDEL", "SNV")
    params_grid = expand.grid(combination, MUT_TYPES)
    colnames(params_grid) = c("comb", "mut")

    i = 1
    df = lapply(1:nrow(params_grid), function(i) {
        mut_type = params_grid[i,]\$mut
        combination = params_grid[i,]\$comb
        tmp = strsplit(as.character(combination), '_') %>% unlist()
        coverage = tmp[[1]]
        purity = tmp[[2]]
        
        file = data[[combination]]
        file =  file[grep(mut_type, file)]
        
        if (file.exists(file)) {
            metrics = readRDS(file)
            
            parsed_metrics = lapply(names(metrics), function(caller_name) {
                caller_metrics = metrics[[caller_name]]
                sample_names = names(caller_metrics)
                lapply(sample_names, function(sample_id) {
                    caller_metrics[[sample_id]]\$performance_table %>% 
                    dplyr::mutate(sample_id = sample_id) %>%  
                    dplyr::mutate(purity = purity, coverage = coverage)
            }) %>% do.call("bind_rows", .) %>% 
                dplyr::mutate(caller = caller_name)
            }) %>% do.call("bind_rows", .)
            
            dplyr::bind_cols(parsed_metrics, params_grid[i,])
        }
    }) %>% do.call("bind_rows", .)

    plot_rep = function(df, spn_id) {
        # 1. Precision comparison with mut distinction via linetype
        p1 = df %>% 
            tidyr::replace_na(list(precision = 0)) %>% 
            dplyr::group_by(CCF_bin, caller, mut, coverage, purity) %>% 
            dplyr::summarise(mean = mean(precision), ylow = min(precision), ymax = max(precision)) %>%
            ggplot(mapping = aes(x = CCF_bin, y = mean, ymin=ylow,ymax=ymax, color = caller, linetype = mut)) +
            geom_line(aes(group = interaction(caller, mut)), size = 1) +
            geom_point(size = 2) +
            facet_grid(coverage ~ purity, labeller = label_both) +
            labs(x = "CCF Bin", y = "Precision",
                color = "Caller", linetype = "Mutation Type") +
            theme_bw() +
            scale_linetype_manual(values = c("INDEL" = "dotdash", "SNV" = "solid")) +
            scale_color_manual(values = method_colors) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                strip.text = element_text(size = 8)) +
            scale_y_continuous(limits = c(0, 1))
        
        TP_plot = df %>% 
            dplyr::group_by(caller, mut, coverage, purity) %>% 
            dplyr::summarise(TP = sum(true_positives, na.rm = TRUE)) %>% 
            ggplot(mapping = aes(x = caller, y = TP, fill = coverage)) +
            geom_col(position = "dodge") +
            ggh4x::facet_nested(mut~"Purity"+purity, scales = "free_y") +
            theme_bw() +
            scale_fill_manual(values = c("#b3cde3", "#8c96c6")) +
            labs(x = "Caller", y = "Number of correctly detected mutations", 
                fill = "Coverage")
        
        # 2. Sensitivity/Recall comparison with mut distinction via linetype
        p2 = df %>% 
            tidyr::replace_na(list(sensitivity = 0)) %>% 
            dplyr::group_by(CCF_bin, caller, mut, coverage, purity) %>% 
            dplyr::summarise(mean = mean(sensitivity), ylow = min(sensitivity), ymax = max(sensitivity)) %>%
            ggplot(mapping = aes(x = CCF_bin, y = mean, ymin=ylow, color = caller, linetype = mut)) +
            geom_line(aes(group = interaction(caller, mut)), size = 1) +
            geom_point(size = 2) +
            facet_grid(coverage ~ purity, labeller = label_both) +
            labs(x = "CCF Bin", y = "Sensitivity",
                color = "Caller", linetype = "Mutation Type") +
            theme_bw() +
            scale_linetype_manual(values = c("INDEL" = "dotdash", "SNV" = "solid")) +
            scale_color_manual(values = method_colors) +
            theme(axis.text.x = element_text(angle = 45, hjust = 1),
                strip.text = element_text(size = 8)) +
            scale_y_continuous(limits = c(0, 1))
        
        # 5. Combined metric plot with mut distinction via shape
        p5 <- ggplot(df, aes(x = sensitivity, y = precision, color = caller, shape = mut)) +
            geom_point(aes(size = ccf_correlation), alpha = 0.7) +
            facet_grid(coverage ~ purity, labeller = label_both) +
            labs(title = "Precision vs Sensitivity",
                x = "Sensitivity", y = "Precision",
                color = "Caller", size = "VAF Corr", shape = "Mutation Type", caption = "Each point represents average across samples") +
            scale_color_manual(values = method_colors) +
            theme_bw() +
            scale_x_continuous(limits = c(0, 1)) +
            scale_y_continuous(limits = c(0, 1)) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5)
        
        p3 = df %>% 
            tidyr::replace_na(list(sensitivity = 0, precision = 0)) %>% 
            dplyr::filter(mut == "SNV", !(CCF_bin %in% c("0-5%","5-10%"))) %>% 
            dplyr::group_by(caller, mut, purity, sample_id, coverage) %>% 
            dplyr::summarise(sensitivity = mean(sensitivity), precision = mean(precision), .groups = "drop") %>%
            # Arrange by purity to ensure proper arrow direction
            dplyr::arrange(caller, sample_id, coverage, purity) %>%
            ggplot(mapping = aes(x = sensitivity, y = precision, color = caller, shape = sample_id)) +
            geom_point(alpha = 1, size = 5) +
            # Add arrows connecting points from low to high purity
            geom_path(aes(group = interaction(caller, sample_id)), 
                    col = "gray30",
                    arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                    alpha = 0.6, size = 0.5) +
            facet_wrap(~ coverage, labeller = label_both, scales = "free") +
            labs(title = "Precision vs Sensitivity", subtitle = "SNV", caption = "Arrows show purity progression",
                x = "Sensitivity", y = "Precision",
                color = "Caller", shape = "Sample") +
            scale_color_manual(values = method_colors) +
            theme_bw() +
            # Ensure square aspect ratio
            scale_x_continuous(limits = c(NA, 1)) +
            scale_y_continuous(limits = c(NA, 1)) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5)
        
        # p3 = df %>% 
        #   tidyr::replace_na(list(sensitivity = 0, precision = 0)) %>% 
        #   dplyr::filter(mut == "SNV", !(CCF_bin %in% c("0-5%","5-10%"))) %>% 
        #   dplyr::group_by(caller, mut, purity, sample_id, coverage) %>% 
        #   dplyr::summarise(sensitivity = mean(sensitivity), precision = mean(precision)) %>%
        #   ggplot(mapping = aes(x = sensitivity, y = precision, color = caller, shape = sample_id)) +
        #   geom_point(alpha = 0.7, size = 3) +
        #   facet_grid(coverage ~ purity, labeller = label_both) +
        #   labs(title = "Precision vs Sensitivity", subtitle = "SNV",
        #        x = "Sensitivity", y = "Precision",
        #        color = "Caller", size = "VAF Corr", shape = "Sample") +
        #   scale_color_manual(values = method_colors) +
        #   theme_bw() +
        #   scale_x_continuous(limits = c(0, 1)) +
        #   scale_y_continuous(limits = c(0, 1)) +
        #   geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5)
        
        p4 = df %>% 
            tidyr::replace_na(list(sensitivity = 0, precision = 0)) %>% 
            dplyr::filter(mut == "INDEL", !(CCF_bin %in% c("0-5%","5-10%"))) %>% 
            dplyr::group_by(caller, mut, purity, sample_id, coverage) %>% 
            dplyr::summarise(sensitivity = mean(sensitivity), precision = mean(precision)) %>%
            ggplot(mapping = aes(x = sensitivity, y = precision, color = caller, shape = sample_id)) +
            geom_point(alpha = 0.7, size = 3) +
            facet_grid(coverage ~ purity, labeller = label_both) +
            labs(title = "INDEL",
                x = "Sensitivity", y = "Precision",
                color = "Caller", size = "VAF Corr", shape = "Sample") +
            scale_color_manual(values = method_colors) +
            theme_bw() +
            scale_x_continuous(limits = c(0, 1)) +
            scale_y_continuous(limits = c(0, 1)) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5)
        
        p4 = df %>% 
            tidyr::replace_na(list(sensitivity = 0, precision = 0)) %>% 
            dplyr::filter(mut == "INDEL", !(CCF_bin %in% c("0-5%","5-10%"))) %>% 
            dplyr::group_by(caller, mut, purity, sample_id, coverage) %>% 
            dplyr::summarise(sensitivity = mean(sensitivity), precision = mean(precision), .groups = "drop") %>%
            # Arrange by purity to ensure proper arrow direction
            dplyr::arrange(caller, sample_id, coverage, purity) %>%
            ggplot(mapping = aes(x = sensitivity, y = precision, color = caller, shape = sample_id)) +
            geom_point(alpha = 1, size = 5) +
            # Add arrows connecting points from low to high purity
            geom_path(aes(group = interaction(caller, sample_id)), 
                    col = "gray30",
                    arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
                    alpha = 0.6, size = 0.5) +
            facet_wrap(~ coverage, labeller = label_both, scales = "free") +
            labs(title = "INDEL", caption = "Arrows show purity progression",
                x = "Sensitivity", y = "Precision",
                color = "Caller", shape = "Sample") +
            scale_color_manual(values = method_colors) +
            theme_bw() +
            # Ensure square aspect ratio
            scale_x_continuous(limits = c(0, 1)) +
            scale_y_continuous(limits = c(0, 1)) +
            geom_abline(intercept = 0, slope = 1, linetype = "dashed", alpha = 0.5)
        
        # 6. Performance summary with mut as additional faceting dimension
        # Reshape data for multiple metrics
        # df_long <- df %>%
        #   select(caller, coverage, purity, mut, CCF_bin, precision, sensitivity, 
        #          vaf_correlation) %>%
        #   tidyr::pivot_longer(cols = c(precision, sensitivity, vaf_correlation),
        #                       names_to = "metric", values_to = "value")
        # 
        # # Alternative version of p6 with mut as fill instead of facet (if too many facets)
        # p6 <- ggplot(df_long, aes(x = caller, y = value, fill = mut)) +
        #   geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) +
        #   facet_grid(metric ~ paste("Cov:", coverage, "Pur:", purity), 
        #              scales = "free_y", labeller = label_value) +
        #   labs(title = "Performance Metrics Distribution by Caller and Mutation Type",
        #        x = "Caller", y = "Metric Value", fill = "Mutation Type") +
        #   theme_bw() +
        #   scale_fill_manual(values = c("INDEL" = "steelblue", "SNV" = "darkorange")) +
        #   theme(axis.text.x = element_text(angle = 45, hjust = 1),
        #         strip.text = element_text(size = 7))
        
        design <- "
        AAABBB
        AAABBB
        AAABBB
        AAABBB
        CCDDEE
        CCDDEE
        CCDDEE
        "
        
        title <- paste0(spn_id, " - somatic mutations")
        report_plot <- patchwork::free(p1) + patchwork::free(p2) +
            patchwork::free(p3) + patchwork::free(p4) + patchwork::free(TP_plot) + 
            patchwork::plot_layout(design = design) +
            patchwork::plot_annotation(title) & 
            ggplot2::theme(text = ggplot2::element_text(size = 12), 
                        legend.position = "bottom", 
                        legend.direction = "horizontal", legend.box = "vertical", legend.spacing.y = unit(1, "pt"))
        
        report_plot
    }

    final_report = plot_rep(df, spn_id)
    ggsave("final_report.pdf", plot = final_report, width = 18, height = 14)

    """
}




