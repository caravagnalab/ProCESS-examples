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
    
    purity_colors = c("#bfd3e6", "#8c96c6", "#810f7c")
    coverage_colors = c("#ccece6", "#66c2a4", "#006d2c")
    
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
            FP = caller_metrics[[sample_id]]\$detection_summary["False Positive"]
            d = caller_metrics[[sample_id]]\$performance_table %>% 
              dplyr::mutate(sample_id = sample_id) %>%  
              dplyr::mutate(purity = as.numeric(purity), coverage = as.numeric(coverage))
            d\$false_positive[1] = FP
            d
          }) %>% do.call("bind_rows", .) %>% 
            dplyr::mutate(caller = caller_name)
        }) %>% do.call("bind_rows", .)
        
        dplyr::bind_cols(parsed_metrics, params_grid[i,])
      }
    }) %>% do.call("bind_rows", .)
    
    plot_rep = function(df, spn_id) {
      # 1. Recall curve of SNVs, colored by tool, faceted by pi and cov
      p1 = df %>% 
        dplyr::filter(mut == "SNV") %>% 
        tidyr::replace_na(list(sensitivity = 0)) %>% 
        dplyr::group_by(CCF_bin, caller, mut, coverage, purity) %>% 
        dplyr::summarise(mean = mean(sensitivity), ylow = min(sensitivity), ymax = max(sensitivity)) %>%
        ggplot(mapping = aes(x = CCF_bin, y = mean, ymin=ylow, color = caller, linetype = mut)) +
        geom_line(aes(group = interaction(caller, mut)), size = 1) +
        geom_point(size = 2) +
        # facet_grid(coverage ~ purity, labeller = label_both) +
        ggh4x::facet_nested("Coverage" + coverage ~ "Purity" + purity, space = "free_x", scales = "free_x") +
        #facet_grid(coverage ~ purity, labeller = label_both, scales = "free_x", space = "free_x") +
        labs(x = "CCF Bin", y = "Recall", color = "Caller", linetype = "Mutation Type") +
        theme_bw() +
        scale_linetype_manual(values = c("INDEL" = "dotdash", "SNV" = "solid")) +
        scale_color_manual(values = method_colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(size = 8)) +
        scale_y_continuous(limits = c(0, 1))
      
      # 2. Recall curve of SNVs, colored by purity, faceted by tool and cov
      p2 = df %>% 
        dplyr::filter(mut == "SNV") %>% 
        tidyr::replace_na(list(sensitivity = 0)) %>% 
        dplyr::group_by(CCF_bin, caller, mut, coverage, purity) %>% 
        dplyr::summarise(mean = mean(sensitivity), ylow = min(sensitivity), ymax = max(sensitivity)) %>%
        ggplot(mapping = aes(x = CCF_bin, y = mean, ymin=ylow, color = as.factor(purity), linetype = mut)) +
        geom_line(aes(group = interaction(purity, mut)), size = 1) +
        geom_point(size = 2) +
        ggh4x::facet_nested("Coverage" + coverage ~ "Caller" + caller, space = "free_x", scales = "free_x") +
        labs(x = "CCF Bin", y = "Recall", color = "Purity", linetype = "Mutation Type") +
        theme_bw() +
        scale_linetype_manual(values = c("INDEL" = "dotdash", "SNV" = "solid")) +
        scale_color_manual(values = purity_colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(size = 8)) +
        scale_y_continuous(limits = c(0, 1))
      
      # 3. Recall curve of SNVs, colored by coverage, faceted by tool and purity
      p3 = df %>% 
        dplyr::filter(mut == "SNV") %>% 
        tidyr::replace_na(list(sensitivity = 0)) %>% 
        dplyr::group_by(CCF_bin, caller, mut, coverage, purity) %>% 
        dplyr::summarise(mean = mean(sensitivity), ylow = min(sensitivity), ymax = max(sensitivity)) %>%
        ggplot(mapping = aes(x = CCF_bin, y = mean, ymin=ylow, color = as.factor(coverage), linetype = mut)) +
        geom_line(aes(group = interaction(coverage, mut)), size = 1) +
        geom_point(size = 2) +
        ggh4x::facet_nested("Purity" + purity ~ "Caller" + caller, space = "free_x", scales = "free_x") +
        labs(x = "CCF Bin", y = "Recall", color = "Coverage", linetype = "Mutation Type") +
        theme_bw() +
        scale_linetype_manual(values = c("INDEL" = "dotdash", "SNV" = "solid")) +
        scale_color_manual(values = coverage_colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(size = 8)) +
        scale_y_continuous(limits = c(0, 1))
      
      # 4. Recall curve of INDELs, colored by tool, faceted by pi and cov
      p4 = df %>% 
        dplyr::filter(mut == "INDEL") %>% 
        tidyr::replace_na(list(sensitivity = 0)) %>% 
        dplyr::group_by(CCF_bin, caller, mut, coverage, purity) %>% 
        dplyr::summarise(mean = mean(sensitivity), ylow = min(sensitivity), ymax = max(sensitivity)) %>%
        ggplot(mapping = aes(x = CCF_bin, y = mean, ymin=ylow, color = caller, linetype = mut)) +
        geom_line(aes(group = interaction(caller, mut)), size = 1) +
        geom_point(size = 2) +
        # facet_grid(coverage ~ purity, labeller = label_both) +
        ggh4x::facet_nested("Coverage" + coverage ~ "Purity" + purity, space = "free_x", scales = "free_x") +
        #facet_grid(coverage ~ purity, labeller = label_both, scales = "free_x", space = "free_x") +
        labs(x = "CCF Bin", y = "Recall", color = "Caller", linetype = "Mutation Type") +
        theme_bw() +
        #scale_linetype_manual(values = c("INDEL" = "dotdash", "SNV" = "solid")) +
        scale_color_manual(values = method_colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(size = 8)) +
        scale_y_continuous(limits = c(0, 1))
      
      # 5. Recall curve of SNVs, colored by purity, faceted by tool and cov
      p5 = df %>% 
        dplyr::filter(mut == "INDEL") %>% 
        tidyr::replace_na(list(sensitivity = 0)) %>% 
        dplyr::group_by(CCF_bin, caller, mut, coverage, purity) %>% 
        dplyr::summarise(mean = mean(sensitivity), ylow = min(sensitivity), ymax = max(sensitivity)) %>%
        ggplot(mapping = aes(x = CCF_bin, y = mean, ymin=ylow, color = as.factor(purity), linetype = mut)) +
        geom_line(aes(group = interaction(purity, mut)), size = 1) +
        geom_point(size = 2) +
        ggh4x::facet_nested("Coverage" + coverage ~ "Caller" + caller, space = "free_x", scales = "free_x") +
        labs(x = "CCF Bin", y = "Recall", color = "Purity", linetype = "Mutation Type") +
        theme_bw() +
        #scale_linetype_manual(values = c("INDEL" = "dotdash", "SNV" = "solid")) +
        scale_color_manual(values = purity_colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(size = 8)) +
        scale_y_continuous(limits = c(0, 1))
      
      # 6. Recall curve of SNVs, colored by coverage, faceted by tool and purity
      p6 = df %>% 
        dplyr::filter(mut == "INDEL") %>% 
        tidyr::replace_na(list(sensitivity = 0)) %>% 
        dplyr::group_by(CCF_bin, caller, mut, coverage, purity) %>% 
        dplyr::summarise(mean = mean(sensitivity), ylow = min(sensitivity), ymax = max(sensitivity)) %>%
        ggplot(mapping = aes(x = CCF_bin, y = mean, ymin=ylow, color = as.factor(coverage), linetype = mut)) +
        geom_line(aes(group = interaction(coverage, mut)), size = 1) +
        geom_point(size = 2) +
        ggh4x::facet_nested("Purity" + purity ~ "Caller" + caller, space = "free_x", scales = "free_x") +
        labs(x = "CCF Bin", y = "Recall", color = "Coverage", linetype = "Mutation Type") +
        theme_bw() +
        #scale_linetype_manual(values = c("INDEL" = "dotdash", "SNV" = "solid")) +
        scale_color_manual(values = coverage_colors) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              strip.text = element_text(size = 8)) +
        scale_y_continuous(limits = c(0, 1))
      
      # TP_plot = df %>% 
      #   dplyr::group_by(caller, mut, coverage, purity) %>% 
      #   dplyr::summarise(TP = sum(true_positives, na.rm = TRUE)) %>% 
      #   ggplot(mapping = aes(x = caller, y = TP, fill = as.factor(coverage))) +
      #   geom_col(position = "dodge") +
      #   ggh4x::facet_nested(mut~"Purity"+purity, scales = "free_y") +
      #   theme_bw() +
      #   scale_fill_manual(values = purity_colors) +
      #   labs(x = "Caller", y = "Number of correctly detected mutations", 
      #        fill = "Coverage")
      
      # FP_plot = df %>% 
      #   dplyr::filter(caller != "freebayes") %>% 
      #   dplyr::group_by(caller, mut, coverage, purity) %>% 
      #   dplyr::summarise(FP = sum(false_positive, na.rm = TRUE)) %>% 
      #   ggplot(mapping = aes(x = caller, y = FP, fill = as.factor(coverage))) +
      #   geom_col(position = "dodge") +
      #   ggh4x::facet_nested(mut~"Purity"+purity, scales = "free_y") +
      #   theme_bw() +
      #   scale_fill_manual(values = purity_colors) +
      #   labs(x = "Caller", y = "Number of false positive mutations", 
      #        fill = "Coverage")
      
      TP_plot = df %>% 
        dplyr::group_by(caller, mut, coverage, purity) %>% 
        dplyr::summarise(TP = sum(true_positives, na.rm = TRUE), total = sum(total_truth, na.rm = T)) %>% 
        ggplot(mapping = aes(x = factor(coverage), y = TP / total, fill = caller)) +
        geom_col(position = "dodge") +
        ggh4x::facet_nested(mut~"Purity"+purity, scales = "free_y") +
        theme_bw() +
        scale_fill_manual(values = method_colors) +
        labs(x = "Coverage", y = "True Positive Rate", 
             fill = "Caller")
      
      FP_plot = df %>% 
        dplyr::filter(!(caller == "freebayes" & mut == "INDEL")) %>% 
        dplyr::group_by(caller, mut, coverage, purity) %>% 
        dplyr::summarise(FP = sum(false_positive, na.rm = TRUE), TP = sum(true_positives, na.rm = TRUE)) %>% 
        ggplot(mapping = aes(x = factor(coverage), y = FP / (TP + FP), fill = caller)) +
        geom_col(position = "dodge") +
        ggh4x::facet_nested(mut~"Purity"+purity, scales = "free_y") +
        theme_bw() +
        scale_fill_manual(values = method_colors) +
        labs(x = "Coverage", y = "False Discovery Rate", fill = "Caller")
      
      design <- "
        AABBCC
        AABBCC
        DDEEFF
        DDEEFF
        #GGHH#
        #GGHH#
        "
      
      title <- paste0(spn_id, " - somatic mutations")
      
      report_plot <- patchwork::free(p1) + patchwork::free(p2) + patchwork::free(p3) +
        patchwork::free(p4) + patchwork::free(p5) + patchwork::free(p6) +
        patchwork::free(TP_plot) + patchwork::free(FP_plot) + 
        patchwork::plot_layout(design = design) +
        patchwork::plot_annotation(title, tag_levels = list(c("SNV", "", "", "INDEL", "", "", "TPR and FDR"))) &
        ggplot2::theme(text = ggplot2::element_text(size = 12), 
                       legend.position = "bottom", 
                       legend.direction = "horizontal", legend.box = "vertical", legend.spacing.y = unit(1, "pt"))
      
      report_plot
    }
    
    final_report = plot_rep(df, spn_id)
    ggsave("final_report.pdf", plot = final_report, width = 18, height = 18)

    """
}




