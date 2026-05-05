absolute_to_relative_coordinates <- function(muts, reference = CNAqc::chr_coordinates_GRCh38, centromere = F){
  vfrom = reference$from
  names(vfrom) = reference$chr
  if (!centromere){
    muts %>%
      mutate(
        start = start + vfrom[chr],
        end = end + vfrom[chr])
  } else if(centromere){
    muts %>%
      mutate(
        start = start + vfrom[chr],
        end = end + vfrom[chr],
        centromere = centromere + vfrom[chr])
  }
}

plot_qc <- function(cnaqc_list, type = 'simple_clonal'){
  if (type == 'simple_clonal'){
    table <- lapply(names(cnaqc_list), function(sample) {
      x <- cnaqc_list[[sample]]
      s <- stringr::str_split(sample, '_') %>% unlist()
      
      xqc = CNAqc:::compute_QC_table(x)
      QC_table = xqc$QC_table
      QC_table$sample = paste0(s[2], '_', s[3])
      QC_table$pPASS = xqc$percentage_PASS
      QC_table$NA_tests = xqc$NA_tests
      
      return(QC_table)
    }) %>% bind_rows()
    
    plot = table %>%
      ggplot2::ggplot(aes(x = karyotype, y = sample, fill = paste(QC))) +
      ggplot2::facet_wrap(type ~ .) +
      ggplot2::geom_tile(aes(width = .8, height = .8)) +
      ggplot2::scale_fill_manual(values = c(
        `PASS` = 'seagreen',
        `FAIL` = 'indianred3',
        `NA` = 'gainsboro'
      )) +
      CNAqc:::my_ggplot_theme() +
      ggplot2::labs(x = NULL, y = NULL, title = "Simple clonal CNAs ") +
      ggplot2::guides(fill = ggplot2::guide_legend('QC (NA not available)')) 
    
  } else if (type == 'complex_clonal'){
    table <- lapply(names(cnaqc_list), function(sample) {
      x <- cnaqc_list[[sample]]
      s <- str_split(sample, '_') %>% unlist()
      if (!is.null(x$peaks_analysis$general)){
        QC_table = x$peaks_analysis$general$expected_peak
        QC_table$sample = paste0(s[2], '_', s[3])
        return(QC_table)
      }
    }) %>% bind_rows()
    
    plot <- ggplot() + 
      ggplot2::geom_tile(
        data = table,
        ggplot2::aes(
          y = karyotype,
          x = multiplicity,
          fill = matched,
          width = .8,
          height = .8
        )
      ) +
      CNAqc:::my_ggplot_theme() +
      ggplot2::scale_fill_manual(values = c(`FALSE` = 'indianred3', `TRUE` = 'seagreen', `NA` = 'gainsboro')) + 
      ggplot2::labs(x = NULL, y = NULL, title = "Complex clonal CNAs ") +
      facet_grid(.~sample)
  }
  return(plot)
}

get_statistics_qc <- function(cnaqc_list, purity) {
  tt <- lapply(names(cnaqc_list), function(sample) {
    x <- cnaqc_list[[sample]]
    s <- sample
    
    tt = dplyr::tibble(
      sample = s,
      n_muts = x$n_mutations, 
      passed_muts = paste0(x$mutations %>% dplyr::filter(QC_PASS == TRUE) %>% nrow(), ' (', (round(x$mutations %>% dplyr::filter(QC_PASS == TRUE) %>% nrow() / x$mutations %>% nrow(), digits = 2))*100, '%)'), 
      passed_muts0 = x$mutations %>% dplyr::filter(QC_PASS == TRUE) %>% nrow(),
      n_clonal_cnas = x$n_cna_clonal,
      n_subclonal_cnas = x$n_cna_subclonal, 
      passed_clonal_cnas = paste0(x$cna %>% dplyr::filter(QC_PASS == TRUE) %>% nrow(), ' (', (round(x$cna %>% dplyr::filter(QC_PASS == TRUE) %>% nrow() / x$cna %>% nrow(), digits = 2))*100, '%)'), 
      passed_clonal_cnas0 = x$cna %>% dplyr::filter(QC_PASS == TRUE) %>% nrow(),
      true_purity = purity,
      caller_purity = x$purity, 
      purity_correction = round(x$peaks_analysis$score, 5),
      QC = x$peaks_analysis$QC,
    ) %>% t %>% 
      as.data.frame()
    
    tt = tt %>% 
      dplyr::mutate(info = rownames(tt)) %>% 
      dplyr::rename(values = V1) %>% 
      dplyr::mutate(gg = dplyr::case_when(info == "sample" ~ "sample", 
                                          str_ends(info, "muts") ~ "mutations", 
                                          str_ends(info, "cnas") ~ "cnas", 
                                          str_detect(info, "purity") ~ "purity")) %>% 
      dplyr::mutate(sample = s)
    
    return(tt)
  }) %>% dplyr::bind_rows() 
    
  tt_wider <- tt %>% tidyr::pivot_wider(values_from = values, names_from = sample)
  
  table <- tt_wider %>%
    dplyr::select(gg, info, everything()) %>% 
    dplyr::mutate(info = ifelse(info == "sample", "", info)) %>% 
    dplyr::mutate(gg = ifelse(info == "QC", "QC", gg)) %>% 
    dplyr::mutate(info = ifelse(info == "QC", " ", info)) %>% 
    dplyr::group_by(gg) %>% 
    gt::gt(row_group_as_column = T, process_md = T) %>% 
    gt::tab_header(title = gt::md("**QC summary**")) %>% 
    gt::tab_options(column_labels.hidden = T,  table.layout = "auto")
  
  return(list(plot = table, table = tt))
}





# qc <- ggplot() + 
#   ggplot2::geom_tile(
#     data = full_table,
#     ggplot2::aes(
#       y = '',
#       x = sample,
#       fill = QC,
#       width = .8,
#       height = .8
#     )
#   ) +
#   ggh4x::facet_nested(coverage+true_purity ~ spn, scales = 'free_x')  +
#   ggplot2::scale_fill_manual(values = c(`FAIL` = 'indianred3', `PASS` = 'seagreen', `NA` = 'gainsboro'))  +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   ylab('')
# 
# 
# 
# df_clean <- full_table %>%
#   mutate(
#     sample_id = sample,
#     n_clonal_cnas = as.numeric(n_clonal_cnas),
#     n_subclonal_cnas = as.numeric(n_subclonal_cnas),
#     passed_clonal_cnas = as.numeric(passed_clonal_cnas0),
#   ) %>%
#   mutate(
#     total_cnas = n_clonal_cnas,
#     passed_cnas = passed_clonal_cnas,
#     failed_cnas = total_cnas - passed_cnas,
#     pct_passed = passed_cnas / total_cnas * 100,
#     pct_failed = 100 - pct_passed
#   )
# 
# df_long <- df_clean %>%
#   select(sample_id, spn, pct_passed, pct_failed, true_purity, coverage) %>%
#   pivot_longer(cols = starts_with("pct_"), names_to = "status", values_to = "percentage") %>%
#   mutate(status = recode(status, pct_passed = "Passed", pct_failed = "Failed"))
# 
# 
# cnas <- ggplot(df_long, aes(x = sample_id, y = percentage, fill = status)) +
#   geom_bar(stat = "identity") +
#   ggh4x::facet_nested(coverage+true_purity ~ spn, scales = 'free_x')  +
#   labs(x = "sample", y = "Percentage of segment", fill = "Segment status") +
#   scale_fill_manual(values = c("Passed" = "seagreen", "Failed" = "gainsboro")) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))
# 


