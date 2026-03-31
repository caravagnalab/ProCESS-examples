setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/assign_signature')
.libPaths("/orfeo/LTS/LADE/LT_storage/lvaleriani/R/x86_64-pc-linux-gnu-library/4.4/")

library(tidyverse)
library(optparse)

classify_expansion <- function(score, is_clonal) {
  case_when(
    is_clonal                          ~ 'Clonal',
    score >= .9                        ~ 'Strong\nExpansion',
    score >= .55 & score <.9           ~ 'Medium\nExpansion',
    score >= .2  & score <.55          ~ 'Low\nExpansion',
    score <  .2                        ~ 'Neutral\nExpansion'
  )
}
expansion_levels <- c('Clonal', 'Strong\nExpansion', 'Medium\nExpansion', 'Low\nExpansion', 'Neutral\nExpansion')
score_types <- c("all", 'no_tail', 'no_driver', 'no_sign')


option_list <- list(make_option(c("--cna_caller"), type = "character", default = 'sequenza'),
                    make_option(c("--vcf_caller"), type = "character", default = 'mutect2')
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

vcf_caller = opt$vcf_caller #"mutect2"
cna_caller = opt$cna_caller #"sequenza"
base = paste0('/orfeo/cephfs/scratch/cdslab/shared/SCOUT/interpretation/interpretation_', vcf_caller, "_", cna_caller, '/')

results <- readRDS(paste0(base, '/performance_cluster.rds'))

results_class <- results %>%
  mutate(across(
    all_of(paste0("score_", score_types, "_process")),
    ~ factor(classify_expansion(.x, cluster_process == 'Clonal'), levels = expansion_levels),
    .names = "{paste0('class_', str_extract(.col, '(?<=score_).+(?=_process)'), '_process')}"
  )) %>%
  mutate(across(
    all_of(paste0("score_", score_types, "_tool")),
    ~ factor(classify_expansion(.x, is_clonal_tool == TRUE), levels = expansion_levels),
    .names = "{paste0('class_', str_extract(.col, '(?<=score_).+(?=_tool)'), '_tool')}"
  ))

classes <- results_class %>%
  pull(class_all_process) %>%
  unique()

groups <- results_class %>%
  filter(!is.na(class_all_process)) %>%
  distinct(coverage, purity, sig_tool, sub_tool)

confusion_per_group_class <- map_dfr(score_types, function(st) {

  tool_col <- paste0("class_", st, "_tool")
  process_col <- "class_all_process"

  map_dfr(1:nrow(groups), function(i) {

    grp <- groups[i, ]

    df <- results_class %>%
      filter(
        !is.na(.data[[process_col]]), !is.na(.data[[tool_col]]),
        coverage == grp$coverage,
        purity   == grp$purity,
        sig_tool == grp$sig_tool,
        sub_tool == grp$sub_tool
      )

    map_dfr(classes, function(cls) {
      tibble(
        score_type = st,
        coverage   = grp$coverage,
        purity     = grp$purity,
        sig_tool   = grp$sig_tool,
        sub_tool   = grp$sub_tool,
        class      = cls,
        TP = sum( df[[process_col]] == cls &  df[[tool_col]] == cls),
        TN = sum( df[[process_col]] != cls &  df[[tool_col]] != cls),
        FP = sum( df[[process_col]] != cls &  df[[tool_col]] == cls),
        FN = sum( df[[process_col]] == cls &  df[[tool_col]] != cls)
      )
    })
  })
}) %>%
  mutate(
    Precision   = TP / (TP + FP),
    Recall      = TP / (TP + FN),
    Specificity = TN / (TN + FP),
    F1          = 2 * (Precision * Recall) / (Precision + Recall),
    Accuracy    = (TP + TN) / (TP + TN + FP + FN)
  )

saveRDS(object = confusion_per_group_class, file = paste0(base, '/confusion_cluster.rds'))


#mutations
mutations <- readRDS(paste0(base, '/performance_mutations.rds'))
mutations_class <- mutations %>%
  mutate(across(
    all_of(paste0("score_", score_types, "_process")),
    ~ factor(classify_expansion(.x, cluster_id_process == 'Clonal'), levels = expansion_levels),
    .names = "{paste0('class_', str_extract(.col, '(?<=score_).+(?=_process)'), '_process')}"
  )) %>%
  mutate(across(
    all_of(paste0("score_", score_types, "_tool")),
    ~ factor(classify_expansion(.x, is_clonal_tool == TRUE), levels = expansion_levels),
    .names = "{paste0('class_', str_extract(.col, '(?<=score_).+(?=_tool)'), '_tool')}"
  ))

saveRDS(object = mutations_class, file = paste0(base, '/performance_mutations_class.rds'))

classes <- mutations_class %>%
  pull(class_all_process) %>%
  unique()

groups <- mutations_class %>%
  filter(!is.na(class_all_process)) %>%
  distinct(coverage, purity, sig_tool, sub_tool, spn)

library(furrr)
wk = 24
plan(multisession, workers = wk)

confusion_per_group_class_mutations <- map_dfr(score_types, function(st) {
  message(st)

  tool_col    <- paste0("class_", st, "_tool")
  process_col <- "class_all_process"

  mutations_class %>%
    filter(!is.na(.data[[process_col]]), !is.na(.data[[tool_col]])) %>%
    select(coverage, purity, sig_tool, sub_tool, spn,
           all_of(c(process_col, tool_col))) %>%
    group_by(coverage, purity, sig_tool, sub_tool, spn) %>%
    group_modify(~ {
      df <- .x
      map_dfr(classes, function(cls) {
        tibble(
          class = cls,
          TP = sum( df[[process_col]] == cls &  df[[tool_col]] == cls),
          TN = sum( df[[process_col]] != cls &  df[[tool_col]] != cls),
          FP = sum( df[[process_col]] != cls &  df[[tool_col]] == cls),
          FN = sum( df[[process_col]] == cls &  df[[tool_col]] != cls)
        )
      })
    }) %>%
    ungroup() %>%
    mutate(score_type = st)

}) %>%
  mutate(
    Precision   = TP / (TP + FP),
    Recall      = TP / (TP + FN),
    Specificity = TN / (TN + FP),
    F1          = 2 * (Precision * Recall) / (Precision + Recall),
    Accuracy    = (TP + TN) / (TP + TN + FP + FN)
  )

saveRDS(object = confusion_per_group_class_mutations, 
        file = paste0(base, '/confusion_mutations.rds'))

