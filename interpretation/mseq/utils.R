signature_colors = c('#f1696bff', 
                     '#8fbd8cff', 
                     '#87c7d6ff', 
                     '#bac3deff', 
                     '#d7bfd9ff', 
                     '#a8a2a1ff', 
                     '#cfadb3ff', 
                     '#3c609aff', 
                     '#9a4564ff', 
                     '#fbcb5bff', 
                     '#c2b280ff', 
                     '#d47e2dff', 
                     '#5f8676ff',
                     'forestgreen',
                     'orange',
                     'brown4')

colors_cluster = c('indianred', 
                   'steelblue', 
                   'forestgreen', 
                   'goldenrod', 
                   'darkorange3', 
                   'palevioletred', 
                   'mediumpurple', 
                   'cornsilk4', 
                   'olivedrab3', 
                   'steelblue4', 
                   'indianred4',
                   'aquamarine3',
                   'saddlebrown',
                   'deeppink2',
                   'cornflowerblue',
                   'black')

names(colors_cluster) = paste0('C',0:15)

plot_exposure <- function(df, sig_type, colors){
  plot <- df %>% 
    ggplot(aes(fill=Signature, y=Exposure, x=as.factor(Samples))) + 
    geom_bar(position="fill", stat="identity")+
    coord_flip()+
    scale_fill_manual(values=colors)+
    theme_bw()+
    xlab("Cluster")+
    ylab("Exposures")+
    ggtitle(label = paste(sig_type))
  return(plot)
}

cosine_similarity <- function(vec1, vec2) {
  sum(vec1 * vec2) / (sqrt(sum(vec1^2)) * sqrt(sum(vec2^2)))
}

compare_signatures <- function(df1, df2) {
  
  df1_f <- df1 %>% filter(Exposure>.05)
  df2_f <- df2 %>% filter(Exposure>.05)
  
  # Check if signatures are identical
  sigs_match <- setequal(df1_f$Signature, df2_f$Signature)
  ndiff_1 <- length(setdiff(df2_f$Signature, df1_f$Signature))
  ndiff_2 <- length(setdiff(df1_f$Signature, df2_f$Signature))
  ndiff = sum(ndiff_1 + ndiff_2)
  ntot <- length(c(df1_f$Signature, df2_f$Signature) %>% unique())
  
  comparison <- full_join(
    df1 %>% select(Signature, Exposure, n_sig),
    df2 %>% select(Signature, Exposure, n_sig),
    by = "Signature",
    suffix = c("_tmp", "_clonal")
  ) %>%
    mutate(
      diff_exposure = Exposure_tmp - Exposure_clonal,
      diff_n_sig = n_sig_tmp - n_sig_clonal
    )%>% 
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  comparison_data <- full_join(
    df1, 
    df2, 
    by = "Signature", 
    suffix = c("_tmp", "_clonal")
  ) %>% 
    select(Signature, Exposure_tmp, Exposure_clonal, n_sig_tmp, n_sig_clonal) %>% 
    mutate(across(where(is.numeric), ~replace_na(., 0)))
  
  # 4. Calculate similarities
  cos_sim_exposure <- cosine_similarity(
    comparison_data$Exposure_tmp, 
    comparison_data$Exposure_clonal
  )
  
  
  
  return(list('df' = comparison,
              'match' = sigs_match, 
              'cs_exp' = cos_sim_exposure, 
              'n_diff' = ndiff,
              'n_tot' = ntot))
}