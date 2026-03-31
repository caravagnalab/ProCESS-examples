# Sankey plot #

generate_sankey_sbs <- function(sankey_df, spn_id, cov = NULL, pur = NULL) {
  color <- COSMIC_color_palette()#, get_colors_for(sankey_df$Signature %>% unique(), pal_name = 'Set3')
  df_spn <- sankey_df %>% filter(SPN == spn_id)
  
  if (!is.null(cov)) df_spn <- df_spn %>% filter(Coverage == paste0('coverage_', cov))
  if (!is.null(pur)) df_spn <- df_spn %>% filter(Purity == paste0('purity_', pur))
  
  # Prepare data
  long_df <- df_spn %>%
    mutate(
      Method = factor(Method, levels = c("ProCESS", "BASCULE", "SigProfiler")),
      alluvium = interaction(Sample_ID, Signature)
    )
  
  # Label data: calculate midpoint y-position for each stratum (Signature) in each Method
  label_data <- long_df %>%
    group_by(Sample_ID, Method, Signature) %>%
    summarise(Exposure = sum(Exposure), .groups = "drop") %>%
    arrange(Sample_ID, Method, desc(Signature)) %>%  # ordering for stacking
    group_by(Sample_ID, Method) %>%
    mutate(
      cum_exposure = cumsum(Exposure),
      ymin = cum_exposure - Exposure,
      ymid = ymin + Exposure / 2,
      Method = factor(Method, levels = c("ProCESS", "BASCULE", "SigProfiler")),
      label = paste0(Signature, "\n(", round(Exposure, 2), ")")
    ) %>%
    ungroup()
  
  title_text <- paste0(spn_id)
  if (!is.null(cov)) title_text <- paste0(title_text, ", cov=", cov)
  if (!is.null(pur)) title_text <- paste0(title_text, ", pur=", pur)
  
  
  ggplot(long_df,
         aes(x = factor(as.factor(Method), levels = c('SigProfiler', 'ProCESS', 'BASCULE')),
             stratum = Signature, alluvium = alluvium,
             y = Exposure, fill = Signature)) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "gray", alpha = 0.6) +
    geom_stratum(width = 1/3, color = "black") +
    # geom_text(
    #   data = label_data,
    #   aes(x = Method, y = ymid, label = label),
    #   inherit.aes = FALSE,
    #   size = 3,
    #   color = "black"
    # ) +
    scale_fill_manual(values = color,drop=FALSE) +
    facet_wrap(~ Sample_ID, nrow = 1, scales = "free_y") +
    theme_minimal() +
    theme(
      legend.position = "bottom"
      # axis.text.x = element_text(size = 12),
      # axis.text.y = element_blank(),
      # axis.title = element_text(size = 12),
      # plot.title = element_text(size = 12, face = "plain")
    ) +
    labs(
      title = title_text,
      x = "Method",
      y = "Exposure"
    )
}

COSMIC_color_palette = function(seed=50,
                                cosmic_process="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/GRCh38/") {
  COSMIC_sbs <- read.table(paste0(cosmic_process,"SBS_signatures.txt"),header = T) %>% colnames()
  COSMIC_indels <-read.table(paste0(cosmic_process,"indel_signatures.txt"),header = T) %>% colnames()
  catalogues = list(COSMIC_sbs[2:length(COSMIC_sbs)],COSMIC_indels[2:length(COSMIC_indels)])
  sigs =  unlist(catalogues) %>% unique()
  set.seed(seed)
  return(Polychrome::createPalette(length(sigs), c("#6B8E23","#4169E1"), M=1000,
                                   target="normal", range=c(15,80)) %>%
           setNames(sigs))
}

get_colors_for <- function(values, pal_name = "Set3") {
  colors <- NULL
  if (length(values) > 0) {
    num_of_values <- values %>% length()
    if (num_of_values < 3) {
      colors <- RColorBrewer::brewer.pal(3, pal_name)
      colors <- colors[seq_len(num_of_values)]
    } else {
      colors <- RColorBrewer::brewer.pal(num_of_values, pal_name)
    }
    names(colors) <- values
  }
  return(colors)
}

generate_sankey_id <- function(sankey_df, spn_id, cov = NULL, pur = NULL) {
  color <- COSMIC_color_palette()#, get_colors_for(sankey_df$Signature %>% unique(), pal_name = 'Set3')
  df_spn <- sankey_df %>% filter(SPN == spn_id)

  if (!is.null(cov)) df_spn <- df_spn %>% filter(Coverage == paste0('coverage_', cov))
  if (!is.null(pur)) df_spn <- df_spn %>% filter(Purity == paste0('purity_', pur))

  # Prepare data
  long_df <- df_spn %>%
    mutate(
      Method = factor(Method, levels = c("BASCULE", "ProCESS","SigProfiler")),
      alluvium = interaction(Sample_ID, Signature)
    )

  # Label data: calculate midpoint y-position for each stratum (Signature) in each Method
  label_data <- long_df %>%
    group_by(Sample_ID, Method, Signature) %>%
    summarise(Exposure = sum(Exposure), .groups = "drop") %>%
    arrange(Sample_ID, Method, desc(Signature)) %>%  # ordering for stacking
    group_by(Sample_ID, Method) %>%
    mutate(
      cum_exposure = cumsum(Exposure),
      ymin = cum_exposure - Exposure,
      ymid = ymin + Exposure / 3,
      Method = factor(Method, levels = c("SigProfiler","ProCESS","BASCULE")),
      label = paste0(Signature, "\n(", round(Exposure, 2), ")")
    ) %>%
    ungroup()

  title_text <- paste0(spn_id)
  if (!is.null(cov)) title_text <- paste0(title_text, ", cov=", cov)
  if (!is.null(pur)) title_text <- paste0(title_text, ", pur=", pur)

  ggplot(long_df,
         aes(x = Method,
             stratum = Signature, alluvium = alluvium,
             y = Exposure, fill = Signature)) +
    geom_flow(stat = "alluvium", lode.guidance = "frontback", color = "gray", alpha = 0.6) +
    geom_stratum(width = 1/3, color = "black") +
    # geom_text(
    #   data = label_data,
    #   aes(x = Method, y = ymid, label = label),
    #   inherit.aes = FALSE,
    #   size = 3,
    #   color = "black"
    # ) +
    scale_fill_manual(values = color) +
    facet_wrap(~ Sample_ID, ncol = 2, scales = "free_y") +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_blank(),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 12, face = "plain")
    ) +
    labs(
      title = title_text,
      x = "Method",
      y = "Exposure"
    )
}

tool_col = c('SigProfiler' = 'steelblue', 'SparseSignatures' = 'orange')
plot_cosine_similarity <- function(data,
                                   x = "Tool",
                                   y = "CosineSimilarity",
                                   color_var = "Sample",
                                   shape_var = "Purity",
                                   facet_row = "Coverage",
                                   facet_col = "SPN",
                                   shape_values = c(16, 17, 15, 18, 3, 4),
                                   title = "Sample-wise Cosine similarity") {

  # Check required columns
  required_cols <- c(x, y, color_var, shape_var, facet_row, facet_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if(length(missing_cols) > 0){
    stop("Data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Convert shape_var to factor
  data <- data %>%
    mutate(!!shape_var := factor(.data[[shape_var]]))

  xq <- sym(x)
  yq <- sym(y)
  colorq <- sym(color_var)
  shapeq <- sym(shape_var)
  facet_row_q <- sym(facet_row)
  facet_col_q <- sym(facet_col)

  # Use as_labeller for Coverage to prepend "Coverage: "
  coverage_labeller <- as_labeller(function(x) paste0("Coverage: ", x))

  ggplot(data, aes(x = !!xq, y = !!yq, color = !!colorq, shape = !!shapeq)) +
    geom_jitter(width = 0.2, height = 0, alpha = 0.9, size = 5) +
    facet_grid(rows = vars(!!facet_row_q), cols = vars(!!facet_col_q),
               labeller = labeller(
                 !!facet_col := label_value,
                 !!facet_row := coverage_labeller
               )) +
    scale_shape_manual(values = shape_values) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 16, face = "plain"),
      strip.text.x = element_text(size = 16, face = "bold"),
      strip.text.y = element_text(size = 16, face = "plain")
    ) +
    labs(
      title = title,
      x = x,
      y = y,
      color = color_var,
      shape = shape_var
    ) +
    guides(color = guide_legend(override.aes = list(size = 3)),
           shape = guide_legend(override.aes = list(size = 3)))
}


plot_mse_per_signature <- function(data,
                                   x = "Signature",
                                   y = "MSE",
                                   color_var = "Tool",
                                   shape_var = "Purity",
                                   facet_row = "Coverage",
                                   facet_col = "SPN",
                                   shape_values = c(16, 17, 15, 18, 3, 4),
                                   title = "MSE per Signature") {

  # Check columns
  required_cols <- c(x, y, color_var, shape_var, facet_row, facet_col)
  missing_cols <- setdiff(required_cols, colnames(data))
  if(length(missing_cols) > 0){
    stop("Data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  data <- data %>%
    mutate(!!shape_var := factor(.data[[shape_var]]))

  xq <- sym(x)
  yq <- sym(y)
  colorq <- sym(color_var)
  shapeq <- sym(shape_var)
  facet_row_q <- sym(facet_row)
  facet_col_q <- sym(facet_col)
  
  if (color_var == 'Tool'){
    color = tool_col
  }

  coverage_labeller <- as_labeller(function(x) paste0("Coverage: ", x))

  ggplot(data, aes(x = !!xq, y = !!yq, color = !!colorq, shape = !!shapeq)) +
    geom_jitter(width = 0.2, height = 0, size = 5, alpha = 0.9) +
    facet_grid(rows = vars(!!facet_row_q), cols = vars(!!facet_col_q),
               labeller = labeller(
                 !!facet_col := label_value,
                 !!facet_row := coverage_labeller
               ), scales = 'free_x') +
    scale_shape_manual(values = shape_values) +
    scale_color_manual(values = color) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 16, face = "plain"),
      strip.text.x = element_text(size = 16, face = "bold"),
      strip.text.y = element_text(size = 16, face = "plain")
    ) +
    labs(
      title = title,
      x = x,
      y = y,
      color = color_var,
      shape = shape_var
    ) +
    guides(color = guide_legend(override.aes = list(size = 4)),
           shape = guide_legend(override.aes = list(size = 4)))
}

sbs_colors = setNames(
  nm = c("SBS1",
         "SBS17b",
         "SBS18",
         "SBS5",
         "SBS88",
         "SBS10b",
         "SBS13",
         "SBS9",
         "SBS25",
         "SBS4",
         "SBS11",
         "SBS3" ,
         "SBS26",
         "SBS2",
         "SBS10a",
         "SBS31"), 
  object = c('#f1696bff', 
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
)
id_colors = setNames(
  nm = c('ID1', 
         'ID2', 
         "ID4",  
         "ID5",  
         "ID7",  
         "ID8",   
         "ID9",  
         "ID3",
         "ID18"), 
  object = c('#0c8281ff', 
             '#f5a55fff', 
             '#7d287eff', 
             '#2e4f4fff', 
             '#c4ddbcff', 
             '#996869ff', 
             '#daa627ff',
             'sienna3',
             '#bc8f8fff')
)
