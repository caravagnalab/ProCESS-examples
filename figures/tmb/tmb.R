library(maftools)
library(tidyverse)

setwd('/orfeo/cephfs/scratch/area/lvaleriani/races/maf')
source('../ProCESS-examples/getters/process_getters.R')
source("/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/report/plotting/utils.R")

spns = paste0('SPN0', 1:7)

mafs <- lapply(spns, FUN = function(spn){
    print(spn)
    samples_name = get_sample_names(spn)
    samples = paste0(spn, '_', samples_name)

    mafs = lapply(samples, FUN = function(s){
      print(s)
      data = read.maf(maf = paste0(s, '.maf'))
      filter = subsetMaf(data, query = 't_alt_count > 0')
      return(filter)
    })

    names(mafs) <- samples
    return(mafs %>% merge_mafs())
})  

names(mafs) = spns


maf = mafs #merge_mafs
capture_size = NULL
tcga_capture_size = 35.8
cohortName = names(mafs)
tcga_cohorts = NULL
primarySite = FALSE
col = c('gray70', 'black')
bg_col = c('peachpuff', 'lightblue3')
medianCol = 'red'
decreasing = FALSE
logscale = TRUE
rm_hyper = FALSE
rm_zero = TRUE
cohortFontSize = 0.8
axisFontSize = 0.8

tcga.cohort = system.file('extdata', 'tcga_cohort.txt.gz', package = 'maftools')
tcga.cohort = data.table::fread(file = tcga.cohort, sep = '\t', stringsAsFactors = FALSE)

if(primarySite){
  tcga.cohort = tcga.cohort[,.(Tumor_Sample_Barcode, total, site)]
  colnames(tcga.cohort)[3] = 'cohort'
}else{
  tcga.cohort = tcga.cohort[,.(Tumor_Sample_Barcode, total, cohort)]
}

if(!is.null(tcga_cohorts)){
  tcga.cohort = tcga.cohort[cohort %in% tcga_cohorts]
  if(nrow(tcga.cohort) == 0){
    stop("Something went wrong. Provide correct names for 'tcga_cohorts' arguments")
  }
}

if(length(maf) == 1){
  maf = list(maf)
}

maf.mutload = lapply(maf, function(m){
  x = maftools::getSampleSummary(m)[,.(Tumor_Sample_Barcode, total)]
  print(x)
  if(rm_zero){
    warning(paste0("Removed ", nrow(x[x$total == 0]), " samples with zero mutations."))
    x = x[!total == 0]
  }
  x
})

if(is.null(cohortName)){
  cohortName = paste0('Input', seq_len(length(maf)))
}else if(length(cohortName) != length(maf)){
  stop("Please provide names for all input cohorts")
}

names(maf.mutload) = cohortName
maf.mutload = data.table::rbindlist(l = maf.mutload, idcol = "cohort")

#maf.mutload[,cohort := cohortName]
tcga.cohort$total = as.numeric(as.character(tcga.cohort$total))
maf.mutload$total = as.numeric(as.character(maf.mutload$total))

if(rm_hyper){
  #Remove outliers from tcga cohort
  tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
  tcga.cohort = lapply(tcga.cohort, function(x){
    xout = boxplot.stats(x = x$total)$out
    if(length(xout) > 0){
      message(paste0("Removed ", length(xout), " outliers from ", x[1,cohort]))
      x = x[!total %in% xout]
    }
    x
  })
  tcga.cohort = data.table::rbindlist(l = tcga.cohort)
  #Remove outliers from input cohort
  xout = boxplot.stats(x = maf.mutload$total)$out
  if(length(xout) > 0){
    message(paste0("Removed ", length(xout), " outliers from Input MAF"))
    maf.mutload = maf.mutload[!total %in% xout]
  }
}

if(is.null(capture_size)){
  tcga.cohort = rbind(tcga.cohort, maf.mutload)
  pt.test = pairwise.t.test(x = tcga.cohort$total, g = tcga.cohort$cohort, p.adjust.method = "fdr")
  pt.test.pval = as.data.frame(pt.test$p.value)
  data.table::setDT(x = pt.test.pval, keep.rownames = TRUE)
  colnames(pt.test.pval)[1] = 'Cohort'
  message("Performing pairwise t-test for differences in mutation burden..")
  pt.test.pval = data.table::melt(pt.test.pval, id.vars = "Cohort")
  colnames(pt.test.pval) = c("Cohort1", "Cohort2", "Pval")
  tcga.cohort[, plot_total := total]
}else{
  names(capture_size) = cohortName
  message(paste0("Capture size [TCGA]:  ", tcga_capture_size))
  message(paste0("Capture size [Input]: ", paste(capture_size, collapse = ", ")))
  #maf.mutload[,total_perMB := total/capture_size]
  maf.mutload.spl = split(maf.mutload, maf.mutload$cohort)
  maf.mutload.spl = lapply(names(maf.mutload.spl), function(m){
    x = maf.mutload.spl[[m]]
    x[, total_perMB := total / capture_size[m]]
    x
  })
  maf.mutload = data.table::rbindlist(l = maf.mutload.spl, use.names = TRUE, fill = TRUE)
  tcga.cohort[,total_perMB := total/tcga_capture_size]
  tcga.cohort = rbind(tcga.cohort, maf.mutload)
  message("Performing pairwise t-test for differences in mutation burden (per MB)..")
  pt.test = pairwise.t.test(x = tcga.cohort$total_perMB, g = tcga.cohort$cohort, p.adjust.method = "fdr")
  pt.test.pval = as.data.frame(pt.test$p.value)
  data.table::setDT(x = pt.test.pval, keep.rownames = TRUE)
  colnames(pt.test.pval)[1] = 'Cohort'
  pt.test.pval = data.table::melt(pt.test.pval, id.vars = "Cohort")
  colnames(pt.test.pval) = c("Cohort1", "Cohort2", "Pval")
  tcga.cohort[, plot_total := total_perMB]
}

#Median mutations
tcga.cohort.med = tcga.cohort[,.(.N, median(plot_total)),cohort][order(V2, decreasing = decreasing)]
tcga.cohort$cohort = factor(x = tcga.cohort$cohort,levels = tcga.cohort.med$cohort)
colnames(tcga.cohort.med) = c('Cohort', 'Cohort_Size', 'Median_Mutations')
tcga.cohort$TCGA = ifelse(test = tcga.cohort$cohort %in% cohortName,
                          yes = 'Input', no = 'TCGA')

#Plotting data
tcga.cohort = split(tcga.cohort, as.factor(tcga.cohort$cohort))
plot.dat = lapply(seq_len(length(tcga.cohort)), function(i){
  x = tcga.cohort[[i]]
  x = data.table::data.table(rev(seq(i-1, i, length.out = nrow(x))),
                             x[order(plot_total, decreasing = TRUE), plot_total],
                             x[,TCGA])
  x
})
names(plot.dat) = names(tcga.cohort)


data = lapply(names(plot.dat), FUN = function(i){
  if (i %in% c('SPN03', 'SPN04')){
    plot.dat[[i]] %>% 
      as_tibble() %>%
      mutate(cohort = 'LAML') %>% 
      mutate(spn = i) %>% 
      mutate(type = 'Normal')
  } else if (i %in% c('SPN01', 'SPN02')){
    if (i == 'SPN01'){
      plot.dat[[i]] %>% 
        as_tibble() %>%
        mutate(cohort = 'COAD') %>% 
        mutate(spn = i) %>% 
        mutate(type = c('Normal', 'WGD','WGD'))
    } else {
      plot.dat[[i]] %>% 
        as_tibble() %>%
        mutate(cohort = 'COAD') %>% 
        mutate(spn = i) %>% 
        mutate(type = c('Hypermutant','Hypermutant'))
    }
  } else if (i %in% c('SPN05')){
    plot.dat[[i]] %>% 
      as_tibble() %>%
      mutate(cohort = 'BRCA') %>% 
      mutate(spn = i) %>% 
      mutate(type = 'Hypermutant')
  } else if (i %in% c('SPN06')){
    plot.dat[[i]] %>% 
      as_tibble() %>%
      mutate(cohort = 'LUAD') %>% 
      mutate(spn = i) %>% 
      mutate(type = c('Normal', 'Normal', 'Normal','WGD','WGD'))
  } else if (i %in% c('SPN07')){
    plot.dat[[i]] %>% 
      as_tibble() %>%
      mutate(cohort = 'GBM') %>% 
      mutate(spn = i)  %>% 
      mutate(type = c('Normal', 'Normal', 'Normal', 'Hypermutant','Hypermutant'))
  } else {
    plot.dat[[i]] %>% 
      as_tibble() %>%
      mutate(cohort = i) %>% 
      mutate(spn = V3) %>% 
      mutate(type = 'Normal')
  }
}) %>% 
  bind_rows() %>% 
  filter(cohort %in% c('LAML', 'BRCA', 'LUAD', 'GBM', 'COAD'))

data$cohort <- factor(
  data$cohort,
  levels = c("LAML", "BRCA", "GBM", "COAD", "LUAD")
)

data$facet_id <- as.numeric(factor(data$cohort))

data <- data %>%
  group_by(cohort) %>%
  mutate(
    V1_scaled = if_else(
      spn != "TCGA",                        # only scale SPN points
      scales::rescale(V1, to = range(V1[spn == "TCGA"])),  # rescale to TCGA range
      V1                                  # TCGA keeps original
    )
  ) %>%
  ungroup()

bg <- data %>%
  dplyr::distinct(cohort, facet_id) %>%
  dplyr::mutate(col = ifelse(facet_id %% 2 == 0, "1", "2"))

med_df <- data %>%
  dplyr::group_by(cohort) %>%
  dplyr::summarise(med = median(log10(V2), na.rm = TRUE))

col_point <- c("SPN01"='steelblue', 
               "SPN02"='seagreen', 
               "SPN03"='goldenrod', 
               "SPN04"='coral', 
               "SPN05"="magenta4",
               "SPN06"='palevioletred', 
               "SPN07"='indianred3', 
               'TCGA' = 'gray70')

cex_opt = getOption('CNAqc_cex', default = 1)
plt_tmb <- data %>% 
  ggplot() +
  geom_rect(
    data = bg,
    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, fill = col),
    inherit.aes = FALSE,
    alpha = 0.1,
    show.legend = F
  ) +
  scale_fill_manual(values = bg_col) +
  geom_hline(
    data = med_df,
    aes(yintercept = med),
    color = "gray60",
    linewidth = 0.7
  ) + 
  geom_hline(
    yintercept = 0:6,
    linetype = "dashed",
    color = "grey",
    linewidth = 0.3
  ) + 
  geom_point(data = data  %>% filter(spn == 'TCGA'), aes(x = V1_scaled, y = log10(V2), col = spn, shape = type), size =.4) +
  geom_point(data = data  %>% filter(spn != 'TCGA'), aes(x = V1_scaled,y = log10(V2), col = spn, shape = type), size = 2) +
  scale_color_manual('', values = col_point)+ 
  scale_shape_manual('', values = c('Normal' = 16, 'WGD' =15, 'Hypermutant' = 17))+ 
  facet_grid(.~cohort, scales = 'free_x', switch = 'x') +
  scale_y_continuous(
    breaks = 0:4.5,
    labels = function(x) 10^x, limits = c(0,4.5)
  ) + 
  ylab('TMB') +
  xlab('') +
  ggplot2::theme_light(base_size = 10 * cex_opt) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.text.y = element_text(size = 8, margin = margin(l = 0, r = 0), colour = 'gray20'),
    panel.background = element_blank(),      # removes grey background
    panel.grid.major = element_blank(),      # removes major grid lines
    panel.grid.minor = element_blank(),      # removes minor grid lines
    panel.spacing = unit(0, "mm"),                       # remove spacing between facets
    panel.border = element_blank(),
    plot.background = element_rect(color = "black", fill = NA),
    strip.background = element_blank(),
    strip.text.x = element_text(
      angle = 0,      # rotation
      vjust = 0.5,     # vertical alignment
      hjust = 0.5,      # horizontal alignment,
      size = 8, 
      margin = margin(t = 0.1, b = 0.1), 
      colour = 'gray20'
    )
  ) +
  guides(
    color = guide_legend(override.aes = list(size = 2)),
    fill  = guide_legend(override.aes = list(size = 2)),
    shape = guide_legend(override.aes = list(size = 2))
  )  
  #my_ggplot_theme()
plt_tmb

# ggsave(plot = plt, 
#        filename = '/orfeo/cephfs/scratch/area/lvaleriani/races/maf/tmb.pdf', 
#        width = 4, height = 6, units = 'in')
