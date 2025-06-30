rm(list = ls())
options(bitmapType='cairo')
library(tidyverse)
library(optparse)
library(caret)
library(dplyr)
library(patchwork)

source('../../getters/process_getters.R')
source('../../getters/tumourevo_getters.R')

option_list <- list( 
  make_option(c("-v", "--variantcaller"), type="character", default='mutect2', help="variant caller"),
  make_option(c("-c", "--cnacaller"), type="character", default='ascat', help="cna caller")
)
param <- parse_args(OptionParser(option_list=option_list))

SPN <- paste0('SPN0', seq(1,7))
COV <- c(50, 100, 150, 200)
PUR <- c(0.3, 0.6, 0.9)

final_table <- tibble()
for (spn in SPN){
  print(spn)
  for (cov in COV){
    for (pur in PUR){
      samples <- get_sample_names(spn)
      
      table <- lapply(samples, FUN = function(sample){
        file <- get_tumourevo_qc(spn = spn, coverage = cov, purity = pur, tool = 'tinc', vcf_caller = param$variantcaller, cna_caller = param$cnacaller, sample = sample)

        if (length(file) > 0){
          if (file.exists(file$fit_rds)){
            data <- readRDS(file$fit_rds)
            tmp <- tibble(TIT = data$TIT, TIN = data$TIN, sample = sample)
            return(tmp)
          }
        }
      }) %>% bind_rows()
      
      if (nrow(table) > 1){
        table$coverage = as.numeric(cov)
        table$purity = as.numeric(pur)
        table$vc = param$variantcaller
        table$cnc = 'ascat'
        table$SPN = spn
        table$N = length(samples)
        table <- table %>% mutate(error = abs(purity-TIT))
      }
      final_table <- bind_rows(final_table, table)
    }
  }
}
#saveRDS(object = table, file = paste0(out, param$cov,'x_',param$pur,'p_',param$variantcaller,'_ascat.rds'))

sp <- ggpubr::ggscatter(final_table, x = "TIT", y = "purity",
                color = "SPN", palette = "jco", position = 'jitter',
                add = "reg.line") 
corr <- sp + ggpubr::stat_cor(aes(color = SPN), label.x = 0.6, label.y.npc = c(0.27, 0.3))


v1 <- final_table %>% 
  ggplot() +
  geom_abline(linewidth = 0.5, col = 'gray') +
  geom_jitter(aes(x = TIT, y = purity, col = SPN), size = 3) +
  geom_point(aes(x = purity, y = purity), size = 3, shape = 8) +
  facet_grid(coverage ~ purity) + 
  xlab('TIT (TINC)') +
  ylab('purity (ProCESS)') + 
  scale_color_manual(values = c('steelblue', 'seagreen', 'goldenrod', 'coral', 'palevioletred', 'indianred3')) +
  xlim(0,1) +
  ylim(0,1) +
  theme_bw()


v2 <- ggplot(final_table, aes(x = SPN, y = error, col = SPN)) +
  geom_boxplot() +
  geom_jitter() +
  ylab('|ProCESS - TINC|') +
  scale_color_manual(values = c('steelblue', 'seagreen', 'goldenrod', 'coral', 'palevioletred', 'indianred3')) +
  facet_grid(coverage ~ purity) +
  theme_bw()
  

ggsave(filename = 'validation_TINC_v1.png', plot = v1, width = 7, height = 4, units = 'in', dpi = 600)
ggsave(filename = 'validation_TINC_v2.png', plot = v2, width = 9, height = 4, units = 'in', dpi = 600)
