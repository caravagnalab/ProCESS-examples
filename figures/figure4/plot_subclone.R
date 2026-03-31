library(dplyr)
library(ggplot2)
library(ggnewscale)
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/figures/figure3/utils_plot.R")
source("/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/validation/SCOUT/colors.R")

fga_df <-  tibble(spn = paste0('SPN0',1:7),
                         class = c('High FGA', 'Low FGA', 'Low FGA', 'Low FGA', 'High FGA', 'Low FGA', 'High FGA'),
                  status = c('Polyclonal', 'Monoclonal', 'Polyclonal', 'Monoclonal', 'Polyclonal', 'Polyclonal', 'Polyclonal'))

data <- readRDS('/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/figures/figure4/performance_subclones.rds') %>% 
  left_join(fga_df)


stat.test <- data %>%
  mutate(purity = factor(purity),
         coverage = factor(coverage)) %>% 
  rstatix::wilcox_test(Accuracy ~ status) %>%
  rstatix::adjust_pvalue(method = "BH") %>%       
  rstatix::add_significance("p.adj") %>%          
  rstatix::add_xy_position(x = "class") %>% 
  mutate(y.position = 1.05)



data %>% 
  ggplot(aes(x = status, y = Accuracy)) +
  
  # FIRST COLOR SCALE (class)
  geom_boxplot(
    aes(color = as.factor(purity), fill = as.factor(purity)),
    outliers = TRUE, alpha = .4
  )  +
  my_ggplot_theme() + 
  scale_fill_manual('Purity', values = c('0.3' = '#A8CCDA', '0.6' = '#759DB9', '0.9' ='#426D98' ))+
  scale_color_manual('Purity',values = c('0.3' = '#A8CCDA', '0.6' = '#759DB9', '0.9' ='#426D98' ))+
  xlab('Clonal Status') +
  ylab('Subclonal Detection\nAccuracy')

data %>% 
  ggplot(aes(x = status, y = Accuracy)) +
  
  # FIRST COLOR SCALE (class)
  geom_boxplot(
    #aes(color = as.factor(purity), fill = as.factor(purity)),
    outliers = TRUE, alpha = .4
  )  +
  my_ggplot_theme() + 
  xlab('Clonal Status') +
  ylab('Subclonal Detection\nAccuracy')
  

plt_mobster <- data %>% 
  ggplot(aes(x = class, y = Accuracy)) +
  
  # FIRST COLOR SCALE (class)
  geom_boxplot(
    aes(color = class),
    outliers = TRUE, alpha = .4
  ) +
  
  scale_color_manual('FGA Class', values = c(
    "High FGA" = "indianred2",
    "Low FGA"  = "dodgerblue3"
  )) +
  
  # 🔴 RESET COLOR SCALE HERE
  ggnewscale::new_scale_color() +
  
  # SECOND COLOR SCALE (spn)
  stat_summary(
    aes(color = spn),
    fun.data = mean_cl_boot,
    position = position_dodge(width = 0.3),
    size = .4, show.legend = F
  ) +
  
  stat_summary(
    aes(color = spn, group = spn),
    fun.data = mean_cl_boot,
    geom = 'line',
    position = position_dodge(width = 0.3),
    linewidth = .35, show.legend = F
  ) +
  
  scale_color_manual(values = SPN_colors, name = "SPN") +
  
  stat_pvalue_manual(
    stat.test,
    label = "p.adj.signif",
    tip.length = 0.01,
    inherit.aes = FALSE,
    color = "gray40"
  ) +
  
  my_ggplot_theme() +
  xlab('FGA class') +
  ylab('Subclonal Detection\nAccuracy')
plt_mobster

