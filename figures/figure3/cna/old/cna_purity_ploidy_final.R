library(ggnewscale)
library(ggplot2)
library(dplyr)
library(tidyr)

df_long <- df_all_combs_SPN_cna %>%
  mutate(correct_purity_estimation_class = case_when(
    abs(delta_purity) < 0.2 ~ "correctly estimated",
    delta_purity > 0.2 ~ "underestimated",
    delta_purity < -0.2 ~ "overestimated"),
    correct_ploidy_estimation_class = case_when(
      abs(delta_ploidy) < 1 ~ "correctly estimated",
      delta_ploidy > 1 ~ "underestimated",
      delta_ploidy < -1 ~ "overestimated"),
    true_purity_label = paste0('Purity =', true_purity)
  ) %>%
  pivot_longer(
    cols = c(delta_purity, delta_ploidy),
    names_to = "delta_type",
    values_to = "delta_value"
  ) %>%
  mutate(
    correct_estimation_class = ifelse(delta_type == "delta_purity",
                                      correct_purity_estimation_class,
                                      correct_ploidy_estimation_class)
  ) %>% 
  mutate(delta_type=case_when(delta_type=="delta_ploidy" ~ "Ploidy",
                              delta_type=="delta_purity" ~ "Purity"))


plt_final <- ggplot(df_long, aes(
  x = reorder(tool, delta_value, \(x) sd(x, na.rm = TRUE)),
  y = delta_value,
  group = interaction(tool, true_purity_label)
)) +
  
  # Boxplots colored by tool
  geom_boxplot(aes(color = tool, fill = tool), alpha = 0.3, outlier.shape = NA) +
  scale_color_manual(name = "CNA caller", values = col_cna_tools) +
  scale_fill_manual(name = "CNA caller", values = col_cna_tools) +
  
  ggnewscale::new_scale_color() +
  ggnewscale::new_scale_fill() +
  
  # Jitter points colored by FGA class
  geom_jitter(aes(color = fga_class, fill = fga_class, shape = correct_estimation_class),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
              size = 2, alpha = 0.8) +
  scale_color_manual(name = "FGA class", values = c("High FGA"="indianred3","Low FGA"="dodgerblue3")) +
  scale_fill_manual(name = "FGA class", values = c("High FGA"="indianred3","Low FGA"="dodgerblue3")) +
  scale_shape_manual(name = "Estimation class",
                     values = c("correctly estimated"=16,
                                "underestimated"=25,
                                "overestimated"=24)) +
  
  # Horizontal lines
  geom_hline(aes(yintercept = 0), linetype = "solid", color = "grey40") +
  geom_hline(data = subset(df_long, delta_type=="Purity"),
             aes(yintercept = 0.2), linetype = "dashed", color = "grey") +
  geom_hline(data = subset(df_long, delta_type=="Purity"),
             aes(yintercept = -0.2), linetype = "dashed", color = "grey") +
  geom_hline(data = subset(df_long, delta_type=="Ploidy"),
             aes(yintercept = 1), linetype = "dashed", color = "grey") +
  geom_hline(data = subset(df_long, delta_type=="Ploidy"),
             aes(yintercept = -1), linetype = "dashed", color = "grey") +
  
  # Facet grid: delta_type as rows, true purity as columns
  facet_grid(delta_type ~ true_purity_label, scales = "free_y") +
  
  labs(y = "Deviation", x = "") +
  theme_bw() +
  theme(axis.text.x = element_blank())
plt_final
