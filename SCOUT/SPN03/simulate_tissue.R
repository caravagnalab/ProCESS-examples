library(ProCESS)
library(dplyr)
library(ggplot2)

rm(list = ls())
base <- '/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN03/process/'
setwd(base)
unlink("SPN03", recursive = TRUE)

set.seed(12345)
sim <- TissueSimulation(name = 'SPN03', seed = 12345, save_snapshots = T)
sim$history_delta <- 1 
sim$death_activation_level <- 50

# Clone 1
sim$add_mutant(name = "Clone 1", growth_rates = 0.3, death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 10000)


# Sample 1.1
bbox <- sim$search_sample(c("Clone 1" = 2000), 50, 50)
sim$sample_cells("SPN03_1.1", bbox$lower_corner, bbox$upper_corner)
t1 <- plot_tissue(sim) + geom_rect(xmin = bbox$lower_corner[1], xmax = bbox$upper_corner[1],
                                   ymin = bbox$lower_corner[2], ymax = bbox$upper_corner[2],
                                   fill = NA, color = "black")

# Clone 2
sim$update_rates(species = "Clone 1", rates = c(growth = 0.1, death = 0.01))
sim$add_mutant(name = "Clone 2", growth_rates = 0.5, death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size(species = 'Clone 2', 20000)


# Sample 2.1
bbox_lower_corner <- c(451,426)
bbox_upper_corner <- c(500,475)
sim$sample_cells("SPN03_2.1", bbox_lower_corner, bbox_upper_corner)
t2 <- plot_tissue(sim) + geom_rect(xmin = bbox_lower_corner[1], xmax = bbox_upper_corner[1],
                             ymin = bbox_lower_corner[2], ymax = bbox_upper_corner[2],
                             fill = NA, color = "black")


# Update 
sim$update_rates(species = "Clone 1", rates = c(growth = 0.05, death = 0.05))
sim$update_rates(species = "Clone 2", rates = c(growth = 1, death = 0.01))
sim$run_up_to_size(species = 'Clone 2', 35000)

sim$update_rates(species = "Clone 1", rates = c(growth = 0.01, death = 0.07))
sim$update_rates(species = "Clone 2", rates = c(growth = 1.5, death = 0.005))
sim$run_up_to_size(species = 'Clone 2', 55000)

sim$update_rates(species = "Clone 1", rates = c(growth = 0.01, death = 0.07))
sim$update_rates(species = "Clone 2", rates = c(growth = 1.8, death = 0.005))
sim$run_up_to_size(species = 'Clone 2', 70000)

# Grow 2
sim$update_rates(species = "Clone 1", rates = c(growth = 0, death = 1))
sim$update_rates(species = "Clone 2", rates = c(growth = 2.1, death = 0.003))
sim$run_up_to_size(species = 'Clone 2', 85000)

# Grow 3
sim$update_rates(species = "Clone 2", rates = c(growth = 2.3, death = 0.005))
sim$run_up_to_size(species = 'Clone 2', 130000)

# Sample 3.1
bbox_lower_corner <- c(500,210)
bbox_upper_corner <- c(550,260)
sim$sample_cells("SPN03_3.1", bbox_lower_corner, bbox_upper_corner)
t3 <- plot_tissue(sim) + geom_rect(xmin = bbox_lower_corner[1], xmax = bbox_upper_corner[1],
                             ymin = bbox_lower_corner[2], ymax = bbox_upper_corner[2],
                             fill = NA, color = "black")



# Clone 3
sim$add_mutant(name = "Clone 3", growth_rates = 6, death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.3, death = 0.1))
cell <- sim$get_cells() %>% 
  filter(mutant == 'Clone 2', position_x > 550,position_x < 600, position_y<300, position_y>260) %>% 
  sample_n(1)
sim$mutate_progeny(cell, "Clone 3")
sim$run_up_to_size("Clone 3", 60000)


# Sample 4.1
bbox_lower_corner <- c(510,130)
bbox_upper_corner <- c(560,180)
sim$sample_cells("SPN03_4.1", bbox_lower_corner, bbox_upper_corner)
t4 <- plot_tissue(sim) + geom_rect(xmin = bbox_lower_corner[1], xmax = bbox_upper_corner[1],
                                   ymin = bbox_lower_corner[2], ymax = bbox_upper_corner[2],
                                   fill = NA, color = "black")

forest <- sim$get_sample_forest()
plt_forest <- plot_forest(forest) %>%
   annotate_forest(forest, MRCAs = T, samples = T)

forest$save("sample_forest.sff")
