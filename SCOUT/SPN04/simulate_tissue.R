rm(list = ls())
library(ProCESS)
library(dplyr)

dir <- getwd()
outdir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN04/process/"

setwd(outdir)
seed = 777
set.seed(seed)

# Prep simulation ####
sim <- TissueSimulation("SPN04", width=1e3, height=1e3, seed = seed, save_snapshots = TRUE)
sim$history_delta <- 1
sim$death_activation_level <- 50

# Clone 1 ####
sim$add_mutant(name = "Clone 1", growth_rates = .1, death_rates = .01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 1000)

# Clone 2 ####
sim$add_mutant("Clone 2",growth_rates = .3, death_rates = .01)
sim$update_rates("Clone 1", rates = c(growth = .033))
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size("Clone 2", 2000)
sim$update_rates("Clone 1", rates = c(growth = 0, death = .025))
sim$run_up_to_size("Clone 2", 4000)
sim$update_rates("Clone 1", rates = c(growth = 0, death = .05))
sim$run_up_to_size("Clone 2", 8000)

# Clone 3 ####
sim$add_mutant("Clone 3",growth_rates = 1, death_rates = .01)
sim$update_rates("Clone 2", rates = c(growth = .1))
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$run_up_to_size("Clone 3", 8000)
sim$update_rates("Clone 2", rates = c(growth = 0, death = .025))
sim$run_up_to_size("Clone 3", 16000)
sim$update_rates("Clone 2", rates = c(growth = 0, death = .05))
sim$run_up_to_size("Clone 3", 50000)
sim$update_rates("Clone 1", rates = c(growth = 0, death = 100))
sim$update_rates("Clone 2", rates = c(growth = 0, death = 100))
sim$run_up_to_size("Clone 3", 70000)

# Sample A ####
x_center = 470
y_center = 470
d = 23
bbox = list(lower_corner=c(x_center-d,y_center-d), upper_corner=c(x_center+d,y_center+d))
sim$sample_cells("SPN04_1.1", bbox$lower_corner, bbox$upper_corner)
t1 <- plot_tissue(sim, num_of_bins = 300)

# Treatment ####
treatment_start <- sim$get_clock()
sim$update_rates("Clone 3",rates = c(growth = .01, death=.3))
v <- sim$var("Clone 3")
#condition <- v <= 5000
condition <- v <= 50
sim$run_until(condition)
sim$update_rates("Clone 3",rates = c(growth = .2, death=.001))
sim$run_up_to_size("Clone 3", 1000)
treatment_end <- sim$get_clock()

# Relapse ####
sim$update_rates("Clone 3",rates = c(growth = .5, death=.001))
sim$run_up_to_size("Clone 3", 60000)

# Sample B ####
x_center = 400
y_center = 490
d = 25
bbox = list(lower_corner=c(x_center-d,y_center-d), upper_corner=c(x_center+d,y_center+d))
sim$sample_cells("SPN04_2.1", bbox$lower_corner, bbox$upper_corner)
t2 <- plot_tissue(sim, num_of_bins = 300)

# Save results ####
# Forest

forest <- sim$get_sample_forest()
plot_forest(forest) %>% 
  annotate_forest(forest, MRCAs = TRUE, samples = TRUE)
forest$save(paste0(outdir,"sample_forest.sff"))

# Treatment Info
treatment_info <- list(treatment_start=treatment_start, treatment_end=treatment_end)
saveRDS(treatment_info, "treatment_info.rds")
