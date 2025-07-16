rm(list = ls())
library(ProCESS)
library(dplyr)


dir <- getwd()
outdir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN01/process/"
setwd(outdir)
set.seed(06117)
sim <- TissueSimulation(name = 'SPN01', seed = 12345, save_snapshot=TRUE)
sim$history_delta <- 1
sim$death_activation_level <- 50



# Clone 1
sim$add_mutant(name = "Clone 1", growth_rates = 0.08, death_rates = 0.01)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 500)

# Clone 2
sim$add_mutant(name = "Clone 2", growth_rates = 0.3, death_rates = 0.01)
sim$update_rates(species = "Clone 1", rates = c(growth = 0.06, death = 0.01))
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size(species = 'Clone 2', 700)
sim$update_rates(species = "Clone 1", rates = c(growth = 0.02, death = 0.03))
sim$update_rates(species = "Clone 2", rates = c(growth = 0.4, death = 0.005))
sim$run_up_to_size(species = 'Clone 2', 1000)


# Grow 2
sim$update_rates(species = "Clone 1", rates = c(growth = 0.001, death = 0.3))
sim$run_up_to_size(species = 'Clone 2', 2000)


# Clone 3
sim$add_mutant(name = "Clone 3", growth_rates = 1, death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.1, death = 0.01))
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$run_up_to_size("Clone 3", 4000)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.05, death = 0.5))
sim$run_up_to_size("Clone 3", 6000)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.01, death = 0.8))
sim$run_up_to_size("Clone 3", 8000) ## orignal 15000

sim$update_rates(species = "Clone 2", rates = c(growth = 0.01, death = 1))
sim$run_up_to_size("Clone 3", 10000)

# Clone 4
sim$add_mutant(name = "Clone 4", growth_rates = 2.5 , death_rates = 0.01)
sim$update_rates(species = "Clone 2", rates = c(growth = 0.005, death = 1))
sim$update_rates(species = "Clone 3", rates = c(growth = 0.5, death = 0.02))
sim$mutate_progeny(sim$choose_cell_in("Clone 3"), "Clone 4")
sim$run_up_to_size("Clone 4", 4000) ## original 6000
sim$update_rates(species = "Clone 3", rates = c(growth = 0.02, death = 0.02))
sim$run_up_to_size("Clone 4", 6000) ##60000 ## orignal 8000
#### try this
sim$update_rates(species = "Clone 3", rates = c(growth = 0.01, death = 0.01))
sim$run_up_to_size("Clone 4", 8000)
sim$run_up_to_size("Clone 4", 15000)

print("End simulation")
## Sample 3
bboxC_lower_corner <- c(300, 400)
bboxC_upper_corner <- c(350, 450)
sim$sample_cells("SPN01_1.3", bboxC_lower_corner, bboxC_upper_corner)

## Sample 1
bboxA <- sim$search_sample(c("Clone 4" = 400, 'Clone 3' = 400), 50,50)
sim$sample_cells("SPN01_1.1", bboxA$lower_corner, bboxA$upper_corner)

bboxB_lower_corner <- c(420,500)
bboxB_upper_corner <- c(470,550)
sim$sample_cells("SPN01_1.2", bboxB_lower_corner, bboxB_upper_corner)

print("SPN01_1.2 has been sampled")

# # Forest
forest <- sim$get_sample_forest()
forest$save(paste0(outdir,"sample_forest.sff"))
