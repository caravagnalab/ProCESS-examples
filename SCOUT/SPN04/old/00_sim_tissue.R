rm(list = ls())
library(ProCESS)
library(dplyr)
library(ggplot2)
library(patchwork)
source("utils.R")

seed <- 777
set.seed(seed)

# Prep simulation ####
sim <- SpatialSimulation("SPN01", width=1e3, height=1e3, seed = seed)
sim$history_delta <- 1
sim$death_activation_level <- 50
#sim <- new(Simulation, seed = seed, save_snapshot = F)
#sim$duplicate_internal_cells <- T

# First and Second mutant ####
sim$add_mutant(name = "Clone 1", growth_rates = .1, death_rates = .01)

####### CLONE 1 #######

sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 1000)
#plot_muller(sim)

####### CLONE 2 #######

sim$add_mutant("Clone 2",growth_rates = .3, death_rates = .01)
sim$update_rates("Clone 1", rates = c(growth = .033))

sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size("Clone 2", 2000)

sim$update_rates("Clone 1", rates = c(growth = 0, death = .025))
sim$run_up_to_size("Clone 2", 4000)

sim$update_rates("Clone 1", rates = c(growth = 0, death = .05))
sim$run_up_to_size("Clone 2", 8000)
#plot_muller(sim)

####### CLONE 3 #######

sim$add_mutant("Clone 3",growth_rates = 1, death_rates = .01)
#sim$update_rates("Clone 2", rates = c(growth = .1, death = .05))
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

#plot_tissue(sim)
#ggsave("tissue/tissue_00.pdf", dpi=300, width = 8, height = 8)
#plot_muller(sim)
#ggsave("tissue/muller_00.pdf", dpi=300, width = 8, height = 8)


####### SAMPLE A #######
n_w <- n_h <- 20
ncells <- .99*n_w*n_h
bbox <- sim$search_sample(c("Clone 3" = ncells), n_w, n_h)
sim$sample_cells("Sample A", bbox$lower_corner, bbox$upper_corner)

t1 <- plot_tissue(sim, num_of_bins = 300)
#ggsave("tissue/tissue_01.pdf", dpi=300, width = 8, height = 8)
#plot_muller(sim)
#ggsave("tissue/muller_01.pdf", dpi=300, width = 8, height = 8)

####### TREATMENT #######
treatment_start <- sim$get_clock()

sim$update_rates("Clone 3",rates = c(growth = 0, death=.3))
v <- sim$var("Clone 3")
condition <- v <= 1000
sim$run_until(condition)

sim$update_rates("Clone 3",rates = c(growth = .1, death=.08))
sim$run_up_to_size("Clone 3", 1100)

treatment_end <- sim$get_clock()

#plot_tissue(sim, num_of_bins = 300)
#ggsave("tissue/tissue_02.pdf", dpi=300, width = 8, height = 8)
#plot_muller(sim) + xlim(20, NA)
#ggsave("tissue/muller_02.pdf", dpi=300, width = 8, height = 8)

####### RELAPSE #######
sim$update_rates("Clone 3",rates = c(growth = .5, death=.001))
sim$run_up_to_size("Clone 3", 60000)

#plot_tissue(sim, num_of_bins = 300)
#ggsave("tissue/tissue_03.pdf", dpi=300, width = 8, height = 8)
#plot_muller(sim) + xlim(20, NA)
#ggsave("tissue/muller_03.pdf", dpi=300, width = 8, height = 8)

####### SAMPLE B #######
n_w <- n_h <- 20
ncells <- .99*n_w*n_h
bbox <- sim$search_sample(c("Clone 3" = ncells), n_w, n_h)
sim$sample_cells("Sample B", bbox$lower_corner, bbox$upper_corner)

t2 <- plot_tissue(sim, num_of_bins = 300)
#ggsave("tissue/tissue_04.pdf", dpi=300, width = 8, height = 8)
#plot_muller(sim) + xlim(20, NA) + geom_vline(xintercept = c(treatment_start, treatment_end), color="indianred", linewidth=.3)
#ggsave("tissue/muller_04.pdf", dpi=300, width = 8, height = 8)

# Save
muller <- plot_muller(sim)
forest <- sim$get_samples_forest()
forest$save("data/samples_forest.sff")

treatment_info <- list(treatment_start=treatment_start, treatment_end=treatment_end)
saveRDS(treatment_info, "data/treatment_info.rds")

plt_forest <- plot_forest(forest) %>%
  annotate_forest(forest)

piechart <- plot_state(sim)
timeseries <- plot_timeseries(sim)

pl <- t1 + t2 + piechart + timeseries + muller + plt_forest + plot_layout(design = 'AABB\nEFGG\nHHHH\nHHHH\nHHHH')
ggsave('plots/SPN04_tissue.png', plot = pl, width = 210, height = 297, units = "mm", dpi = 300)
ggsave('plots/SPN04_tissue.pdf', plot = pl, width = 210, height = 297, units = "mm", dpi = 300)

