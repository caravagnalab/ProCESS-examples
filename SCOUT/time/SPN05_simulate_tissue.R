### Setup ###
library(ProCESS)
library(dplyr)
#library(ggplot2)
#library(patchwork)

seed <- 9
set.seed(seed)
p_info <- ps::ps_handle()
start_time <- Sys.time()
initial_cpu <- ps::ps_cpu_times(p_info)
initial_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3

### Simulation ###
#outdir <- "/scratch/salvatore.milite/SNP05_building_scripts/SPN05/process/"
outdir <- '/orfeo/cephfs/scratch/area/lvaleriani/races/ProCESS-examples/SCOUT/time/'


sim <- TissueSimulation(seed = seed, save_snapshots = F, name = "SPN05")
sim$history_delta <- 1

# Set the death activation level to avoid drift
sim$death_activation_level <- 50
# Set initial growth rates
Clone1_gr <- 0.1
Clone2_gr <- 0.3
Clone3_gr <- 0.6
Clone4_gr <- 0.8
Clone5_gr <- 1.1

# Add first mutant
sim$add_mutant(name = "Clone 1",
               growth_rates = Clone1_gr,
               death_rates = 0.01)

sim$place_cell("Clone 1", 500, 500)

# Let the simulation evolve until "Clone 1" consists of 1300 cells
sim$run_up_to_size("Clone 1", 100)
clone1_tissue <-plot_tissue(sim)
clone1_muller <- plot_muller(sim)

# Add second mutant
sim$add_mutant("Clone 2",growth_rates = Clone2_gr,
               death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")


# Let the simulation evolve until "Clone 2" consists of 5000 cells
sim$run_up_to_size("Clone 2", 100)
sim$update_rates("Clone 1",rates = c(growth = 0, death=0.05))
sim$run_up_to_size("Clone 2", 200)
sim$update_rates("Clone 2",rates = c(growth = Clone2_gr, death=0.05))
sim$run_up_to_size("Clone 2", 500)
sim$update_rates("Clone 1",rates = c(growth = 0, death=0.1))

clone2_tissue <- plot_tissue(sim)
clone2_muller <- plot_muller(sim)


# Add third mutant
sim$add_mutant("Clone 3",growth_rates = Clone3_gr,
               death_rates = 0.01)
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
clone3_clock = sim$get_clock()
clone3_clock %>% saveRDS(file.path(outdir, "clone3_clock.rds"))
sim$update_rates("Clone 2",rates = c(growth = 0, death=0.05))
sim$run_up_to_size("Clone 3", 100)
sim$update_rates("Clone 2",rates = c(growth = 0, death=0.1))
sim$run_up_to_size("Clone 3", 1e2)

clone3_tissue <- plot_tissue(sim)
clone3_muller <- plot_muller(sim)

sim$add_mutant("Clone 4",growth_rates = Clone4_gr,
               death_rates = 1e-18)

sim$mutate_progeny(sim$choose_cell_in("Clone 3"), "Clone 4")

sim$run_up_to_size("Clone 4", 800)


sim$add_mutant("Clone 5",growth_rates = Clone5_gr,
               death_rates = 1e-18)
sim$update_rates("Clone 3",rates = c(growth = Clone3_gr - 0.2, death=0.1))

sim$update_rates("Clone 2",rates = c(growth = 0, death=0.5))
sim$update_rates("Clone 1",rates = c(growth = 0, death=0.5))

sim$mutate_progeny(sim$choose_cell_in("Clone 3"), "Clone 5")
sim$run_up_to_size("Clone 5",1e4)



#final_muller <- plot_muller(sim)
#final_tissue <- plot_tissue(sim)
#piechart <- plot_state(sim)
#timeseries_plot <- plot_timeseries(sim)

#tissue_plots <- list(muller = final_muller, 
#                     issue = final_tissue, 
#                     piechart = piechart, 
#                     timeseries = timeseries_plot)

#saveRDS(tissue_plots, "tissue_plots.rds")



#final_tissue + geom_rect(xmin = 505 - 20, xmax = 525 + 30,
#                         ymin = 410 - 15, ymax = 430 + 25,
#                         fill = NA, color = "black")

### Sampleing ###

sim$sample_cells("SPN05_1.1",c(505 - 20,410 - 15),c(525+ 30,430+ 25), num_of_cells= 2000)

#final_tissue + geom_rect(xmin = 420, xmax = 480,
#                         ymin = 350, ymax = 400,
#                         fill = NA, color = "black")

sim$sample_cells("SPN05_1.2", c(420,350), c(480,400), num_of_cells= 2000)

#tissue_sampled_B <- plot_tissue(sim)


#final_tissue + geom_rect(xmin = 405 - 15, xmax = 425 + 65,
#                         ymin = 535, ymax = 550 + 10,
#                         fill = NA, color = "black")

sim$sample_cells("SPN05_1.3", c(405 - 15 ,535), c(425 + 65,550 + 10),  num_of_cells= 2000)

tissue_sampled_C <- plot_tissue(sim)

### Forest ###
forest <- sim$get_sample_forest()
forest$save(paste0(outdir,"sample_forest.sff"))

end_time <- Sys.time()
final_cpu <- ps::ps_cpu_times(p_info)
final_mem <- ps::ps_memory_info(p_info)["rss"] / 1024^3
elapsed_time <- end_time - start_time
elapsed_time <- as.numeric(elapsed_time, units = "mins")
cpu_used <- (final_cpu["user"] + final_cpu["system"]) - (initial_cpu["user"]+ initial_cpu["system"])
mem_used <- final_mem - initial_mem
resource_usage <- data.frame(
  elapsed_time_mins =  elapsed_time,
  cpu_time_secs = cpu_used,
  memory_used_MB = mem_used
)
saveRDS(object = resource_usage, file = 'SPN05_tissue.rds')


