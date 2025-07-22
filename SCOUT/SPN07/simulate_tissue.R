#########################
######### SPN07 ######### 
#########################
rm(list = ls())
library(ProCESS)
library(dplyr)
library(ggplot2)

dir <- getwd()
outdir <- "/orfeo/scratch/cdslab/shared/SCOUT/SPN07/process" #"/orfeo/scratch/cdslab/antonelloa/ProCESS-examples/SCOUT/SPN07/"
setwd(outdir)
set.seed(12345)
sim <- TissueSimulation(name = 'SPN07', seed = 12345, save_snapshot=T, width = 1e3, height = 1e3) #save_snapshot=T
sim$history_delta <- .1
sim$death_activation_level <- 50
sim$border_growth_model <- F

# Clone 1
sim$add_mutant(name = "Clone 1", growth_rates = .1, death_rates = .03)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1",10) #1e2
# plot_tissue(sim)

# Clone 2
sim$add_mutant(name = "Clone 2", growth_rates = .2, death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size("Clone 2",1e2) #1e2
# plot_tissue(sim, num_of_bins=100)

# Clone 3
sim$update_rates("Clone 1", c(death = 1))
sim$add_mutant(name = "Clone 3", growth_rates = .3, death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$run_up_to_size("Clone 3",1e2) #1e2
# plot_tissue(sim, num_of_bins=100)

# Clone 4
sim$add_mutant(name = "Clone 4", growth_rates = .5, death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 4")
sim$run_up_to_size("Clone 4",1e4)
sim$get_clock()
# plot_tissue(sim, num_of_bins=100)

# ProCESS::plot_muller(sim)
# plot_state(sim)

# Sampling 1
sample_a <- sim$search_sample(min_num_of_cells=c("Clone 2" = 100, 'Clone 3' = 900), width = 50,height = 50)  # 10% Clone 2, 90% Clone 3 
sample_b <- sim$search_sample(min_num_of_cells=c("Clone 2" = 600, 'Clone 3' = 200, 'Clone 4' = 200), width = 50,height = 50)  # 60% Clone 2, 20% Clone 3, 20% Clone 4
sample_c <- sim$search_sample(min_num_of_cells=c("Clone 4" = 1000), width = 50,height = 50)  # 100% Clone 4

sample1_plot = plot_tissue(sim, num_of_bins=100) +
  ggplot2::geom_rect(ggplot2::aes(xmin=sample_a$lower_corner[1],
                xmax=sample_a$upper_corner[1],
                ymin=sample_a$lower_corner[2],
                ymax=sample_a$upper_corner[2]),
            color= 'white', fill='white')+
  ggplot2::geom_rect(ggplot2::aes(xmin=sample_b$lower_corner[1],
                xmax=sample_b$upper_corner[1],
                ymin=sample_b$lower_corner[2],
                ymax=sample_b$upper_corner[2]),
            color= 'white', fill='white')+
  ggplot2::geom_rect(ggplot2::aes(xmin=sample_c$lower_corner[1],
                xmax=sample_c$upper_corner[1],
                ymin=sample_c$lower_corner[2],
                ymax=sample_c$upper_corner[2]),
            color= 'white', fill='white')

ggsave(sample1_plot, filename = "/orfeo/cephfs/scratch/cdslab/antonelloa/ProCESS-examples/SCOUT/SPN07/sample_1.png")

sim$sample_cells("SPN07_1.1", sample_a$lower_corner, sample_a$upper_corner)
sim$sample_cells("SPN07_1.2", sample_b$lower_corner, sample_b$upper_corner)
sim$sample_cells("SPN07_1.3", sample_c$lower_corner, sample_c$upper_corner)

forest <- sim$get_sample_forest()
f1 = annotate_forest(forest = forest,
                     tree_plot = plot_forest(forest),
                     MRCAs = T,
                     samples = T)
ggsave(f1, filename="/orfeo/cephfs/scratch/cdslab/antonelloa/ProCESS-examples/SCOUT/SPN07/forest_1.png")

### Chemotherapy
sim$update_rates("Clone 2", c(death = 3))
sim$update_rates("Clone 3", c(death = 3))
sim$update_rates("Clone 4", c(death = 2))
sim$run_up_to_time(sim$get_clock() + 2)
# plot_tissue(sim)
sim$get_clock()

### Resistant clones birth
sim$add_mutant(name = "Clone 5",
               growth_rates = .27,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("Clone 4"), "Clone 5")
sim$run_up_to_size("Clone 5",1e2)
# plot_tissue(sim)

sim$add_mutant(name = "Clone 6",
               growth_rates = .5,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("Clone 5"), "Clone 6")
sim$run_up_to_size("Clone 6",1e4)
# plot_tissue(sim)

# plot_state(sim)
muller_plot = ProCESS::plot_muller(sim)
ggsave(muller_plot, filename="/orfeo/cephfs/scratch/cdslab/antonelloa/ProCESS-examples/SCOUT/SPN07/muller_plot.png")

### Sampling 2
sample_d <- sim$search_sample(c("Clone 5" = 1000), 50,50)
sample_e <- sim$search_sample(c("Clone 6" = 1000), 50,50)

plot_sampling2=plot_tissue(sim, num_of_bins=100) +
  ggplot2::geom_rect(ggplot2::aes(xmin=sample_d$lower_corner[1],
                xmax=sample_d$upper_corner[1],
                ymin=sample_d$lower_corner[2],
                ymax=sample_d$upper_corner[2]),
            color= 'white', fill='white')+
  ggplot2::geom_rect(ggplot2::aes(xmin=sample_e$lower_corner[1],
                xmax=sample_e$upper_corner[1],
                ymin=sample_e$lower_corner[2],
                ymax=sample_e$upper_corner[2]),
            color= 'white', fill='white')

ggsave(plot_sampling2, filename="/orfeo/cephfs/scratch/cdslab/antonelloa/ProCESS-examples/SCOUT/SPN07/sample_2.png")

sim$sample_cells("SPN07_2.1", sample_d$lower_corner, sample_d$upper_corner)
sim$sample_cells("SPN07_2.2", sample_e$lower_corner, sample_e$upper_corner)

forest <- sim$get_sample_forest()
f2 = annotate_forest(forest = forest,tree_plot = plot_forest(forest),MRCAs = T)
ggsave(f2, filename='/orfeo/cephfs/scratch/cdslab/antonelloa/ProCESS-examples/SCOUT/SPN07/forest_2.png')
forest$save(paste0(outdir, '/sample_forest.sff'))


