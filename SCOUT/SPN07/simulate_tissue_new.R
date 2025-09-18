rm(list = ls())
library(ProCESS)
library(dplyr)
library(ggplot2)

dir <- getwd()
outdir <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPN07_last/process" 
dir.create(outdir, recursive = T)
setwd(outdir)
unlink("SPN07", recursive = TRUE)
set.seed(1)

sim <- TissueSimulation(name = 'SPN07', seed = 1, save_snapshot=T, width = 1e3, height = 1e3) #save_snapshot=T
sim$history_delta <- .1
sim$death_activation_level <- 50

# Clone 1
sim$add_mutant(name = "Clone 1", growth_rates = .2, death_rates = .03)
sim$place_cell("Clone 1", 500, 500)
sim$run_up_to_size("Clone 1", 500) #1e2
plot_tissue(sim)
sim$get_clock()

# Clone 2
sim$add_mutant(name = "Clone 2", growth_rates = .3, death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("Clone 1"), "Clone 2")
sim$run_up_to_size("Clone 2", 500) #1e2
plot_tissue(sim, num_of_bins=100)
sim$get_clock()

# Clone 3
sim$update_rates("Clone 1", c(death = 1))
sim$add_mutant(name = "Clone 3", growth_rates = .5, death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 3")
sim$run_up_to_size("Clone 3", 500) #1e2
plot_tissue(sim, num_of_bins=100)
sim$get_clock()

# Clone 4
sim$add_mutant(name = "Clone 4", growth_rates = .9, death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("Clone 2"), "Clone 4")
sim$run_up_to_size("Clone 4",6e3)
sim$get_clock()
plot_tissue(sim, num_of_bins=100)

ProCESS::plot_muller(sim)

sample_a <- sim$search_sample(min_num_of_cells=c("Clone 2" = 100, 'Clone 3' = 900), width = 50,height = 50)  # 10% Clone 2, 90% Clone 3
sim$sample_cells("SPN07_1.1", sample_a$lower_corner, sample_a$upper_corner)

sample_b <- sim$search_sample(min_num_of_cells=c("Clone 2" = .6*200, 'Clone 3' = .2*200, 'Clone 4' = .2*200), width = 50,height = 50)  # 60% Clone 2, 20% Clone 3, 20% Clone 4
sim$sample_cells("SPN07_1.2", sample_b$lower_corner, sample_b$upper_corner)

sample_c <- sim$search_sample(min_num_of_cells=c("Clone 4" = 1600), width = 45,height = 45)  # 100% Clone 4
sim$sample_cells("SPN07_1.3", sample_c$lower_corner, sample_c$upper_corner)

sample1_plot = plot_tissue(sim, num_of_bins=100) +
  ggplot2::geom_rect(ggplot2::aes(xmin=sample_a$lower_corner[1],
                                  xmax=sample_a$upper_corner[1],
                                  ymin=sample_a$lower_corner[2],
                                  ymax=sample_a$upper_corner[2]),
                     color= 'black', fill=NA)+
  ggplot2::geom_rect(ggplot2::aes(xmin=sample_b$lower_corner[1],
                                  xmax=sample_b$upper_corner[1],
                                  ymin=sample_b$lower_corner[2],
                                  ymax=sample_b$upper_corner[2]),
                     color= 'black', fill=NA)+
  ggplot2::geom_rect(ggplot2::aes(xmin=sample_c$lower_corner[1],
                                  xmax=sample_c$upper_corner[1],
                                  ymin=sample_c$lower_corner[2],
                                  ymax=sample_c$upper_corner[2]),
                     color= 'black', fill=NA)


forest <- sim$get_sample_forest()
f1 = annotate_forest(forest = forest,
                     tree_plot = plot_forest(forest),
                     MRCAs = T,
                     samples = T)

print(paste0('Before treatment ', sim$get_clock()))

### Chemotherapy
sim$update_rates("Clone 2", c(growth = 0, death = 2))
sim$update_rates("Clone 3", c(growth = 0, death = 2))
sim$update_rates("Clone 4", c(growth = 0, death = 2))
sim$run_up_to_time(sim$get_clock() + 2)
plot_tissue(sim)


### Resistant clones birth
sim$add_mutant(name = "Clone 5",
               growth_rates = .4,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("Clone 4"), "Clone 5")
sim$run_up_to_size("Clone 5",500)
plot_tissue(sim)
print(paste0('After treatment ', sim$get_clock()))

sim$add_mutant(name = "Clone 6",
               growth_rates = .6,
               death_rates = .03)
sim$mutate_progeny(sim$choose_cell_in("Clone 5"), "Clone 6")
sim$run_up_to_size("Clone 6",8000)
plot_tissue(sim)

muller_plot = ProCESS::plot_muller(sim)
muller_plot

### Sampling 2
#sample_d <- sim$search_sample(c("Clone 5" = 1000, 'Clone 6' = 0), 50,50)
sample_d_lower_corner = c(390, 400)
sample_d_upper_corner = c(440, 450)

sample_e_lower_corner = c(350, 280)
sample_e_upper_corner = c(400, 330)

sim$sample_cells("SPN07_2.1", sample_d_lower_corner, sample_d_upper_corner)
sim$sample_cells("SPN07_2.2", sample_e_lower_corner, sample_e_upper_corner)

plot_sampling2 = plot_tissue(sim, num_of_bins=100) +
  ggplot2::geom_rect(ggplot2::aes(xmin=sample_d_lower_corner[1],
                                  xmax=sample_d_upper_corner[1],
                                  ymin=sample_d_lower_corner[2],
                                  ymax=sample_d_upper_corner[2]),
                     color= 'black', fill=NA) +
  ggplot2::geom_rect(ggplot2::aes(xmin=sample_e_lower_corner[1],
                                  xmax=sample_e_upper_corner[1],
                                  ymin=sample_e_lower_corner[2],
                                  ymax=sample_e_upper_corner[2]),
                     color= 'black', fill=NA)


forest <- sim$get_sample_forest()
f2 = annotate_forest(forest = forest, tree_plot = plot_forest(forest), MRCAs = T)
forest$save(paste0(outdir, '/sample_forest.sff'))
# 
# samples_table <- function(sim, sample_forest) {
#   info = sim$get_samples_info() ## requested from either the simulation recovery or as saved table
# 
#   nodes = sample_forest$get_nodes()
#   clones = nodes %>%
#     dplyr::filter(!is.na(sample)) %>%
#     dplyr::group_by(sample, mutant) %>%
#     dplyr::pull(mutant) %>%
#     unique()
#   clones_of_origin = nodes %>%
#     dplyr::filter(!is.na(sample)) %>%
#     dplyr::group_by(sample, mutant) %>%
#     # dplyr::mutate(mutant = gsub(" ", "_", mutant)) %>%
#     dplyr::count(mutant) %>%
#     tidyr::pivot_wider(values_from = n, names_from = mutant, values_fill = 0) %>%
#     dplyr::rowwise() %>%
#     dplyr::mutate(Sample_Type = sum(c_across(clones) == 0)) %>%
#     dplyr::mutate(Sample_Type = ifelse(Sample_Type == (length(clones)-1), "Monoclonal", "Polyclonal")) %>%
#     dplyr::mutate(Total_Cells = sum(c_across(clones), na.rm = T)) %>%
#     dplyr::mutate(across(all_of(clones), ~ round(.x/Total_Cells,2), .names = "{.col} proportion"))
# 
#   samples_tb = dplyr::full_join(info, clones_of_origin, by = join_by("name" == "sample")) %>%
#     dplyr::select(!c("xmin","xmax","ymin","ymax","id","tumour_cells","tumour_cells_in_bbox")) %>%
#     dplyr::rename(Sample_ID=name) %>%
#     dplyr::rename(Samping_Time=time) %>%
#     dplyr::mutate(Samping_Time=round(Samping_Time,2)) %>%
#     dplyr::arrange(Sample_ID)
# 
# 
# 
#   return(samples_tb)
# }
# samples_table(sim, forest)
