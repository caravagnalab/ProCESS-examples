library(dplyr)

COV = c(50, 100, 150, 200)
PUR = c(0.3, 0.6, 0.9)
NORMAL_COV = 30
########################
# change with SPN name
spn <- "SPN07"
########################
workdir  <- "/orfeo/scratch/cdslab/antonelloa/ProCESS-examples/report"
output <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"
input <- "/orfeo/cephfs/scratch/cdslab/shared/SCOUT"

dir.create(paste0(output,'/', spn, '/report/'), recursive = T, showWarnings = F)
setwd(workdir)

print('Rendering SPN report:')
params <- list(
  spn = spn,
  files = list(
    base_dir = paste0(input, '/', spn, '/process'),
    sample_forest = paste0(input, '/', spn, '/process/sample_forest.sff'),
    phylo_forest = paste0(input, '/', spn, '/process/phylo_forest.sff'),
    seq_res = paste0(input, '/', spn, '/sequencing')
  ),
  cna_dir = paste0(input, '/', spn, '/process/cna_data')
)
rmarkdown::render("Report.Rmd", 
                  params = params, 
                  output_file = paste0(output, '/', spn, '/report/', spn, '.html'))

print('Done')