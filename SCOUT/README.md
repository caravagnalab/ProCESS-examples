# SPN simulations
To ensure uniformity in the creation of each SPN, certain requirements must be met.  

## Simulation parameters
### Number of cells to sample 
Remember to sample around **2k cells** per sample.

### Mutations rates
Set mutation rates in MutationEngine to:
- SNV: `1e-8` (except for the SPN02 that is hypermutant)
- ID: `1e-9`
- CNA: depends on tumour type (eg. SPN01 simulation have a lot of CN events and 1e-10 rate, low number of CN with 1e-12 rate)

Eg.
```
m_engine$add_mutant(mutant_name = "Clone 1",
                    passenger_rates = c(SNV = 1e-8, CNA = 1e-12, indel = 1e-9))
```


### Reference path
The references to generated the `MutationEngine` are stored in `/orfeo/cephfs/scratch/cdslab/shared/ProCESS/GRCh38`. Please when setting the `MutationEngine` do:

```{r}
setwd("/orfeo/cephfs/scratch/cdslab/shared/ProCESS/GRCh38")
```

### Set tumour type in MutationEngine
```
# check available tumour type
get_available_tumours_in(setup_code = "GRCh38") %>% head()


m_engine <- MutationEngine(setup_code = "GRCh38",
                           tumour_type = "COAD",
                           tumour_study = "US")
```

### Set gender of simulation
```
#return available germline subject with gender
m_engine$get_germline_subjects()

#select the one that you want based on gender that you need
m_engine$set_germline_subject("")
```

Male samples:
- `SPN01`
- `SPN03`
- `SPN04`
- `SPN06`

Female samples:
- `SPN02`
- `SPN05`
- `SPN07`

### Set pre-neoplastic mutations
```
phylo_forest <- m_engine$place_mutations(samples_forest, num_of_preneoplatic_SNVs = 800, num_of_preneoplatic_indels = 200)
```


## Saving snapshot and set seed
Remember to save the snapshot of the simulations and set the seed by running:  

```{r}
set.seed(12345)
sim <- SpatialSimulation(name = 'SPN01', seed = 12345, save_snapshot=TRUE)
```

## Naming requirments

### Clone nomenclature
The standardize clone names should be `Clone N` (e.g.  `Clone 1`) for all SPNs.

### Sample nomenclature
To standardize sample names, the following naming convention must be used:
```
SPN{X}_{T}.{N}
```
where:
- `X` is the SPN number (eg, 01, 02, ...)
- `T` is the sampling time (1 = first time point, 2 = second time point, etc.)
- `N` is the sample number for a given time point T (1 = first sample at T time point, 2 = second sample at T time point)

If you are working on SPN01 and are sampling your third sample at the second time point, you must name it: `SPN01_2.3`

## Required script and files
For each SPN, the following files must be present:
- `simulate_tissue.R`
- `simulate_mutation.R`

These files must be pushed to the ProCESS-examples GitHub repository in the corresponding folder SCOUT/SPNX.

## Required output files
The following output files must be saved in the directory:
`/orfeo/cephfs/scratch/cdslab/shared/SCOUT/SPNX/process` (replace X with the number of your SPN):

- `samples_forest.sff`
- `phylogenetic_forest.sff`
- `<sample>_cna.rds`
- snapshot folder

The `<sample>_cna.rds` has to be obtained by running the following command for each sample:
```{r}
phylo_forest$get_bulk_allelic_fragmentation()
```

Eg.
```
sample_names <- phylo_forest$get_samples_info()[["name"]]
lapply(sample_names,function(s){
    cna <- phylo_forest$get_bulk_allelic_fragmentation(s)
    saveRDS(file=paste0(outdir,"cna_data/",s,"_cna.rds"),object=cna)
})
```

## How to run Rscript with singularity image
To run the R scripts, a Singularity image with the latest version of ProCESS must be used. You can do this with the following script:

```{sh}
#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 50gb
#SBATCH --time=1:00:00
#SBATCH --output=ProCESS.out
#SBATCH --error=ProCESS.err

module load singularity
image="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/process_1.1.0.sif"
# change with your path to the simulate_tissue.R and simulate_mutation.R scripts
base="/orfeo/scratch/area/lvaleriani/races/ProCESS-examples/SCOUT/SPN01"

singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_tissue.R
singularity exec --bind /orfeo:/orfeo --no-home $image Rscript $base/simulate_mutation.R
```
