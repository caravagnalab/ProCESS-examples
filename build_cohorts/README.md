
# How to Build The Cohort from The Phylogenetic Forests
For every SPN, the cohort contains 13 sets of reads: 12 simulate tumor sequencing with coverages of 50X, 100X, 150X, and 200X, and purities of 0.3, 0.6, and 0.9. The last set corresponds to a 30X simulated sequencing of normal samples. The cells in the normal sample share the same mutations, i.e., the germline mutations.

For efficiency, not all 12 tumor read sets were directly produced by ProCESS. Instead, for each purity level, the generation of the 200X-coverage read set was split into 40 lots of 5X coverage each. The lots are pairwise different due to their sequencing simulation random seeds: the random seed for the $i$-th lot is $i$.

The 50X, 100X, 150X, and 200X coverage sets were assembled by grouping 10, 20, 30, and 40 lots, respectively.

The normal data were generated similarly, with each 30X-coverage read set constructed from 6 lots of 5X coverage.

Using the same scripts, the `.rds` files produced by `simulate_seq()` for each lot will be merged into a single file.

Finally, a sample sheet and a bash script for running Sarek will be created.

<!-- For every SPN, the cohort contains 13 sets of reads: 12 of them simulate 
the tumour sequencing with coverage in 50X, 100X, 150X, and 200X and purity
varying among 0.3, 0.6, and 0.9; the last one is the 30X simulated
sequencing of normal samples. The cells in the normal sample share the same
mutations, i.e., the germline mutations.

For efficiency purposes, not all the 12 tumour read sets were directly
produced by ProCESS. Instead, for each purity, we split the generation 
of the 200X-coverage read set into 40 5X-coverage lots.
The lots are pairwise different because of their sequencing simulation
random seeds: the $i$-th lot random seed is $i$.

The 50X, 100X, 150X, 200X coverage sets were assembled by grouping 10, 20, 
30, and 40 lots, respectively.

The normal data were produced in an analogous way by splitting the computation 
of each 30X-coverage read set into lots of 5X coverage. -->

<!-- The lot-based approach allowed us to distribute the read set generation
on a cluster: the more lots, the faster the generation. 
The lot generation were distributed on a computational cluster managed by SLURM.
The cluster consists of 8 nodes equipped with an AMD EPYC 7H12 64-core
processor with 512GB of RAM and 352GB of local scratch.
The number of lots, i.e., 40, is tailored on this computational facility. -->

## Before sumbitting the sequencing

<!-- Before submitting the main script, you **MUST** create a temporary folder in your `fast` directory that will be used to store the singularity temporary files. So please do:
```{sh}
mkdir /orfeo/cephfs/fast/cdslab/${USER}/tmp
```-->

You must configure your `.bashrc` to point to:

1. The path to the shared singularity cache directory, that will be used later to run nextflow pipelines (`/orfeo/cephfs/scratch/cdslab/shared/containers/singularity/sarek_tumourevo`)
2. The path to the `work/` folder that will contain all intermediate files of nextflow pipelines (`/orfeo/cephfs/fast/cdslab/${USER}/work`)

So add the following line to your `.bashrc`:

```{sh}
export NXF_SINGULARITY_CACHEDIR="/orfeo/cephfs/scratch/cdslab/shared/containers/singularity/sarek_tumourevo"
export NXF_WORK="/orfeo/cephfs/fast/cdslab/${USER}/work"
```
> [!IMPORTANT]  
> Remember to clean up your tmp folders once the main sequencing script as well as nextflow pipelines successfully finished.

## Generating the Tumour Cohort
This directory contains the Python script `benchmark_build_cohort.py`, which builds the cohort of an SPN using its phylogenetic forest file on Orfeo. The sequencing results will be stored in `scratch/` and once the sequencining is done, you should copy the results in Long Term Storage (`LTS`) in order to have them properly backuped. Please see [After sequencing](#after-sequencing) section

To run the script, create a `bash` executable file like the example below.  
`run_build_cohort.sbatch` with SPN01 parameters is provided as a reference:

```{sh}
#!/bin/bash
#SBATCH --partition=EPYC
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem 20gb
#SBATCH --time=96:00:00
#SBATCH --output=seq.out
#SBATCH --error=seq.err

module load singularity/3.10.4

# change them accordingly
user="cdslab"
spn="SPN01"

# change with your own absolute path
path="/orfeo/cephfs/scratch/cdslab/ggandolfi/Github/ProCESS-examples/build_cohorts"

# keep them as they are
partition=EPYC
phylo="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/process/phylo_forest.sff"
tmp="/orfeo/cephfs/fast/cdslab/${USER}/tmp_files"
image="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/process_1.1.0.sif"
out="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sequencing"
sarek_output_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/sarek"
tumourevo_output_dir="/orfeo/cephfs/scratch/cdslab/shared/SCOUT/${spn}/tumourevo"

$path/benchmark_build_cohort.py -P $partition -A $user -s $tmp -I $image $spn $phylo $out -C $path/orfeo.config -SD $sarek_output_dir -TD $tumourevo_output_dir
```


You need to modify the following variables in the script:  

- `user`: your Orfeo group.  
- `spn`: the name of the SPN.  
- `path`: absolute path to the `benchmark_build_cohort.py` script.  

Leave the following variables unchanged:  

- `tmp`: absolute path to a directory where temporary files will be written and later deleted.
- `partition`: the Orfeo HPC partition where the script will run.
- `phylo`: absolute path to the phylogenetic forest.  
- `image`: absolute path to the Singularity image.  
- `config`: absolute path to the Sarek config.  
- `out`: absolute path where sequencing results will be written.  
- `sarek_output_dir`: absolute path to the directory where Sarek output will be stored.  

After replacing the variables, ensure that `benchmark_build_cohort.py` is executable by running:  

```{sh}
chmod +x benchmark_build_cohort.py
```
Finally, submit the script using:
```{sh}
sbatch run_build_cohort.sbatch
```
## Output files

<!-- To generate the cohort, copy the scripts into a _writable_
directory in your home on Orfeo, log in on Orfeo, and execute the following commands.

```{sh}
# clean the node local scratch directory
./clean_local_scratches.sh EPYC <USER>

# build the cohort
nohup ./build_cohort.py -P EPYC -A <USER> <SPN> <PHYLO_FOREST> <OUTPUT_DIR> &
```

`<USER>` is your SLURM user, `<SPN>` is the name of the SPN (e.g., SPN01), `<PHYLO_FOREST>` is the absolute
path of the phylogenetic forest, and `<OUTPUT_DIR>` is the absolute path of the output
directory.

For faster executions, both `<PHYLO_FOREST>` and `<OUTPUT_DIR>` should be
located in one of the `/fast` subdirectories.

> [!CAUTION]
> The content of the `/fast` subdirectories is deleted every 30 days -->

<!-- The Python script creates three directory into the `<OUTPUT_DIR>`: 
`tumour`, `normal`, and `sarek`.

The directories `tumour` and `normal` organizes the produced data by purity.
The former contains three sub-directory: `purity_0.3`, `purity_0.6`, and
`purity_0.9`. The latter has only one sub-directory, i.e., `purity_1`, meaning 
that the normal data has no contaminants.
Each purity directory stores the sequencing dataframe output (sub-directory 
`data`), the BAM file (sub-directory `BAM`), the correspoding compressed
fastq file (sub-directory `FASTQ`), and the log files of each lot.

Finally, the directory `sarek` contains X files: one for each combination of coverage 
among 50X, 100X, 150X, and 200X and purity among 0.3, 0.6, and 0.9.
There are `.csv` file with sarek datasheet for an experiment whose tumour coverage and 
purity are those declared in the file name and whose normal coverage is 30X for both 
the alignement part and the variant calling one.
In addition there will be `.sh` files that are required for running `sarek`. -->

The Python script creates four directories inside `<OUTPUT_DIR>`:  
`tumour`, `normal`, `sarek` and `tumourevo`.

The `tumour` and `normal` directories organize the produced data by purity.  
- The `tumour` directory contains three subdirectories: `purity_0.3`, `purity_0.6`, and `purity_0.9`.  
- The `normal` directory has only one subdirectory, `purity_1`, meaning that the normal data has no contaminants.  

Each purity directory includes:  
- The sequencing dataframe output (`data` subdirectory).  
- The corresponding compressed FASTQ file (`FASTQ` subdirectory).  
- The log files for each lot (`log` subdirectory);
- The resource usage for command line processes (`TIME` subdirectory).

The `sarek` directory contains multiple files, one for each combination of:  
- Tumour coverage: 50X, 100X, 150X, and 200X.  
- Tumour purity: 0.3, 0.6, and 0.9.  

It includes:  
- `.csv` files containing the Sarek datasheet for experiments with the specified tumour coverage and purity, with a normal coverage fixed at 30X for both alignment and variant calling.  
- `.sh` files required for running `sarek`.  

The `tumourevo` directory contains multiple files, one for each combination of:  
- Tumour coverage: 50X, 100X, 150X, and 200X;  
- Tumour purity: 0.3, 0.6, and 0.9;
- Copy number caller and somatic mutation caller.

It includes:  
- `.csv` files containing the Sarek datasheet for experiments with the specified tumour coverage and purity, with a normal coverage fixed at 30X for both alignment and variant calling.  
- `.sh` files required for running `sarek`.  

```
SCOUT/SPN01/sequencing
├── normal
│   └── purity_1
│       ├── data                                            # normal data
│       │   ├── mutations                                   # the normal purity is 1
│       │   │   ├── seq_results_muts_merged_coverage_30x.rds
│       │   │   ├── seq_results_muts_SPN01_n0.rds
│       │   │   └── ...
│       │   ├── parameters
│       │   │   ├── seq_results_params_SPN01_n0.rds
│       │   │   └── ...
│       │   └── resources
│       │       ├── seq_results_resources_SPN01_n0.rds
│       │       └── ...
│       ├── FASTQ
│       │   ├── n0_normal_sample.R1.fastq.gz
│       │   ├── n0_normal_sample.R2.fastq.gz
│       │   ├── n0_normal_sample.singleton.fastq.gz
│       │   ├── n0_normal_sample.unpaired.fastq.gz
│       │   └── ...
│       ├── log
│       │   ├── lot_n0.log                                  # lot 00 log 
│       │   └── ...
│       ├── n0_BAM.done
│       ├── n0_final.done
│       ├── ...
│       └── TIME                                            # resource usage for bash commands
│           ├── out_fastq_1_n0_normal_sample
│           ├── ...
│           ├── out_samtools_merge_1_n0_chr_Y
│           ├── ...
│           ├── out_samtools_split_1_n0.bam
│           └── ...
├── sarek
│   ├── sarek_50x_0.3p.csv                                  # sarek datasheet for tumour coverage 50X and purity 0.3 
│   ├── ...
│   ├── sarek_mapping_50x_0.3p.sh                           # .sh file for running sarek mapping for tumour coverage 50X and purity 0.3 
│   ├── ...
│   ├── sarek_mapping_normal.sh                             # .sh file for running sarek mapping for normal sample
│   ├── sarek_normal.csv                                    # sarek datasheet for normal sample
│   ├── sarek_variant_calling_50x_0.3p.csv                  # sarek datasheet for variant calling for tumour coverage 50X and purity 0.3
│   ├── ...    
│   ├── sarek_variant_calling_50x_0.3p.sh                   # .sh file for running sarek variant calling for tumour coverage 50X and purity 0.3 
│   └── ...
├── tumour                                                  # tumour data
│   ├── purity_0.3                                          # tumour data for purity 0.3
│   │   ├── data
│   │   │   ├── mutations                                   # ProCESS sequencing result for lot 00
│   │   │   │   ├── seq_results_muts_SPN01_t0.rds
│   │   │   │   └── ...
│   │   │   │   ├── seq_results_muts_merged_50x_SPN01.rds   # ProCESS merged sequencing results for coverage 50X 
│   │   │   │   └── ...
│   │   │   ├── parameters
│   │   │   │   ├── seq_results_params_SPN01_t0.rds         # ProCESS sequencing parameters for lot 00
│   │   │   │   └── ...
│   │   │   └── resources               
│   │   │       ├── seq_results_resources_SPN01_t0.rds      # resources required for sequencing lot 00
│   │   │       └── ...
│   │   ├── FASTQ
│   │   │   ├── t0_SPN01_1.1.R1.fastq.gz                    # 5' reads of sample "SPN01_1.1" for lot 00
│   │   │   ├── t0_SPN01_1.1.R2.fastq.gz                    # 3' reads of sample "SPN01_1.1" for lot 00
│   │   │   ├── t0_SPN01_1.1.singleton.fastq.gz
│   │   │   ├── t0_SPN01_1.1.unpaired.fastq.gz
│   │   │   ├── t0_SPN01_1.2.R1.fastq.gz                    # 5' reads of sample "SPN01_1.2" for lot 00
│   │   │   ├── t0_SPN01_1.2.R2.fastq.gz                    # 3' reads of sample "SPN01_1.2" for lot 00
│   │   │   ├── t0_SPN01_1.2.singleton.fastq.gz
│   │   │   ├── t0_SPN01_1.2.unpaired.fastq.gz
│   │   │   ├── t0_SPN01_1.3.R1.fastq.gz                    # 5' reads of sample "SPN01_1.3" for lot 00
│   │   │   ├── t0_SPN01_1.3.R2.fastq.gz                    # 3' reads of sample "SPN01_1.3" for lot 00
│   │   │   ├── t0_SPN01_1.3.singleton.fastq.gz
│   │   │   ├── t0_SPN01_1.3.unpaired.fastq.gz
│   │   │   └── ...
│   │   ├── log
│   │   │   ├── lot_t0.log                                  # lot 00 log 
│   │   │   └── ...
│   │   ├── t0_BAM.done
│   │   ├── t0_final.done
│   │   ├── ...
│   │   └── TIME                                            # resource usage for bash commands
│   │       ├── out_fastq_0.3_t0_SPN01_1.1                  
│   │       ├── out_fastq_0.3_t0_SPN01_1.2
│   │       ├── out_fastq_0.3_t0_SPN01_1.3
│   │       ├── ...
│   │       ├── out_samtools_merge_0.3_t0_chr_Y
│   │       ├── ...
│   │       ├── out_samtools_split_0.3_t0.bam
│   │       └── ...
│   ├── purity_0.6                                           # tumour data for purity 0.6
│       └── ...
│   ├── purity_0.9                                           # tumour data for purity 0.9 
│       └── ...
└── tumourevo
    ├── tumourevo_50x_0.3p_freebayes_ascat.csv              # tumourevo datasheet for freebayes and ascat for tumour coverage 50X and purity 0.3      
    ├── tumourevo_50x_0.3p_freebayes_ascat.sh               # .sh file for running tumourevo for freebayes and ascat for tumour coverage 50X and purity 0.3      
    ├── tumourevo_50x_0.3p_mutect2_ascat.csv               # tumourevo datasheet for mutect2 and ascat for tumour coverage 50X and purity 0.3 
    ├── tumourevo_50x_0.3p_mutect2_ascat.sh                # .sh file for running tumourevo for mutect2 and ascat for tumour coverage 50X and purity 0.3   
    ├── tumourevo_50x_0.3p_strelka_ascat.csv               # tumourevo datasheet for strelka and ascat for tumour coverage 50X and purity 0.3 
    ├── tumourevo_50x_0.3p_strelka_ascat.sh                # .sh file for running tumourevo for strelka and ascat for tumour coverage 50X and purity 0.3
    └── ...           
```
<!--
## After sequencing
Once sequencing is terminated, please copy the entire folder in LTS folder using rsync. (ADD MORE DETAILS, LIKE THE CMD, THE OUTPUT FOLDER,
THE NEED FOR RESOURCE).

> [!IMPORTANT]
> If the sequencing went straight without any issue, please remove the `BAM`
> folder in the `SCOUT/{SPNID}/sequencing` folder.
-->

## Run nextflow pipelines

Once the raw sequencing data have been generated and saved in the correct place,
it is time to run two nextflow pipelines: nf-core/sarek and nf-core/tumourevo.
Everyone involved in the SCOUT cohort generation will use the same `nextflow` executable available
in a shared folder at `/orfeo/cephfs/scratch/cdslab/shared/SCOUT/nextflow`,
ensuring that the same version is used (nextflow version 24.10.3.5933).
Both nextflow pipelines will be run with the same configuration file available at [orfeo.config](orfeo.config).

In summary you MUST run:

1. Mapping and preprocessing of raw normal sequencing reads (see [Mapping section](#mapping));
2. Mapping and preprocessing of raw tumour sequencing reads (see [Mapping section](#mapping)) for each coverage-purity combination;
3. Variant calling for selected tools starting from mapped reads (see [Variant Calling section](#variant-calling)) for each coverage-purity combination;
4. nf-core/tumourevo for three different combinations of SNVs caller (strelka, mutect2, freebayes) and CNA caller (ASCAT) and for each coverage-purity combination (see [Tumourevo section](#tumourevo)).

### Before running nextflow pipelines

Before running any nextflow pipeline, you need to configure the `SINGULARITY_TMPDIR` in a proper way, which in this case is specific for Orfeo cluster only. You will need to generate two files, `prolog.sh` and `epilog.sh` in your `/orfeo/cephfs/fast/cdslab/${USER}/` folder. Create and add to the `prolog.sh` the following lines:

```{bash}
#!/bin/bash
echo "questo è il prolog $(hostname), jobid ${SLURM_JOB_ID}"
mkdir -p /tmp/${SLURM_JOB_ID}
export TMPDIR="/tmp/${SLURM_JOB_ID}"
export NXF_TMP=$TMPDIR
export SINGULARITY_TMPDIR=$TMPDIR
export TMP=$TMPDIR
```

Create and add to the `epilog.sh` the following lines:

```{bash}
echo "questo è il epilog $(hostname), jobid ${SLURM_JOB_ID}"
echo $SINGULARITY_TMPDIR
echo $TMPDIR
rm -rf /tmp/${SLURM_JOB_ID}
```

In this way, all temporary files generated by `singularity` will be saved
into a job specific folder into the `/tmp` of each node. As soon as each process finishes
the `/tmp/${SLURM_JOB_ID}` will be removed, in order not to fill the inodes of each node.

### Mapping  
As described previously, the `sarek` directory will store both `.csv` and `.sh` files for running Sarek steps.  

> [!IMPORTANT]
> Before running `nf-core/sarek` install version 3.5.1 by running:
> `nextflow pull nf-core/sarek -r 3.5.1`

The first files generated will be `sarek_normal.csv` and `sarek_mapping_normal.sh`. Once these files are created, you can start the mapping of the normal sample by running:  

```{sh}
sbatch sarek_mapping_vc_normal.sh
```

Once the first purity combination (0.3) is completed, the following files will be generated:
- `sarek_50x_0.3p.csv`, `sarek_100x_0.3p.csv`, `sarek_150x_0.3p.csv`, `sarek_200x_0.3p.csv`
- `sarek_mapping_50x_0.3p.sh`, `sarek_mapping_100x_0.3p.sh`, `sarek_mapping_150x_0.3p.sh`, `sarek_mapping_200x_0.3p.sh`

You can then submit the following jobs:
```
sbatch sarek_mapping_50x_0.3p.sh
sbatch sarek_mapping_100x_0.3p.sh
sbatch sarek_mapping_150x_0.3p.sh
sbatch sarek_mapping_200x_0.3p.sh
```

At the same time, `sarek_variant_calling_{C}x_0.3p.csv` and `sarek_variant_calling_{C}x_0.3p.sh` will also be generated, but they should **NOT** be run yet.

After the 0.3 purity mapping is completed, the 0.6 and 0.9 purity mappings will be processed in sequence. After each purity level is completed, you can submit the corresponding mapping scripts in the same way.

### Variant Calling
Once the mapping for a given purity level is completed, you can start the variant calling. For example, after the mapping of 50X_0.3p is finished, run:

```{sh}
sbatch sarek_variant_calling_50x_0.3p.sh 
```

Repeat this step for each combination of coverage and purity.

### Tumourevo
Once the variant calling for a given purity level is completed, you can start running tumourevo with different tools. For example, after the variant calling of 50X_0.3p with ASCAT and mutect2 is finished, run:

```{sh}
sbatch tumourevo_50x_0.3p_mutect2_ascat.sh
```

Repeat this step for each combination of coverage and purity and different somatic mutation callers (mutect2, strelka and freebayes).
