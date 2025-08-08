process BATTENBERG {
    tag "${meta.id}"
    label 'process_high_memory'
    label 'error_retry'
    
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://opengenomics/battenberg:2.2.9' :
        'docker.io/opengenomics/battenberg:2.2.9' }"
 
    input:
    tuple val(meta), path(normal_bam), path(tumour_bam), path(normal_bai), path(tumour_bai)

    output:
    tuple val(meta), path("plot/*"), emit: plot
    tuple val(meta), path("log/*"), emit: log
    tuple val(meta), path("impute/*"), emit: impute
    tuple val(meta), path("haplotype/*"), emit: haplotype
    tuple val(meta), path("alleleFrequencies/*"), emit: allelefreq
    tuple val(meta), path("beagle/*"), emit: beagle
    tuple val(meta), path("*.txt"), emit: files

    script:
    """
    #!/usr/bin/env Rscript
    library(Battenberg)

    analysis = "paired"
    SAMPLENAME = "${meta.id}"
    NORMALNAME = "${meta.patient}"
    NORMALBAM = "$normal_bam"
    SAMPLEBAM = "$tumour_bam"
    GENDER = "${meta.gender}"

    if (GENDER == 'XY'){
      IS.MALE = T
    } else {
      IS.MALE = F
    }

    RUN_DIR = "${meta.id}"

    SKIP_ALLELECOUNTING = FALSE
    SKIP_PREPROCESSING = FALSE
    SKIP_PHASING = FALSE
    NTHREADS = as.integer("${task.cpus}")
    PRIOR_BREAKPOINTS_FILE = NULL
    GENOMEBUILD = "hg38"

    supported_analysis = c("paired", "cell_line", "germline")
    if (!analysis %in% supported_analysis) {
            stop(paste0("Requested analysis type ", analysis, " is not available. Please provide either of ", paste(supported_analysis, collapse=" ")))
    }

    supported_genome_builds = c("hg19", "hg38")
    if (!GENOMEBUILD %in% supported_genome_builds) {
            stop(paste0("Provided genome build ", GENOMEBUILD, " is not supported. Please provide either of ", paste(supported_genome_builds, collapse=" ")))
    }


    BEAGLE_BASEDIR = "${params.battenber_ref}"

    ALLELECOUNTER = "alleleCounter"
    IMPUTE_EXE = "impute2"
    USEBEAGLE = T
    JAVAJRE = file.path(BEAGLE_BASEDIR, "beagle_chr/beagle.08Feb22.fa4.jar")

    IMPUTEINFOFILE = file.path(BEAGLE_BASEDIR, "impute_info.txt")
    G1000PREFIX_AC = file.path(BEAGLE_BASEDIR, "1000G_loci_hg38_chr/1kg.phase3.v5a_GRCh38nounref_allele_index_")
    G1000PREFIX = file.path(BEAGLE_BASEDIR, "1000G_loci_hg38_chr/1kg.phase3.v5a_GRCh38nounref_loci_")
    GCCORRECTPREFIX = file.path(BEAGLE_BASEDIR, "GC_correction_hg38_chr/1000G_GC_")
    REPLICCORRECTPREFIX = file.path(BEAGLE_BASEDIR, "RT_correction_hg38_chr/1000G_RT_")
    PROBLEMLOCI = file.path(BEAGLE_BASEDIR, "probloci_chr/probloci.txt.gz")

    beagleref.template = file.path(BEAGLE_BASEDIR, "beagle_chr/CHROMNAME.1kg.phase3.v5a_GRCh38nounref.vcf.gz")
    beagleplink.template = file.path(BEAGLE_BASEDIR, "beagle_chr/plink.CHROMNAME.GRCh38.map")
    beaglejar = file.path(BEAGLE_BASEDIR, "beagle_chr/beagle.08Feb22.fa4.jar")

    BEAGLEJAR = beaglejar
    BEAGLEREF.template = beagleref.template
    BEAGLEPLINK.template = beagleplink.template

    externalHaplotypeFile = NA
    externalhaplotypefile = NA
    snp6_reference_info_file = NA

    CHROM_COORD_FILE = file.path(BEAGLE_BASEDIR,"chromosome_coordinates_hg38_chr.txt")

    PLATFORM_GAMMA = 1
    PHASING_GAMMA = 1
    SEGMENTATION_GAMMA = 10
    SEGMENTATIIN_KMIN = 3
    PHASING_KMIN = 1
    CLONALITY_DIST_METRIC = 0
    ASCAT_DIST_METRIC = 1
    MIN_PLOIDY = 1.6
    MAX_PLOIDY = 4.8
    MIN_RHO = 0.1
    MIN_GOODNESS_OF_FIT = 0.63
    BALANCED_THRESHOLD = 0.51
    MIN_NORMAL_DEPTH = 10
    MIN_BASE_QUAL = 20
    MIN_MAP_QUAL = 35
    #CALC_SEG_BAF_OPTION = 1
    CALC_SEG_BAF_OPTION = 3
    USEBEAGLE=TRUE
    BEAGLE_MAX_MEM=15
    BEAGLENTHREADS=1
    BEAGLEWINDOW=40
    BEAGLEOVERLAP=4

    # Enable cairo device (needed to prevent 'X11 not available' errors)
    options(bitmapType='cairo')

    # Change to work directory and load the chromosome information
    battenberg (analysis=analysis,
              tumourname=SAMPLENAME,
              normalname=NORMALNAME,
              tumour_data_file=SAMPLEBAM,
              normal_data_file=NORMALBAM,
              ismale=IS.MALE,
              imputeinfofile=IMPUTEINFOFILE,
              g1000prefix=G1000PREFIX,
              g1000allelesprefix=G1000PREFIX_AC,
              gccorrectprefix=GCCORRECTPREFIX,
              repliccorrectprefix=REPLICCORRECTPREFIX,
              problemloci=PROBLEMLOCI,
              data_type="wgs",
              impute_exe=IMPUTE_EXE,
              allelecounter_exe=ALLELECOUNTER,
              usebeagle=USEBEAGLE, ##set to TRUE to use beagle
              beaglejar=BEAGLEJAR, ##path
              beagleref=BEAGLEREF.template, ##pathtemplate
              beagleplink=BEAGLEPLINK.template, ##pathtemplate
              beaglemaxmem=BEAGLE_MAX_MEM,
              beaglenthreads=BEAGLENTHREADS,
              beaglewindow=BEAGLEWINDOW,
              beagleoverlap=BEAGLEOVERLAP,
              javajre=JAVAJRE,
              nthreads=NTHREADS,
              platform_gamma=PLATFORM_GAMMA,
              phasing_gamma=PHASING_GAMMA,
              segmentation_gamma=SEGMENTATION_GAMMA,
              segmentation_kmin=SEGMENTATIIN_KMIN,
              phasing_kmin=PHASING_KMIN,
              clonality_dist_metric=CLONALITY_DIST_METRIC,
              ascat_dist_metric=ASCAT_DIST_METRIC,
              min_ploidy=MIN_PLOIDY,
              max_ploidy=MAX_PLOIDY,
              min_rho=MIN_RHO,
              min_goodness=MIN_GOODNESS_OF_FIT,
              uninformative_BAF_threshold=BALANCED_THRESHOLD,
              min_normal_depth=MIN_NORMAL_DEPTH,
              min_base_qual=MIN_BASE_QUAL,
              min_map_qual=MIN_MAP_QUAL,
              calc_seg_baf_option=CALC_SEG_BAF_OPTION,
              skip_allele_counting=SKIP_ALLELECOUNTING,
              skip_preprocessing=SKIP_PREPROCESSING,
              skip_phasing=SKIP_PHASING,
              prior_breakpoints_file=PRIOR_BREAKPOINTS_FILE,
              GENOMEBUILD=GENOMEBUILD,
              chrom_coord_file=CHROM_COORD_FILE,
              write_battenberg_phasing = F
    )

    plot_dir = "plot"
    plot_file = list.files(pattern = "*png\$", full.names = TRUE)
    if (!dir.exists(plot_dir)) dir.create(plot_dir, showWarnings=F, recursive=T)
    dest_files = file.path(plot_dir, basename(plot_file))
    file.rename(plot_file, dest_files)

    log_dir = "log"
    log_file = list.files(pattern = "*log\$", full.names = TRUE)
    if (!dir.exists(log_dir)) dir.create(log_dir, showWarnings=F, recursive=T)
    dest_files = file.path(log_dir, basename(log_file))
    file.rename(log_file, dest_files)

    beagle_dir = "beagle"
    beagle_file = list.files(pattern = "*beagle*", full.names = TRUE)
    if (!dir.exists(beagle_dir)) dir.create(beagle_dir, showWarnings=F, recursive=T)
    dest_files = file.path(beagle_dir, basename(beagle_file))
    file.rename(beagle_file, dest_files)

    allele_dir = "alleleFrequencies"
    allele_file = list.files(pattern = "*alleleFrequencies*", full.names = TRUE)
    if (!dir.exists(allele_dir)) dir.create(allele_dir, showWarnings=F, recursive=T)
    dest_files = file.path(allele_dir, basename(allele_file))
    file.rename(allele_file, dest_files)

    impute_dir = "impute"
    impute_file = list.files(pattern = "*impute*", full.names = TRUE)
    if (!dir.exists(impute_dir)) dir.create(impute_dir, showWarnings=F, recursive=T)
    dest_files = file.path(impute_dir, basename(impute_file))
    file.rename(impute_file, dest_files)

    haplotype_dir = "haplotype"
    haplotype_file = list.files(pattern = "*haplotyped*", full.names = TRUE)
    if (!dir.exists(haplotype_dir)) dir.create(haplotype_dir, showWarnings=F, recursive=T)
    dest_files = file.path(haplotype_dir, basename(haplotype_file))
    file.rename(haplotype_file, dest_files)

    """
}
