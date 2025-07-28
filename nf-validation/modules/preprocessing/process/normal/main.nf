process PROCESS_GERMLINE_PREPROCESS {

    tag { "${row.sample}_${caller}_chr${chromosome}" }

    input:
    tuple val(row), val(caller), path(rds_file)

    output:
    tuple val(row), val(caller),path("*.rds"),	emit: rds
    publishDir "${params.outdir}/germline/${row.sample}/${caller}/", mode: 'copy'

    script:
    """
    #!/usr/bin/env Rscript
    library(ProCESS)
    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/getters/sarek_getters.R")
    source("${projectDir}/bin/somatic/utils/process_utils.R")

    process_normal <- readRDS("$rds_file") %>% 
        filter(classes =='germinal') %>% 
        dplyr::mutate(chr = paste0('chr', chr),
                      mutationID = paste(chr,chr_pos, sep = ':')) %>% 
        dplyr::rename(BAF = normal_sample.VAF,
               DP = normal_sample.coverage,
               NV = normal_sample.occurrences)
    saveRDS(process_normal,paste0("${row.sample}",".rds"))
    """
}
