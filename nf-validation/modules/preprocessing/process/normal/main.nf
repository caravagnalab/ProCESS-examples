process PROCESS_GERMLINE_PREPROCESS {

    tag { "${meta.spn}_${caller}" }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v1' :
        'docker.io/lvaleriani/process_validation:v1' }"

    input:
    tuple val(meta), val(caller), path(rds_file)

    output:
    tuple val(meta), val(caller), path("*.rds"), emit: rds
    publishDir "${params.outdir}/${meta.spn}/germline/${meta.sample}/${caller}/", mode: 'copy'

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
    saveRDS(process_normal,paste0("${meta.spn}",".rds"))
    """
}
