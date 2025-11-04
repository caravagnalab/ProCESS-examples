process PROCESS_PREPROCESS {
    tag { "${meta.sample}_${meta.caller}_chr${chromosome}" }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v2' :
        'docker.io/lvaleriani/process_validation:v2' }"
    input:
    tuple val(meta), val(chromosome), path(rds_file)

    output:
    tuple val(meta), val(chromosome), path("chr*.rds"),	emit: rds
    publishDir "${params.outdir}/${meta.spn}/somatic/${meta.coverage}x_${meta.purity}p/${meta.caller}/${meta.sample}/", mode: 'copy'

    script:
    """
    #!/usr/bin/env Rscript
    library(ProCESS)
    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/getters/sarek_getters.R")
    source("${projectDir}/bin/somatic/utils/process_utils.R")

    chromosome <- "$chromosome"
    process_seq_results("$meta.spn", "$meta.purity", "$meta.coverage", 
        chromosome, 
        outdir = "${params.outdir}",
        rds_path="$rds_file",
        sample_id="$meta.sample")
    """
}
