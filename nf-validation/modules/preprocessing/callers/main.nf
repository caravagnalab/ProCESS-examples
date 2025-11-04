process SOMATIC_PREPROCESS {

    tag { "${meta.sample}_${meta.caller}_chr${chromosome}" }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v2' :
        'docker.io/lvaleriani/process_validation:v2' }"

    input:
    tuple val(meta), val(chromosome), path(vcf_file)

    output:
    tuple val(meta), val(chromosome), path("chr*.rds"),	emit: rds
    publishDir "${params.outdir}/${meta.spn}/somatic/${meta.coverage}x_${meta.purity}p/${meta.caller}/${meta.sample}/", mode: 'copy'

    script:
    """
    #!/usr/bin/env Rscript
    library(ProCESS)
    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/getters/sarek_getters.R")
    source("${projectDir}/bin/somatic/utils/mutect_utils.R")
    source("${projectDir}/bin/somatic/utils/strelka_utils.R")
    source("${projectDir}/bin/somatic/utils/freeBayes_utils.R")

    chromosome <- "$chromosome"
    caller="${meta.caller}"
    if (caller=="mutect2"){
      process_mutect2_results("$meta.spn", "$meta.purity", "$meta.coverage",
         chromosome, 
         outdir = "${params.outdir}",
         vcf_path="$vcf_file",
         sample_id="$meta.sample")

    } else if (caller=="strelka"){
          process_strelka_results("$meta.spn", "$meta.purity", "$meta.coverage",
             chromosome, 
             outdir = "${params.outdir}",
             vcf_path="$vcf_file",
             sample_id="$meta.sample")

    } else if (caller=="freebayes"){
          process_freebayes_results("$meta.spn", "$meta.purity", "$meta.coverage",
             chromosome,
             outdir = "${params.outdir}",
             vcf_path="$vcf_file",
             sample_id="$meta.sample",
             pass_quality = 20, 
             min_vaf = 0.01, 
             max_normal_vaf = 0.02)
    }
    """
}
