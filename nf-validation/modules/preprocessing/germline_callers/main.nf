process GERMLINE_PREPROCESS {

    tag { "${meta.spn}-${meta.sample}-${meta.caller}-chr${chromosome}" }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v1' :
        'docker.io/lvaleriani/process_validation:v1' }"

    input:
    tuple val(meta), val(chromosome), path(vcf_file)

    output:
    tuple val(meta), val(chromosome), path("chr*.rds"),	emit: rds
    publishDir "${params.outdir}/${meta.spn}/germline/${meta.sample}/${meta.caller}/rds", mode: 'copy'

    script:
    """
    #!/usr/bin/env Rscript
    library(ProCESS)
    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/getters/sarek_getters.R")
    source("${projectDir}/bin/Germline/vcf_parser.R")

    chromosome <- "$chromosome"
    caller="${meta.caller}"
    outfile = paste0("chr",chromosome,"_${meta.sample}.rds")

    if (caller=="haplotypecaller"){
      rds_vcf = parse_HaplotypeCaller(file="$vcf_file",out_file=outfile, save = T)
    } else if (caller=="freebayes"){
      rds_vcf = parse_freebayes(file="$vcf_file",out_file=outfile, save = T,cutoff=0.3)
    } else if (caller=="strelka"){
      rds_vcf = parse_strelka(file="$vcf_file",out_file=outfile, save = T)
    } 
    """
}
