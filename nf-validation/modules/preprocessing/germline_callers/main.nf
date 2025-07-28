process GERMLINE_PREPROCESS {

    tag { "${row.sample}_${caller}_chr${chromosome}" }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v1' :
        'docker.io/lvaleriani/process_validation:v1' }"

    input:
    tuple val(row), val(chromosome), val(caller), val(mut_type), path(vcf_file)

    output:
    // tuple val(meta), path("*.rds"),                            emit: rds
    // path "${row.spn}_${row.coverage}_${row.purity}_chr${chromosome}_${caller}_output.txt"
    tuple val(row), val(chromosome), val(caller),path("chr*.rds"),	emit: rds
    publishDir "${params.outdir}/germline/${row.sample}/${caller}/rds", mode: 'copy'

    script:
    """
    #!/usr/bin/env Rscript
    library(ProCESS)
    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/getters/sarek_getters.R")
    source("${projectDir}/bin/Germline/vcf_parser.R")

    chromosome <- "$chromosome"
    caller="${caller}"
    outfile = paste0("chr",chromosome,"_${row.sample}.rds")

    if (caller=="haplotypecaller"){
      rds_vcf = parse_HaplotypeCaller(file="$vcf_file",out_file=outfile, save = T)
    } else if (caller=="freebayes"){
      rds_vcf = parse_freebayes(file="$vcf_file",out_file=outfile, save = T,cutoff=0.3)
    } else if (caller=="strelka"){
      rds_vcf = parse_strelka(file="$vcf_file",out_file=outfile, save = T)
    } 
    """
}
