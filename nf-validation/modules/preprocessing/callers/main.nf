process SOMATIC_PREPROCESS {

    tag { "${row.sample}_${caller}_chr${chromosome}" }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v1' :
        'docker.io/lvaleriani/process_validation:v1' }"

    input:
    tuple val(row), val(chromosome), val(caller), val(mut_type),path(vcf_file)

    output:
    // tuple val(meta), path("*.rds"),                            emit: rds
    // path "${row.spn}_${row.coverage}_${row.purity}_chr${chromosome}_${caller}_output.txt"
    tuple val(row), val(chromosome), val(caller),path("**/chr*.rds"),	emit: rds
    publishDir "${params.outdir}/somatic/${row.coverage}x_${row.purity}p/${caller}/${row.sample}/", mode: 'copy'

    script:
    """
    #!/usr/bin/env Rscript
    library(ProCESS)
    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/getters/sarek_getters.R")
    source("${projectDir}/bin/somatic/utils/mutect_utils.R")
    source("${projectDir}/bin/somatic/utils/strelka_utils.R")
    source("${projectDir}/bin/somatic/utils/freeBayes_utils.R")

    #write("$row.spn", file = "${row.spn}_${row.coverage}_${row.purity}_chr${chromosome}_${caller}_output.txt")
    chromosome <- "$chromosome"
    caller="${caller}"
    if (caller=="mutect2"){

      process_mutect2_results("$row.spn", "$row.purity", "$row.coverage",
         chromosome, base_path = "$row.directory", outdir = "${params.outdir}",
         vcf_path="$vcf_file",sample_id="$row.sample")
    } else if (caller=="strelka"){
          process_strelka_results("$row.spn", "$row.purity", "$row.coverage",
             chromosome, base_path = "$row.directory", outdir = "${params.outdir}",
             vcf_path="$vcf_file",sample_id="$row.sample",mut_type="${mut_type}")
    } else if (caller=="freebayes"){
          process_freebayes_results("$row.spn", "$row.purity", "$row.coverage",
             chromosome,base_path = "$row.directory", outdir = "${params.outdir}",
             vcf_path="$vcf_file",sample_id="$row.sample",
             pass_quality = 20, min_vaf = 0.01, max_normal_vaf = 0.02)
    }
    """
}
