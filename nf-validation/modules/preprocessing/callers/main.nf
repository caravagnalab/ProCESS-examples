process SOMATIC_PREPROCESS {

    tag { "${row.sample}_${caller}_chr${chromosome}" }

    input:
    tuple val(row), val(chromosome), val(caller), path(vcf_file)

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
    source("${projectDir}/bin/somatic/utils/process_utils.R")

    #write("$row.spn", file = "${row.spn}_${row.coverage}_${row.purity}_chr${chromosome}_${caller}_output.txt")
    chromosome <- "$chromosome"
    caller="${caller}"
    if (caller=="mutect2"){

      process_mutect2_results("$row.spn", "$row.purity", "$row.coverage",
         chromosome, base_path = "$row.directory", outdir = "${params.outdir}",
         vcf_path="$vcf_file",sample_id="$row.sample")
    } else if (caller=="strelka"){
         process_strelka_results(()
    }
    """
}
