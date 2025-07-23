process GERMLINE_REPORT {
    tag { "${caller} chr${chromosome}" }

    input:
    tuple val(row), val(chromosome), val(caller), path(vcf_file)

    output:
    path "chr${chromosome}_normal_sample.${caller}.vcf.gz"
    // path "results/${row.spn}/validation/germline/vcf/chr${chromosome}_normal_sample.${caller}.vcf.gz"

    publishDir "results/${row.spn}/validation/germline/vcf", mode: 'copy'
    script:
    """

    cd ${projectDir}/bin/Germline
    Rscript ${DIRECTORY}/Germline/compare_caller.R -s ${SPN}

    bcftools view \$vcf_path --regions chr${chromosome} --threads ${task.cpus} -o chr${chromosome}_normal_sample.${caller}.vcf.gz -Oz
    """
}
