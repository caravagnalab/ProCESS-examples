process BCFTOOLS_CONCAT {
    tag { "${meta.sample}" }
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bcftools:1.20--h8b25389_0':
        'biocontainers/bcftools:1.20--h8b25389_0' }"

    input:
    tuple val(meta), path(vcf_path1), path(vcf_tbi_path1),  path(vcf_path2), path(vcf_tbi_path2)

    output:
    tuple val(meta), path("*.vcf.gz"),	emit: vcf

    publishDir "${params.outdir}/${meta.spn}/somatic/${meta.coverage}x_${meta.purity}p/${meta.caller}/", mode: 'copy'

    script:
    """
    bcftools concat "${vcf_path1}" "${vcf_path2}" --allow-overlaps --threads ${task.cpus} -o ${meta.sample}.vcf.gz -Oz
    """
}