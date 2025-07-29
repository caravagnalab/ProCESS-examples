process BCFTOOLS_VIEW {
    tag { "${meta.spn}-${meta.caller}-chr${chromosome}" }

    input:
    tuple val(meta), val(chromosome), file(vcf_path)

    output:
    tuple val(meta), val(chromosome), path("*.vcf.gz"),	emit: chr_vcf

    publishDir { 
        if ("${meta.type}" == "1") {
            return "${params.outdir}/${meta.spn}/somatic/${meta.coverage}x_${meta.purity}p/${meta.caller}/vcf"
        } else if ("${meta.type}" == "0") {
            return "${params.outdir}/${meta.spn}/germline/${meta.sample}/${meta.caller}/vcf"
        }
    },     mode: 'copy'

    script:

    if ("${meta.caller}"=="mutect2"){
      out_vcf = "chr${chromosome}_${meta.spn}.vcf"
    } else {
      out_vcf = "chr${chromosome}_${meta.sample}.vcf"
    }

    """
    zgrep ^# ${vcf_path} > header_vcf.tmp
    zgrep -v ^# ${vcf_path} | grep ^chr${chromosome} > filtered_vcf.tmp
    cat header_vcf.tmp filtered_vcf.tmp > ${out_vcf}
    gzip ${out_vcf}
    """
}
