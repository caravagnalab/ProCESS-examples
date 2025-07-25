process BCFTOOLS_VIEW {
    tag { "${status} ${caller} chr${chromosome}" }

    input:
    tuple val(row), val(chromosome), val(caller), val(status), val(mut_type), file(vcf_path), file(vcf_tbi_path)

    output:
    tuple val(row), val(chromosome), val(caller), val(mut_type), path("*.vcf.gz"),	emit:chr_vcf
    // path ""
    // path "results/${row.spn}/validation/germline/vcf/chr${chromosome}_normal_sample.${caller}.vcf.gz"

    // publishDir "${params.outdir}/${status}/vcf", mode: 'copy'
    publishDir { 
        if (status == 'somatic') {
            return "${params.outdir}/somatic/${row.coverage}x_${row.purity}p/${caller}/vcf"
        } else if (status == 'germline') {
            return "${params.outdir}/germline/${row.sample}/${caller}/vcf"
        }
    },     mode: 'copy'

    script:

    if (caller=="mutect2"){
      //out_vcf = "chr${chromosome}_${row.spn}.vcf.gz"
      out_vcf = "chr${chromosome}_${row.spn}.vcf"
    } else if (caller=="strelka" && status == "somatic"){
      if (mut_type=="SNV"){
        //out_vcf = "chr${chromosome}_${row.sample}.snv.vcf.gz"
        out_vcf = "chr${chromosome}_${row.sample}.snv.vcf"
      } else{
        //out_vcf = "chr${chromosome}_${row.sample}.indel.vcf.gz"
        out_vcf = "chr${chromosome}_${row.sample}.indel.vcf"
      }
    } else {
      //out_vcf = "chr${chromosome}_${row.sample}.vcf.gz"
      out_vcf = "chr${chromosome}_${row.sample}.vcf"
    }

    """
    zgrep ^# ${vcf_path} > header_vcf.tmp
    zgrep -v ^# ${vcf_path} | grep ^chr${chromosome} > filtered_vcf.tmp
    cat header_vcf.tmp filtered_vcf.tmp > ${out_vcf}
    gzip ${out_vcf}
    #bcftools view ${vcf_path} --regions chr${chromosome} --threads ${task.cpus} -o ${out_vcf} -Oz
    """
}
