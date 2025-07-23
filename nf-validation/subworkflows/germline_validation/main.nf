nextflow.enable.dsl=2

include { BCFTOOLS_VIEW } from '../../modules/bcftools/view/main.nf'
// include { run_rscript }    from '../modules/run_rscript.nf'

workflow GERMLINE_VALIDATION {

    take:
    normal_sample_ch

    main:

    normal_sample_ch
        .flatMap { row ->
            def chromosomes = (1..22).collect { it.toString() } + ['X']
            if (row.sex == 'XY') {
              chromosomes += 'Y' 
            }
            def callers = ['haplotypecaller', 'freebayes', 'strelka']
            def status = "germline"
            chromosomes.collectMany { chr ->
                callers.collect { caller ->
                    [row, chr, caller,status]
                }
            }
        }
        .map { row, chr, caller,status ->
            def vcf_path
            if (caller == 'strelka') {
                vcf_path = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.sample}/${row.sample}.${caller}.variants.vcf.gz"
                vcf_tbi_path = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.sample}/${row.sample}.${caller}.variants.vcf.gz.tbi"
            } else if (caller == 'haplotypecaller') {
                vcf_path = "${row.directory}/${row.spn}/sarek/normal/variant_calling/${caller}/${row.sample}/${row.sample}.${caller}.filtered.vcf.gz"
                vcf_tbi_path = "${row.directory}/${row.spn}/sarek/normal/variant_calling/${caller}/${row.sample}/${row.sample}.${caller}.filtered.vcf.gz.tbi"
            } else if (caller == 'freebayes') {
                vcf_path = "${row.directory}/${row.spn}/sarek/normal/variant_calling/${caller}/${row.sample}/${row.sample}.${caller}.vcf.gz"
                vcf_tbi_path = "${row.directory}/${row.spn}/sarek/normal/variant_calling/${caller}/${row.sample}/${row.sample}.${caller}.vcf.gz.tbi"
            }
            return [row, chr, caller, status,file(vcf_path), file(vcf_tbi_path)]
        }
        .set { normal_jobs_ch }

    BCFTOOLS_VIEW(normal_jobs_ch)

    viewed_vcfs = BCFTOOLS_VIEW.out.chr_vcf

    // r_output = run_rscript(viewed_vcfs)

    emit:
    viewed_vcfs
    // r_output
}
