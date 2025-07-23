nextflow.enable.dsl=2

include { BCFTOOLS_VIEW } from '../../modules/bcftools/view/main.nf'
include { SOMATIC_PREPROCESS } from '../../modules/somatic/preprocessing/main.nf'
workflow SOMATIC_VALIDATION {

    take:
    t_sample_ch

    main:

    t_sample_ch
        .flatMap { row ->
            def chromosomes = (1..22).collect { it.toString() } + ['X']
            if (row.sex == 'XY') {
                chromosomes += 'Y'
            }
            def callers = ['mutect2'] // , 'strelka']
            def status = "somatic"
            chromosomes.collectMany { chr ->
                callers.collect { caller ->
                    [row, chr, caller,status]
                }
            }

        }
        .map { row, chr, caller, status ->
            def vcf_path
            if (caller == 'mutect2') {
                vcf_path = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.spn}/${row.spn}.mutect2.filtered.vcf.gz"
                vcf_tbi_path = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.spn}/${row.spn}.mutect2.filtered.vcf.gz.tbi"
            } else if (caller == 'freebayes') {
                vcf_path = "${row.directory}/${row.spn}/sarek/normal_sample/variant_calling/${caller}/normal_sample/normal_sample.freebayes.vcf.gz"
            }
            return [row, chr, caller, status,file(vcf_path), file(vcf_tbi_path)]
        }
        .set { somatic_jobs_ch }
    
    //viewed_vcfs=BCFTOOLS_VIEW(somatic_jobs_ch)

    BCFTOOLS_VIEW(somatic_jobs_ch)
    SOMATIC_PREPROCESS(BCFTOOLS_VIEW.out.chr_vcf)
    
    r_output = SOMATIC_PREPROCESS.out.rds

    emit:
    // viewed_vcfs
    r_output
}
