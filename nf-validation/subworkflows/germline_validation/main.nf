nextflow.enable.dsl=2

include { BCFTOOLS_VIEW } from '../../modules/bcftools/view/main.nf'
include { GERMLINE_PREPROCESS } from '../../modules/preprocessing/germline_callers/main.nf'
include { PROCESS_GERMLINE_PREPROCESS } from '../../modules/preprocessing/process/normal/main.nf'
//include { GENERATE_GERMLINE_REPORT }    from '../../modules/germline_report/main.nf'

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
            def mut_type = "all"
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
            return [row, chr, caller, status,mut_type,file(vcf_path), file(vcf_tbi_path)]
        }
        .set { normal_jobs_ch }

    normal_sample_ch
        .flatMap { row ->
            def callers = ['process']
            def status = "germline"
            callers.collect { caller ->
                [row, caller, status]
            }
        }
        .map { row, caller, status ->
            def rds_path = "${row.directory}/${row.spn}/sequencing/normal/purity_1/data/mutations/seq_results_muts_merged_coverage_30x.rds"
            return [row, caller, file(rds_path)]
        }
        .set { process_germline_jobs_ch }
    BCFTOOLS_VIEW(normal_jobs_ch)
    GERMLINE_PREPROCESS(BCFTOOLS_VIEW.out.chr_vcf)
    
    PROCESS_GERMLINE_PREPROCESS(process_germline_jobs_ch)
    processed_output = PROCESS_GERMLINE_PREPROCESS.out.rds
    r_output = GERMLINE_PREPROCESS.out.rds 
    
    // Ensure both sides are done before generating report
    normal_sample_ch
        .map { row -> row.spn }
        .distinct()
        .set { spn_ch }
    
    //GENERATE_GERMLINE_REPORT(spn_ch)
    //report_germline = GENERATE_GERMLINE_REPORT.out.metrics_germline_report
    emit:
    // viewed_vcfs
    r_output
}
