nextflow.enable.dsl=2

include { BCFTOOLS_VIEW } from '../../modules/bcftools/view/main.nf'
include { GERMLINE_PREPROCESS } from '../../modules/preprocessing/germline_callers/main.nf'
include { PROCESS_GERMLINE_PREPROCESS } from '../../modules/preprocessing/process/normal/main.nf'
include { GENERATE_GERMLINE_REPORT }    from '../../modules/germline_report/main.nf'

workflow GERMLINE_VALIDATION {

    take:
    normal_sample_ch

    main:

    normal_sample_ch
        .flatMap { meta,path ->
            def chromosomes = (20..22).collect { it.toString() } + ['X']
            if (meta.sex == 'XY') {
              chromosomes += 'Y' 
            }
            def callers = ['haplotypecaller', 'freebayes', 'strelka']
            def status = "germline"
            chromosomes.collectMany { chr ->
                callers.collect { caller ->
                    [meta, path, chr, caller, status]
                }
            }
        }.map { meta, path, chr, caller, status ->
            def vcf_path
            def mut_type = "all"
            if (caller == 'strelka') {
                vcf_path = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${caller}/${meta.sample}/${meta.sample}.${caller}.variants.vcf.gz"
                vcf_tbi_path = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${caller}/${meta.sample}/${meta.sample}.${caller}.variants.vcf.gz.tbi"
            } else if (caller == 'haplotypecaller') {
                vcf_path = "${path}/${meta.spn}/sarek/normal/variant_calling/${caller}/${meta.sample}/${meta.sample}.${caller}.filtered.vcf.gz"
                vcf_tbi_path = "${path}/${meta.spn}/sarek/normal/variant_calling/${caller}/${meta.sample}/${meta.sample}.${caller}.filtered.vcf.gz.tbi"
            } else if (caller == 'freebayes') {
                vcf_path = "${path}/${meta.spn}/sarek/normal/variant_calling/${caller}/${meta.sample}/${meta.sample}.${caller}.vcf.gz"
                vcf_tbi_path = "${path}/${meta.spn}/sarek/normal/variant_calling/${caller}/${meta.sample}/${meta.sample}.${caller}.vcf.gz.tbi"
            }
            return [meta, chr, caller, status, mut_type, file(vcf_path), file(vcf_tbi_path)]
        }
        .set { normal_jobs_ch }

    normal_sample_ch
        .flatMap { meta, path ->
            def callers = ['process']
            def status = "germline"
            callers.collect { caller ->
                [meta, path, caller, status]
            }
        }
        .map { meta, path, caller, status ->
            def rds_path = "${path}/${meta.spn}/sequencing/normal/purity_1/data/mutations/seq_results_muts_merged_coverage_30x.rds"
            return [meta, caller, file(rds_path)]
        }
        .set { process_germline_jobs_ch }

    BCFTOOLS_VIEW(normal_jobs_ch) // split vcf
    GERMLINE_PREPROCESS(BCFTOOLS_VIEW.out.chr_vcf) // from vcf to rds
   
    GERMLINE_PREPROCESS.out.rds.map{meta, chr, caller, rds -> 
            meta = meta + [caller: caller]
            [meta, rds]
           }
           | groupTuple
           | map{ meta, rds -> 
            caller = meta.caller
            [meta.subMap('spn', 'sample', 'coverage', 'purity', 'type', 'sex'), rds, caller]
           }
           | groupTuple
           | set { vcf_caller_join }
     
    PROCESS_GERMLINE_PREPROCESS(process_germline_jobs_ch) // read process input, 1 unique rds
    PROCESS_GERMLINE_PREPROCESS.out.rds.map{
        meta, type, rds -> [meta, rds]
    }.set{ process_rds }
    
    
    GENERATE_GERMLINE_REPORT(vcf_caller_join.join(process_rds))
    report_germline = GENERATE_GERMLINE_REPORT.out.metrics_germline_report
    
    emit:
    report_germline
    // viewed_vcfs
    // r_output
}
