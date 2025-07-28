nextflow.enable.dsl=2

include { BCFTOOLS_VIEW } from '../../modules/bcftools/view/main.nf'
include { SOMATIC_PREPROCESS } from '../../modules/preprocessing/callers/main.nf'
include { PROCESS_PREPROCESS } from '../../modules/preprocessing/process/main.nf'
include { GENERATE_SOMATIC_REPORT_SINGLE_COMBINATION } from '../../modules/somatic_report_single_combination/main.nf'
include { GENERATE_SOMATIC_REPORT_ALL_CALLER } from '../../modules/somatic_report_all_caller/main.nf'


workflow SOMATIC_VALIDATION {

    take:
    t_sample_ch

    main:
    t_sample_ch
        .flatMap { meta,path ->
            def chromosomes = (20..22).collect { it.toString() } + ['X']
            if (meta.sex == 'XY') {
                chromosomes += 'Y'
            }
            def callers = ['mutect2', 'strelka','freebayes']
            def status = "somatic"
            chromosomes.collectMany { chr ->
                callers.collect { caller ->
                    [meta, path, chr, caller, status]
                }
            }
        }
        .flatMap { meta, path, chr, caller, status ->
            if (caller == 'mutect2') {
                def vcf_path = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${caller}/${meta.spn}/${meta.spn}.mutect2.filtered.vcf.gz"
                def vcf_tbi_path = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${caller}/${meta.spn}/${meta.spn}.mutect2.filtered.vcf.gz.tbi"
                def mut_type = "all"
                return [[meta, chr, caller, status, mut_type, file(vcf_path), file(vcf_tbi_path)]]
            } else if (caller == 'strelka') { //SPN01_1.3_vs_normal_sample
                def snv_vcf = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${caller}/${meta.sample}_vs_normal_sample/${meta.sample}_vs_normal_sample.strelka.somatic_snvs.vcf.gz"
                def snv_tbi = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${caller}/${meta.sample}_vs_normal_sample/${meta.sample}_vs_normal_sample.strelka.somatic_snvs.vcf.gz.tbi"
                def indel_vcf = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${caller}/${meta.sample}_vs_normal_sample/${meta.sample}_vs_normal_sample.strelka.somatic_indels.vcf.gz"
                def indel_tbi = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${caller}/${meta.sample}_vs_normal_sample/${meta.sample}_vs_normal_sample.strelka.somatic_indels.vcf.gz.tbi"
                def mut_type_indels = "INDEL"
                def mut_type_snv = "SNV" 
                return [
                    [meta, chr, caller, status, mut_type_snv, file(snv_vcf), file(snv_tbi)],
                    [meta, chr, caller, status, mut_type_indels, file(indel_vcf), file(indel_tbi)]
                ]
            } else if (caller=="freebayes") { //freebayes/SPN01_1.2_vs_normal_sample/SPN01_1.2_vs_normal_sample.freebayes.vcf.gz
                def vcf_path = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${caller}/${meta.sample}_vs_normal_sample/${meta.sample}_vs_normal_sample.freebayes.vcf.gz"
                def vcf_tbi_path = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${caller}/${meta.sample}_vs_normal_sample/${meta.sample}_vs_normal_sample.freebayes.vcf.gz.tbi"
                def mut_type = "all"
                return [[meta, chr, caller, status, mut_type, file(vcf_path), file(vcf_tbi_path)]]
            }
        }
        .set { somatic_jobs_ch }

    // t_sample_ch
    //     .flatMap { meta ->
    //         def chromosomes = (1..22).collect { it.toString() } + ['X']
    //         if (meta.sex == 'XY') {
    //             chromosomes += 'Y'
    //         }
    //         def callers = ['process']
    //         def status = "somatic"
    //         chromosomes.collectMany { chr ->
    //             callers.collect { caller ->
    //                 [meta, chr, caller,status]
    //             }
    //         }

    //     }
    //     .map { meta, chr, caller, status ->
    //         rds_path = "${path}/${meta.spn}/sequencing/tumour/purity_${meta.purity}/data/mutations/seq_results_muts_merged_coverage_${meta.coverage}x.rds"
    //         return [meta, chr, caller, file(rds_path)]
    //     }
    //     .set { process_somatic_jobs_ch }

    BCFTOOLS_VIEW(somatic_jobs_ch)
    // SOMATIC_PREPROCESS(BCFTOOLS_VIEW.out.chr_vcf)

    // PROCESS_PREPROCESS(process_somatic_jobs_ch)
    
    // processed_output = PROCESS_PREPROCESS.out.rds
    // r_output = SOMATIC_PREPROCESS.out.rds
    
    
    // Prepare input for final report
    // t_sample_ch
    //     .map { meta -> tuple(meta.spn, meta.coverage, meta.purity) }
    //     .distinct()
    //     .set { spn_coverage_purity_ch }
        
    // GENERATE_SOMATIC_REPORT_SINGLE_COMBINATION(spn_coverage_purity_ch)
    // report_somatic_single_combination = GENERATE_SOMATIC_REPORT_SINGLE_COMBINATION.out.metrics_somatic_report

    // GENERATE_SOMATIC_REPORT_ALL_CALLER(spn_coverage_purity_ch)
    // report_somatic_all_callers = GENERATE_SOMATIC_REPORT_ALL_CALLER.out.metrics_somatic_all_caller_report
    
    emit:
    null
    //viewed_vcfs
    // report_somatic_single_combination
    // report_somatic_all_callers
}
