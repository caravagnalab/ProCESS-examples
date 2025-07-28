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
        .flatMap { row ->
            def chromosomes = (1..22).collect { it.toString() } + ['X']
            if (row.sex == 'XY') {
                chromosomes += 'Y'
            }
            def callers = ['mutect2', 'strelka','freebayes']
            def status = "somatic"
            chromosomes.collectMany { chr ->
                callers.collect { caller ->
                    [row, chr, caller, status]
                }
            }
        }
        .flatMap { row, chr, caller, status ->
    
            if (caller == 'mutect2') {
                def vcf_path = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.spn}/${row.spn}.mutect2.filtered.vcf.gz"
                def vcf_tbi_path = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.spn}/${row.spn}.mutect2.filtered.vcf.gz.tbi"
                def mut_type = "all"
                return [[row, chr, caller, status, mut_type, file(vcf_path), file(vcf_tbi_path)]]
            } else if (caller == 'strelka') { //SPN01_1.3_vs_normal_sample
                def snv_vcf = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.sample}_vs_normal_sample/${row.sample}_vs_normal_sample.strelka.somatic_snvs.vcf.gz"
                def snv_tbi = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.sample}_vs_normal_sample/${row.sample}_vs_normal_sample.strelka.somatic_snvs.vcf.gz.tbi"
                def indel_vcf = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.sample}_vs_normal_sample/${row.sample}_vs_normal_sample.strelka.somatic_indels.vcf.gz"
                def indel_tbi = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.sample}_vs_normal_sample/${row.sample}_vs_normal_sample.strelka.somatic_indels.vcf.gz.tbi"
                def mut_type_indels = "INDEL"
                def mut_type_snv = "SNV" 
                return [
                    [row, chr, caller, status, mut_type_snv, file(snv_vcf), file(snv_tbi)],
                    [row, chr, caller, status, mut_type_indels,file(indel_vcf), file(indel_tbi)]
                ]
            } else if (caller=="freebayes") { //freebayes/SPN01_1.2_vs_normal_sample/SPN01_1.2_vs_normal_sample.freebayes.vcf.gz
                def vcf_path = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.sample}_vs_normal_sample/${row.sample}_vs_normal_sample.freebayes.vcf.gz"
                def vcf_tbi_path = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.sample}_vs_normal_sample/${row.sample}_vs_normal_sample.freebayes.vcf.gz.tbi"
                def mut_type = "all"
                return [[row, chr, caller, status, mut_type, file(vcf_path), file(vcf_tbi_path)]]
            }
        }
        .set { somatic_jobs_ch }
    
//    t_sample_ch
//        .flatMap { row ->
//            def chromosomes = (1..22).collect { it.toString() } + ['X']
//            if (row.sex == 'XY') {
//                chromosomes += 'Y'
//            }
//            def callers = ['mutect2'] // , 'strelka']
//            def status = "somatic"
//            chromosomes.collectMany { chr ->
//                callers.collect { caller ->
//                    [row, chr, caller,status]
//                }
//            }
//
//        }
//        .map { row, chr, caller, status ->
//            def vcf_path
//            if (caller == 'mutect2') {
//                vcf_path = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.spn}/${row.spn}.mutect2.filtered.vcf.gz"
//                vcf_tbi_path = "${row.directory}/${row.spn}/sarek/${row.coverage}x_${row.purity}p/variant_calling/${caller}/${row.spn}/${row.spn}.mutect2.filtered.vcf.gz.tbi"
//            } else if (caller == 'strelka') {
//                //vcf_path = "${row.directory}/${row.spn}/sarek/normal_sample/variant_calling/${caller}/normal_sample/normal_sample.freebayes.vcf.gz"
//            }
//            return [row, chr, caller, status,file(vcf_path), file(vcf_tbi_path)]
//        }
//        .set { somatic_jobs_ch }
//
    t_sample_ch
        .flatMap { row ->
            def chromosomes = (1..22).collect { it.toString() } + ['X']
            if (row.sex == 'XY') {
                chromosomes += 'Y'
            }
            def callers = ['process']
            def status = "somatic"
            chromosomes.collectMany { chr ->
                callers.collect { caller ->
                    [row, chr, caller,status]
                }
            }

        }
        .map { row, chr, caller, status ->
            rds_path = "${row.directory}/${row.spn}/sequencing/tumour/purity_${row.purity}/data/mutations/seq_results_muts_merged_coverage_${row.coverage}x.rds"
            return [row, chr, caller, file(rds_path)]
        }
        .set { process_somatic_jobs_ch }


    //viewed_vcfs=BCFTOOLS_VIEW(somatic_jobs_ch)

    BCFTOOLS_VIEW(somatic_jobs_ch)
    SOMATIC_PREPROCESS(BCFTOOLS_VIEW.out.chr_vcf)

    PROCESS_PREPROCESS(process_somatic_jobs_ch)
    
    processed_output = PROCESS_PREPROCESS.out.rds
    r_output = SOMATIC_PREPROCESS.out.rds
    
    
    // Prepare input for final report
    t_sample_ch
        .map { row -> tuple(row.spn, row.coverage, row.purity) }
        .distinct()
        .set { spn_coverage_purity_ch }
        
    GENERATE_SOMATIC_REPORT_SINGLE_COMBINATION(spn_coverage_purity_ch)
    report_somatic_single_combination = GENERATE_SOMATIC_REPORT_SINGLE_COMBINATION.out.metrics_somatic_report

    GENERATE_SOMATIC_REPORT_ALL_CALLER(spn_coverage_purity_ch)
    report_somatic_all_callers = GENERATE_SOMATIC_REPORT_ALL_CALLER.out.metrics_somatic_all_caller_report
    
    emit:
    //viewed_vcfs
    report_somatic_single_combination
    report_somatic_all_callers
}
