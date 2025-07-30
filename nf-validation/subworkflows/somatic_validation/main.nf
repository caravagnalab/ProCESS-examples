nextflow.enable.dsl=2

include { BCFTOOLS_CONCAT } from '../../modules/bcftools/concat/main.nf'
include { BCFTOOLS_VIEW } from '../../modules/bcftools/view/main.nf'
include { SOMATIC_PREPROCESS } from '../../modules/preprocessing/callers/main.nf'
include { PROCESS_PREPROCESS } from '../../modules/preprocessing/process/main.nf'
include { GENERATE_SOMATIC_REPORT_SINGLE_COMBINATION } from '../../modules/somatic_report_single_combination/main.nf'
include { GENERATE_SOMATIC_REPORT_ALL_CALLER } from '../../modules/somatic_report_all_caller/main.nf'
include { GENERATE_SOMATIC_REPORT_ALL_COMBINATION } from '../../modules/somatic_report_all_combination/main.nf'


workflow SOMATIC_VALIDATION {
    take:
    t_sample_ch

    main:
    t_sample_ch.map{ meta, path -> 
        meta = meta + [caller:'strelka']
        def snv_vcf = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${meta.caller}/${meta.sample}_vs_normal_sample/${meta.sample}_vs_normal_sample.strelka.somatic_snvs.vcf.gz"
        def snv_tbi = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${meta.caller}/${meta.sample}_vs_normal_sample/${meta.sample}_vs_normal_sample.strelka.somatic_snvs.vcf.gz.tbi"
        def indel_vcf = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${meta.caller}/${meta.sample}_vs_normal_sample/${meta.sample}_vs_normal_sample.strelka.somatic_indels.vcf.gz"
        def indel_tbi = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${meta.caller}/${meta.sample}_vs_normal_sample/${meta.sample}_vs_normal_sample.strelka.somatic_indels.vcf.gz.tbi"
        [meta, snv_vcf, snv_tbi, indel_vcf, indel_tbi]
    }.set{ vcf_strelka }

    BCFTOOLS_CONCAT(vcf_strelka).flatMap { meta, vcf ->
        def chromosomes = (20..22).collect { it.toString() } + ['X']
            if (meta.sex == 'XY') {
                chromosomes += 'Y'
            }
            chromosomes.collect { chr ->
             [meta, chr, vcf]
            }
        }.set {strelka_somatic_ch}


    t_sample_ch
        .flatMap { meta, path ->
            def chromosomes = (20..22).collect { it.toString() } + ['X']
            if (meta.sex == 'XY') {
                chromosomes += 'Y'
            }
            def callers = ['freebayes', 'mutect2']
            chromosomes.collectMany { chr ->
                callers.collect { caller ->
                    [meta, path, chr, caller]
                }
            }
        }
        .flatMap { meta, path, chr, caller ->
            if (caller == 'mutect2') {
                meta = meta + [caller:'mutect2']
                def vcf_path = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${meta.caller}/${meta.spn}/${meta.spn}.mutect2.filtered.vcf.gz"
                return [[meta, chr, file(vcf_path)]]

            } else if (caller=="freebayes") { 
                meta = meta + [caller:'freebayes']
                def vcf_path = "${path}/${meta.spn}/sarek/${meta.coverage}x_${meta.purity}p/variant_calling/${meta.caller}/${meta.sample}_vs_normal_sample/${meta.sample}_vs_normal_sample.freebayes.vcf.gz"
                return [[meta, chr, file(vcf_path)]]
            }
        }
        .set { mutect_freeb_somatic_ch }

    somatic_ch = strelka_somatic_ch.concat(mutect_freeb_somatic_ch)
    BCFTOOLS_VIEW(somatic_ch)
    SOMATIC_PREPROCESS(BCFTOOLS_VIEW.out.chr_vcf)

    t_sample_ch.flatMap { meta, path ->
            def chromosomes = (20..22).collect { it.toString() } + ['X']
            if (meta.sex == 'XY') {
                chromosomes += 'Y'
            }
            def callers = ['process']
            chromosomes.collectMany { chr ->
                callers.collect { caller ->
                    [meta, path, chr, caller]
                }
            }
        }
        .map { meta, path, chr, caller ->
            meta = meta + [caller:caller]
            def rds_path = "${path}/${meta.spn}/sequencing/tumour/purity_${meta.purity}/data/mutations/seq_results_muts_merged_coverage_${meta.coverage}x.rds"
            return [meta, chr, file(rds_path)]
        }
        .set { process_somatic_ch }
    PROCESS_PREPROCESS(process_somatic_ch)

    SOMATIC_PREPROCESS.out.rds.map { meta, chr, rds -> 
        [meta, rds] }
        | groupTuple
        | map{ meta, rds -> 
            [meta.subMap('spn', 'coverage', 'purity', 'type', 'sex', 'caller'), rds, meta.sample]} 
        | groupTuple
        | map {meta, rds, sample -> 
            def caller = meta.caller
            [meta.subMap('spn', 'coverage', 'purity', 'type', 'sex'), rds, sample, caller]
        }
        | set {spn_coverage_purity_ch}

    PROCESS_PREPROCESS.out.rds.map { meta, chr, rds ->
        [meta, rds]}
        | groupTuple
        | map{ meta, rds -> 
            [meta.subMap('spn', 'coverage', 'purity', 'type', 'sex', 'caller'), rds, meta.sample]} 
        | groupTuple
        | map {meta, rds, sample -> 
            def caller = meta.caller
            [meta.subMap('spn', 'coverage', 'purity', 'type', 'sex'), rds, sample, caller]
        }
        | set { process_coverage_purity_ch }

    spn_coverage_purity_ch
    | groupTuple
    | set { spn_coverage_purity_ch_caller }

    GENERATE_SOMATIC_REPORT_SINGLE_COMBINATION(process_coverage_purity_ch.combine(spn_coverage_purity_ch, by:0))
    GENERATE_SOMATIC_REPORT_ALL_CALLER(process_coverage_purity_ch.combine(spn_coverage_purity_ch_caller, by:0))

    report_comb_caller = GENERATE_SOMATIC_REPORT_SINGLE_COMBINATION.out.metrics_somatic_report
    report_comb = GENERATE_SOMATIC_REPORT_ALL_CALLER.out.metrics_somatic_all_caller_report

    metric_comb_caller = GENERATE_SOMATIC_REPORT_SINGLE_COMBINATION.out.metrics_somatic_rds
    metric_comb = GENERATE_SOMATIC_REPORT_ALL_CALLER.out.metrics_somatic_all_caller_rds
    
    metric_comb.map{ meta, rds -> 
        def combination = "${meta.coverage}_${meta.purity}"
        [meta.subMap('spn', 'type', 'sex'), rds, combination]
    }
    | groupTuple
    | set { spn_ch }

    GENERATE_SOMATIC_REPORT_ALL_COMBINATION(spn_ch)
    plot_spn = GENERATE_SOMATIC_REPORT_ALL_COMBINATION.out.report

    emit:
    metric_comb_caller
    report_comb_caller
    metric_comb 
    report_comb 
    plot_spn
}
