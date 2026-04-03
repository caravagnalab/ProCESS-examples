nextflow.enable.dsl=2
include { CNAQC_VALIDATION_COMBINATION } from '../../modules/cnaqc_validation/main.nf'

workflow QC_VALIDATION {

    take:
    t_sample_ch

    main:

    t_sample_ch.map{ meta, path -> 
        [meta.subMap('spn', 'coverage', 'purity', 'sex', 'type'), path]
    }
    .unique()
    .flatMap { meta, path ->
        def vcf_callers = params.variant_caller.split(',') 
        def cna_caller = params.cna_caller.split(',') 

        vcf_callers.collectMany { vcf ->
            cna_caller.collect { cna ->
                [meta, path, vcf, cna]
            }
        }
    }.set{ signature_ch }

    CNAQC_VALIDATION_COMBINATION(signature_ch)
    cnaqc_validate = CNAQC_VALIDATION_COMBINATION.out.comb_cnaqc_validate
    cnaqc_stats = CNAQC_VALIDATION_COMBINATION.out.comb_cnaqc_stats
    cnaqc_report = CNAQC_VALIDATION_COMBINATION.out.comb_cnaqc_validate_report
    

    emit:
    cnaqc_validate
    cnaqc_stats
    cnaqc_report
}
