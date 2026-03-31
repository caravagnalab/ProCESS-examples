nextflow.enable.dsl=2
include { SIGNATURE_VALIDATION_COMBINATION } from '../../modules/signature_validation/main.nf'

workflow SIGNATURE_VALIDATION {

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

    SIGNATURE_VALIDATION_COMBINATION(signature_ch)
    metrics_spn = SIGNATURE_VALIDATION_COMBINATION.out.metrics_spn
    metrics_sample = SIGNATURE_VALIDATION_COMBINATION.out.metrics_sample
    cosine_mse = SIGNATURE_VALIDATION_COMBINATION.out.cosine_mse
    sankey_plot_png = SIGNATURE_VALIDATION_COMBINATION.out.sankey_plot_png
    sankey_plot_rds = SIGNATURE_VALIDATION_COMBINATION.out.sankey_plot_rds
    

    emit:
    cosine_mse
    metrics_spn
    metrics_sample
    sankey_plot_png
    sankey_plot_rds
}
