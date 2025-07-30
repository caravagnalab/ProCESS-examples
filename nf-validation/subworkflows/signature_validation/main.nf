nextflow.enable.dsl=2

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
    }.view()

    emit:
    null
}
