nextflow.enable.dsl=2

include { CNA_REPORT_COMBINATION } from '../../modules/cna_report_combination/main.nf'
include { CNA_REPORT_ALL } from '../../modules/cna_report_all/main.nf'

workflow CNA_VALIDATION {

    take:
    t_sample_ch

    main:

    CNA_REPORT_COMBINATION(t_sample_ch)

    CNA_REPORT_COMBINATION.out.rds_metric.map{meta, rds ->
        def sample=meta.sample
        [meta.subMap('spn', 'coverage','purity', 'type', 'sex'), rds, sample]
    }
    |groupTuple
    |map{ meta, rds, sample -> 
        def combination="${meta.coverage}_${meta.purity}"
        [meta.subMap('spn', 'type', 'sex'), rds, sample, combination]
    }
    |groupTuple
    |set{spn_ch}

    CNA_REPORT_ALL(spn_ch)

    emit:
    rds_all = CNA_REPORT_ALL.out.rds
    plot_all = CNA_REPORT_ALL.out.plot
}
