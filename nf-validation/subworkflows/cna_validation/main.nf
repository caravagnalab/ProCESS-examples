nextflow.enable.dsl=2

// include { BCFTOOLS_VIEW } from '../../modules/bcftools/view/main.nf'
// include { GERMLINE_PREPROCESS } from '../../modules/preprocessing/germline_callers/main.nf'
// include { PROCESS_GERMLINE_PREPROCESS } from '../../modules/preprocessing/process/normal/main.nf'
// include { GENERATE_GERMLINE_REPORT }    from '../../modules/germline_report/main.nf'

workflow CNA_VALIDATION {

    take:
    tumour_sample_ch

    main:
    tumour_sample_ch.view()

    emit:
    null
}
