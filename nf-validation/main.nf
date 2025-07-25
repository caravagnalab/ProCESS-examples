nextflow.enable.dsl=2

include { GERMLINE_VALIDATION } from './subworkflows/germline_validation/main.nf'
include { SOMATIC_VALIDATION } from './subworkflows/somatic_validation/main.nf'

workflow {

    samplesheet = params.input ?: 'samples.csv'

    // Load samplesheet into channel
    samples_ch_normal = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .filter { row -> row.type == '0' } // select normal_sample (status==0)

    samples_ch_tumour = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .filter { row -> row.type == '1' } // select normal_sample (status==0)
    // Trigger germline subworkflow
    GERMLINE_VALIDATION(samples_ch_normal)
    SOMATIC_VALIDATION(samples_ch_tumour)
}
