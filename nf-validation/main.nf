nextflow.enable.dsl=2
include { samplesheetToList } from 'plugin/nf-schema'

include { GERMLINE_VALIDATION } from './subworkflows/germline_validation/main.nf'
include { SOMATIC_VALIDATION } from './subworkflows/somatic_validation/main.nf'

workflow {

    //samplesheet = params.input ?: 'samples.csv'
    samplesheet = params.input ? Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json")) : Channel.empty()
    
    // Load samplesheet into channel
    samples_ch_normal = samplesheet.filter { meta,path -> meta.type == 0 }
    samples_ch_tumour = samplesheet.filter { meta,path -> meta.type == 1 }

    // Trigger germline subworkflow
    //GERMLINE_VALIDATION(samples_ch_normal)
    SOMATIC_VALIDATION(samples_ch_tumour)
    
}
