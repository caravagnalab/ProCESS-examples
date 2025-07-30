nextflow.enable.dsl=2
include { samplesheetToList } from 'plugin/nf-schema'

include { GERMLINE_VALIDATION } from './subworkflows/germline_validation/main.nf'
include { SOMATIC_VALIDATION } from './subworkflows/somatic_validation/main.nf'
include { CNA_VALIDATION } from './subworkflows/cna_validation/main.nf'
include { SIGNATURE_VALIDATION } from './subworkflows/signature_validation/main.nf'


workflow {
    samplesheet = params.input ? Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json")) : Channel.empty()
    
    samples_ch_normal = samplesheet.filter { meta, path -> meta.type == 0 }
    samples_ch_tumour = samplesheet.filter { meta, path -> meta.type == 1 }

    // SAREK
    if (params.step && params.step.split(',').contains('germline')){
        GERMLINE_VALIDATION(samples_ch_normal)
    } 
    
    if (params.step && params.step.split(',').contains('somatic')){
        SOMATIC_VALIDATION(samples_ch_tumour)
    } 

    if (params.step && params.step.split(',').contains('cna')){
        CNA_VALIDATION(samples_ch_tumour)
    }

    // TUMOUR
    //SIGNATURE_VALIDATION(samples_ch_tumour)
    //DRIVER_VALIDATION()
    //SUBCLONAL_VALIDATION()
    //QC_VALIDATION()

}
