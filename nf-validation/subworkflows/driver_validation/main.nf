nextflow.enable.dsl=2
include { DRIVER_VALIDATION_COMBINATION } from '../../modules/driver_validation/main.nf'

workflow DRIVER_VALIDATION {
  
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
  
  DRIVER_VALIDATION_COMBINATION(signature_ch)
  confusion_matrix = DRIVER_VALIDATION_COMBINATION.out.confusion_matrix
  report = DRIVER_VALIDATION_COMBINATION.out.report
  heatmap = DRIVER_VALIDATION_COMBINATION.out.heatmap
  driver_comparsion = DRIVER_VALIDATION_COMBINATION.out.driver_comparsion
  
  
  
  emit:
    confusion_matrix
    driver_comparsion
    report
    heatmap
}
