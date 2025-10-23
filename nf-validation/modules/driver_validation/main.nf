process DRIVER_VALIDATION_COMBINATION {
    tag "${meta.spn}_${meta.coverage}_${meta.purity}-driver_validation"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v2' :
        'docker.io/lvaleriani/process_validation:v2' }"

    input:
    tuple val(meta), val(path), val(vcf_caller), val(cna_caller)

    output:
    tuple val(meta), path("confusion_matrix.rds"),             emit:confusion_matrix
    tuple val(meta), path("all_driver_comparsion.rds"),        emit:driver_comparsion
    tuple val(meta), path("real_drivers_vaf_comparison.rds"),          emit:heatmap
    tuple val(meta), path("report.png"),	            emit:report

    
    publishDir "${params.outdir}/${meta.spn}/driver/${meta.coverage}x_${meta.purity}/${vcf_caller}_${cna_caller}", mode: 'copy'

    script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    library(patchwork)
    library(ProCESS)
    library(ComplexHeatmap)
    library(circlize)
    library(viridis)
    library(awtools)
    library(cowplot)
    



    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/getters/tumourevo_getters.R")
    source("${projectDir}/bin/drivers/drivers_validation_functions.R")




    data_dir = "${path}"
    spn_id = "${meta.spn}"
    coverage = as.numeric("${meta.coverage}")
    purity = as.numeric("${meta.purity}")

    vcf_caller <- "${vcf_caller}"
    cna_caller <- "${cna_caller}"
    callers <- paste0(vcf_caller,"_",cna_caller)
    
    samples = get_sample_names(spn = spn_id)
    
    # load the phylogenetic forest
    phylo_forest = load_phylogenetic_forest(get_phylo_forest(spn_id, data_dir))

    tumourevo_mutations = get_drivers_results(spn = spn_id, 
                                              samples = samples, 
                                              purity = purity, 
                                              coverage = coverage, 
                                              callers = callers, 
                                              cohort = 'SCOUT', 
                                              path = data_dir)
    
    
    # load process data and select its drivers
    process_seq_res = get_process_seq(spn = spn_id, 
                                      purity = purity, 
                                      coverage = coverage,
                                      path = data_dir)
                                      
                                      
    # load process data and select its drivers
    process_seq_res = get_process_seq(spn = spn_id, 
                                      purity = purity, 
                                      coverage = coverage,
                                      path = data_dir)
                                      

    process_drivers_ids = get_process_drivers_codes(phylo_forest)
    
    # filter the sequencing results to keep only driver mutations
    process_drivers = get_process_drivers(process_seq_res, phylo_forest)
    
    true_drivers_table = create_true_drivers_table(process_drivers,
                                                   tumourevo_mutations, 
                                                   process_drivers_ids)
    
    # plot the heatmap
    ht = plot_drivers_heatmap(true_drivers_table)
    saveRDS(object = ht,file = paste0("real_drivers_vaf_comparison.rds"))
    
    
    
    process_drivers_ids = get_process_drivers_codes(phylo_forest)
    
    # get drivers codes in tumourevo

    tumourevo_drivers_ids = get_tumourevo_drivers_codes(tumourevo_mutations)
    
    all_drivers = c(process_drivers_ids, tumourevo_drivers_ids) %>% unique
    
    process_drivers = get_process_drivers(process_seq_res, phylo_forest) 
    process_all_drs = get_all_drivers_process(process_drivers, all_drivers)
    
    tumourevo_all_drs = get_all_drivers_tumourevo(tumourevo_mutations, all_drivers)
    
    all_drivers_table = merge_drivers(all_dr_process = process_all_drs,all_drivers_tumourevo =  tumourevo_all_drs)
    metrics_df <- all_drivers_table %>%
      group_by(sample) %>%
      summarise(
        TP = sum(driver_class == "Process True - Tumourevo True"),
        FP = sum(driver_class == "Process False - Tumourevo True"),
        FN = sum(driver_class == "Process True - Tumourevo False"),
        .groups = "drop"
      ) %>%
      mutate(
        precision = ifelse(TP + FP > 0, TP / (TP + FP), NA),
        recall    = ifelse(TP + FN > 0, TP / (TP + FN), NA)
      ) %>%
      select(sample, precision, recall) %>% 
      mutate(coverage=coverage) %>% 
      mutate(purity=purity) %>% 
      mutate(spn=spn_id)
    plt = plot_drivers(all_drivers_table, colors = colors)
    
    ggsave(filename = "report.png",plot = plt,width = 10,height = 5)
    
    saveRDS(object=metrics_df,file="confusion_matrix.rds")
    saveRDS(object=all_drivers_table,file="all_driver_comparsion.rds")

    """
}