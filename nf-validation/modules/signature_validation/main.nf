process SIGNATURE_VALIDATION_COMBINATION {
    tag "${meta.spn}_${meta.coverage}_${meta.purity}-signature_validation"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://lvaleriani/process_validation:v1' :
        'docker.io/lvaleriani/process_validation:v1' }"

    input:
    tuple val(meta), val(path), val(vcf_caller), val(cna_caller)

    output:
    tuple val(meta), path("metrics_*_spn.rds"),             emit:metrics_spn
    tuple val(meta), path("metrics_*_sample.rds"),          emit:metrics_sample
    tuple val(meta), path("cosine_mse_*.rds"),	            emit:cosine_mse
    tuple val(meta), path("sankey_plot_*.png"),             emit:sankey_plot

    
    publishDir "${params.outdir}/${meta.spn}/signature/${meta.coverage}x_${meta.purity}/${vcf_caller}_${cna_caller}", mode: 'copy'

    script:
    """
    #!/usr/bin/env Rscript
    library(dplyr)
    library(ProCESS)
    library(tidyverse)
    library(lsa)
    library(Metrics)
    library(ggplot2)
    library(ggalluvial)


    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/getters/tumourevo_getters.R")
    source("${projectDir}/bin/signatures/utils_getters.R")
    source("${projectDir}/bin/signatures/utils_validation.R")
    source("${projectDir}/bin/signatures/utils.R")
    source("${projectDir}/bin/signatures/utils_plots.R")
    source("${projectDir}/bin/signatures/utils_sparsesig.R")




    data_dir = "${path}"
    spn_id = "${meta.spn}"
    coverage = as.numeric("${meta.coverage}")
    purity = as.numeric("${meta.purity}")

    vcf_caller <- "${vcf_caller}"
    cna_caller <- "${cna_caller}"
    contexts_all <- c("SBS96","ID83")
    context_classes <-  gsub('[[:digit:]]+', '', contexts_all) 

    #### COSMIC path
    cosmic_path <- '${projectDir}/bin/data/COSMIC_v3.4_SBS_GRCh38.txt'

    message("Reading signature exposures from ProCESS")
    gt_exposure <-  get_process_exposures(spn = spn_id, 
                                        coverage = coverage,
                                        purity = purity)

    message("Reading inferred signatures")                                    
    all_exposures <- list()
    metrics_sample_all <- list()
    metrics_spn_all <- list()
    cosine_mse_all <- list()
    for (context in contexts_all){
        context_classes <- gsub('[[:digit:]]+', '', context)
        ground_truth_nested <- gt_exposure[[context_classes]]  %>%
            tibble::column_to_rownames("Sample_ID") %>% 
            as.matrix()
        ground_truth_nested[is.na(ground_truth_nested)] <- 0
        
        
        result <- tryCatch({
            # Get SparseSignatures paths
            # sparsesig <- get_tumourevo_signatures(
            #    spn = spn_id,
            #    coverage = coverage,
            #    purity = purity,
            #    vcf_caller = vcf_caller,
            #    cna_caller = cna_caller,
            #    tool = "SparseSignatures"
            #)
            
            # Get SigProfiler paths
            sigprofiler <- get_tumourevo_signatures(
                spn = spn_id,
                coverage = coverage,
                purity = purity,
                vcf_caller = vcf_caller,
                cna_caller = cna_caller,
                tool = "SigProfiler",
                context = context
            )
            
            bascule <- get_tumourevo_signatures(
                spn = spn_id,
                coverage = coverage,
                purity = purity,
                vcf_caller = vcf_caller,
                cna_caller = cna_caller,
                tool = "BASCULE",
                context = context_classes
            )
            # Combine and load paths
            paths <- c(
                #sparsesig\$nmf_Lasso_out,
                #sparsesig\$cv_means_mse,
                #sparsesig\$best_params_config,
                #sparsesig\$mut_counts,
                sigprofiler\$COSMIC_exposure,
                sigprofiler\$COSMIC_signatures,
                sigprofiler\$denovo_exposure,
                sigprofiler\$denovo_signatures,
                bascule\$refined_fit
            )
            
            data <- load_signature_data(paths)
            
            names(data) <- c(
            #"SparseSig_nmf_Lasso_out",
            #"SparseSig_cv_means_mse",
            #"SparseSig_best_params_config",
            #"SparseSig_mut_counts",
            "SigProfiler_COSMIC_exposure",
            "SigProfiler_COSMIC_signatures",
            "SigProfiler_denovo_exposure",
            "SigProfiler_denovo_signatures",
            "BASCULE_refined_fit"
            )
            data
            
        })
            
        if (!is.null(result)) {
            tumourevo_signature_res <- result
        }
            
        # sig_out <- tumourevo_signature_res[["SparseSig_nmf_Lasso_out"]]
        # if (is.null(sig_out)) {
        #   warning(paste("Missing sig_out for", spn_id, coverage, purity))
        #   mapped_res <- NULL
        #   next
        # }
        # if (context=="SBS96"){
        #   remapped <- tryCatch({
        #     map_sparsesig_to_cosmic(
        #       sparsesig_out = sig_out,
        #       cosmic_path = cosmic_path,
        #       threshold = 0.5
        #     )
        #   }, error = function(e) {
        #     message(paste("Error in", spn_id, coverage, purity, ":", e\$message))
        #     return(NULL)
        #   })
        #   mapped_res <- remapped
        #   sparsig_cosmic <- mapped_res
        #   sparsesig_aligned <- align_callers(tumourevo_signature_res = sparsig_cosmic,tool = "SparseSignatures",spn = spn_id) %>% as.data.frame()
        #   sparsesig_long <- reshape_exposures_long(exposures_mat=sparsesig_aligned,spn = spn_id,
        #                                            coverage = coverage,purity = purity,method_name = "SparseSignatures") %>% as.data.frame()
        # } else {
        #   sparsesig_aligned <- NULL
        #   sparsesig_long <- NULL 
        # }



        ### Summarize Signatures across combinations ###
        ground_truth_nested <- as.data.frame(ground_truth_nested)
        sigprof_aligned <- align_callers(tumourevo_signature_res = tumourevo_signature_res,tool = "SigProfiler",spn = spn_id) %>% as.data.frame()
        bascule_aligned <- align_callers(tumourevo_signature_res = tumourevo_signature_res,tool = "BASCULE",spn = spn_id) %>% as.data.frame()

        df_long <- bind_rows(
          ground_truth_nested %>% mutate(sample=row.names(.), Method="ProCESS"),
          bascule_aligned %>% mutate(sample=row.names(.), Method="BASCULE"),
          sigprof_aligned %>% mutate(sample=row.names(.), Method="SigProfiler")
        ) %>%
          pivot_longer(cols = starts_with(context_classes),
                       names_to = "Signature",
                       values_to = "Exposure")
        p_sankey <- ggplot(df_long,
                        aes(x = Method, y = Exposure,
                            stratum = Signature, alluvium = Signature,
                            fill = Signature, label = Signature)) +
          geom_flow(stat = "alluvium", lode.guidance = "forward", color = "black") +
          geom_stratum(color = "black") +
          scale_y_continuous(expand = c(0,0)) +
          facet_wrap(~sample, nrow = 1) +
          labs(title = paste0(spn_id,", cov=",coverage,", pur=", purity,", Signature class=",context_classes),
               y = "Exposure", x = "Method") +
          theme_minimal() +
          theme(legend.position = "bottom")
        ggsave(filename = paste0("sankey_plot_",context,".png"),plot =p_sankey ,width = 10,height = 5)

        ### Compare exposure of estimated and true signatures  ###
        
        gt_long <- reshape_exposures_long(exposures_mat = ground_truth_nested,spn = spn_id,
                                            coverage = coverage,purity = purity,method_name = "ProCESS")
        bascule_long <- reshape_exposures_long(exposures_mat=bascule_aligned,spn = spn_id,
                                                coverage = coverage,purity = purity,method_name = "BASCULE")
        sigprof_long <- reshape_exposures_long(exposures_mat=sigprof_aligned,spn = spn_id,
                                                coverage = coverage,purity = purity,method_name = "SigProfiler") %>% as.data.frame()
        
        combined <- bind_rows(gt_long, sigprof_long, bascule_long)
        # Filter to keep only exposures > 0
        all_exposures[[context]] <- combined %>% filter(Exposure > 0)
        
        
        
        # Signature and Exposure Validation
        metrics_sample <- tibble()
        metrics_spn <- tibble()
        cosine_mse <- tibble()

        sample_metrics_bascule <- per_sample_metrics(inferred_df = bascule_long, simulated_df = gt_long) %>% 
            mutate(caller = 'BASCULE', spn = spn_id, coverage = coverage, purity = purity)
        #sample_metrics_sparsesig <- per_sample_metrics(inferred_df = sparsesig_long, simulated_df = gt_long) %>%
            #mutate(caller = 'SparseSignatures', spn = spn_id, coverage = coverage, purity = purity)
        sample_metrics_sigprofiler <- per_sample_metrics(inferred_df = sigprof_long, simulated_df = gt_long) %>% 
            mutate(caller = 'SigProfiler', spn = spn_id, coverage = coverage, purity = purity)
        #summary_sparsesig <- summary_stats(sample_metrics_sparsesig) %>% mutate(caller = 'SparseSignatures', spn = spn_id, coverage = coverage, purity = purity)
        summary_sigprofiler <- summary_stats(sample_metrics_sigprofiler) %>% mutate(caller = 'SigProfiler', spn = spn_id, coverage = coverage, purity = purity)
        summary_bascule <- summary_stats(sample_metrics_bascule) %>% mutate(caller = 'BASCULE', spn = spn_id,coverage = coverage, purity = purity)
        # if(context=="SBS96"){
        #   cosine_mse_sparsesig <- compute_cosine_mse(inferred = sparsesig_aligned, simulated = ground_truth_nested) %>% 
        #     mutate(caller = 'SparseSignatures', spn = spn_id, coverage = coverage, purity = purity)
        # } else{
        #   cosine_mse_sparsesig<- NULL
        # } 
        cosine_mse_sigprofiler <- compute_cosine_mse(inferred = sigprof_aligned, simulated = ground_truth_nested) %>% 
            mutate(caller = 'SigProfiler', spn = spn_id, coverage = coverage, purity = purity)

        cosine_mse_bascule <- compute_cosine_mse(inferred = bascule_aligned, simulated = ground_truth_nested) %>% 
            mutate(caller = 'BASCULE', spn = spn_id, coverage = coverage, purity = purity)
        
        metrics_sample <- bind_rows(metrics_sample,  sample_metrics_sigprofiler,sample_metrics_bascule)
        metrics_spn <- bind_rows(metrics_spn, summary_sigprofiler,summary_bascule)
        cosine_mse <- bind_rows(cosine_mse, cosine_mse_sigprofiler, cosine_mse_bascule)
        saveRDS(object = metrics_spn,file = paste0("metrics_",context,"_spn.rds"))
        saveRDS(object = metrics_sample,file = paste0("metrics_",context,"_sample.rds"))
        saveRDS(object = cosine_mse,file = paste0("cosine_mse_",context,".rds"))
    }


    """

}
