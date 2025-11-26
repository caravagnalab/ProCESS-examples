library(ProCESS)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(colormap)
library(RColorBrewer)
library(patchwork)
library(ComplexHeatmap)
library(GenomicRanges)



map_driver_to_cna_process <- function(simulated_drivers_with_vaf,cna_process,coverage,purity){
  mut_gr <- GRanges(
    seqnames = simulated_drivers_with_vaf$chr,
    ranges = IRanges(
      start = simulated_drivers_with_vaf$from,
      end = simulated_drivers_with_vaf$end
    ),
    mutationID = simulated_drivers_with_vaf$mutationID,
    sample_name = simulated_drivers_with_vaf$sample_name
  )
  cn_gr <- GRanges(
    seqnames = cna_process$chr,
    ranges = IRanges(
      start = cna_process$begin,
      end = cna_process$end
    ),
    major = cna_process$major,
    minor = cna_process$minor,
    ratio = cna_process$ratio
  )
  
  hits <- findOverlaps(mut_gr, cn_gr)
  
  # Extract indices
  mut_idx <- queryHits(hits)
  cn_idx <- subjectHits(hits)
  
  # Copy CN info into the mutations dataframe
  
  mut_with_cn <- simulated_drivers_with_vaf[mut_idx, ]
  mut_with_cn$major <- cn_gr$major[cn_idx]
  mut_with_cn$minor <- cn_gr$minor[cn_idx]
  mut_with_cn$ratio <- cn_gr$ratio[cn_idx]
  
  mut_with_cn <- mut_with_cn %>% 
    mutate(coverage=coverage,
           purity=purity) %>% 
    mutate(karyotype=paste(major,minor,sep=":"))
    # select(mutationID,code,SPN,sample_name,VAF,karyotype,ratio,coverage,purity)
  if (nrow(mut_with_cn)>1){
    mut_with_cn <- mut_with_cn %>% 
      mutate(type="subclonal")
  } else {
    mut_with_cn <- mut_with_cn %>% 
      mutate(type="clonal")
  }
  return(mut_with_cn)
}




map_driver_to_cna_ascat <- function(simulated_drivers_with_vaf_mutect2,cna_caller,coverage,purity){
  mut_gr <- GRanges(
    seqnames = simulated_drivers_with_vaf_mutect2$chr,
    ranges = IRanges(
      start = simulated_drivers_with_vaf_mutect2$from,
      end = simulated_drivers_with_vaf_mutect2$end
    ),
    mutationID = simulated_drivers_with_vaf_mutect2$mutationID
  )
  cn_gr <- GRanges(
    seqnames = cna_caller$chr,
    ranges = IRanges(
      start = cna_caller$from,
      end = cna_caller$to
    ),
    major = cna_caller$major,
    minor = cna_caller$minor
  )
  
  hits <- findOverlaps(mut_gr, cn_gr)
  
  # Extract indices
  mut_idx <- queryHits(hits)
  cn_idx <- subjectHits(hits)
  
  # Copy CN info into the mutations dataframe
  
  mut_with_cn <- simulated_drivers_with_vaf_mutect2[mut_idx, ]
  mut_with_cn$major <- cn_gr$major[cn_idx]
  mut_with_cn$minor <- cn_gr$minor[cn_idx]
  mut_with_cn$ratio <- cn_gr$ratio[cn_idx]
  
  mut_with_cn <- mut_with_cn %>% 
    mutate(coverage=coverage,
           purity=purity) %>% 
    mutate(karyotype=paste(major,minor,sep=":"))
    # select(mutationID,code,SPN,sample_name,VAF,karyotype,ratio,coverage,purity)
  if (nrow(mut_with_cn)>1){
    mut_with_cn <- mut_with_cn %>% 
      mutate(type="subclonal")
  } else {
    mut_with_cn <- mut_with_cn %>% 
      mutate(type="clonal")
  }
  return(mut_with_cn)
}


expected_vaf<-function(m,ccf,purity,karyotype){
  n=strsplit(x = karyotype,split = ":") %>% unlist()
  na=as.numeric(n[1])
  nb=as.numeric(n[2])
  vaf = (m*purity*ccf)/((2*(1-purity))+(purity*(na+nb)))
  return(vaf)
}


expected_ccf<-function(m,vaf,purity,karyotype){
  n=strsplit(x = karyotype,split = ":") %>% unlist()
  na=as.numeric(n[1])
  nb=as.numeric(n[2])
  ccf=((vaf*((na+nb-2)*purity+2))/(m*purity))
  return(ccf)
}


process_ccf <- function(sample_forest,phylo_forest){
  driver_sampled_cell <- phylo_forest$get_sampled_cell_mutations() %>% filter(class=="driver")
  samples_info <-  sample_forest$get_samples_info() %>% 
    dplyr::rename(sample=name)
  ccf_table <- sample_forest$get_nodes() %>% 
    inner_join(samples_info) %>% 
    inner_join(driver_sampled_cell) %>% 
    select(cell_id,sample, cause,tumour_cells_in_bbox) %>% 
    unique() %>% 
  group_by(sample, cause,tumour_cells_in_bbox) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(
      ccf = n / tumour_cells_in_bbox
    ) %>%
    ungroup() %>% 
    dplyr::rename(mutant=cause)
  return(ccf_table)
}

