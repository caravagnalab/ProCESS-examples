library(IRanges)
library(GenomicRanges)


reciprocal_overlap <- function(pred, truth, threshold = 0.5) {
  hits <- findOverlaps(pred, truth)
  keep <- c()
  for (i in seq_along(hits)) {
    q <- queryHits(hits)[i]
    s <- subjectHits(hits)[i]
    inter_len <- width(pintersect(pred[q], truth[s]))
    ro1 <- inter_len / width(pred[q])
    ro2 <- inter_len / width(truth[s])
    if (ro1 >= threshold && ro2 >= threshold) {
      keep <- c(keep, i)
    }
  }
  hits[keep]
}


calculate_metrics <- function(gt, df_caller){
  gt_gr <- GRanges(seqnames = gt$chrom,
                      ranges = IRanges(start = gt$start,
                                       end = gt$end),
                      cn = gt$cn)
  
  caller_gr <- GRanges(seqnames = df_caller$chrom,
                     ranges = IRanges(start = df_caller$start,
                                      end = df_caller$end),
                     cn = df_caller$cn)
  
  # Find overlaps
  matched <- reciprocal_overlap(caller_gr, gt_gr, threshold = 0.5)
  TP <- length(matched)                     
  FP <- length(caller_gr) - length(unique(queryHits(matched)))
  FN <- length(gt_gr) - length(unique(subjectHits(matched))) 
  precision <- TP / (TP + FP)
  recall    <- TP / (TP + FN)
  f1_score  <- 2 * precision * recall / (precision + recall)
  
  return(list("precision"=precision,"recall"=recall,"f1_score"=f1_score))
  
}


color_by_state = c("INFERRED_Major1"='steelblue',
                   "INFERRED_minor1"='indianred',
                   "INFERRED_CN"='seagreen')

fill_by_match = c('match'= alpha('gainsboro', .03),
                  'no match subclone'= alpha('cadetblue3', .08),
                  'match subclone'= alpha('goldenrod', .08),
                  'no match' = alpha('indianred', .08),
                  'subclonal' = alpha('goldenrod', .08))

color_type <- c('clonal' = 'gainsboro', 'sub-clonal' = 'lightsalmon')
color_process <- c('Major1'='steelblue', 'minor1'='indianred', 'Major2'='steelblue1', 'minor2'="orangered")

absolute_to_relative_coordinates <- function(muts, reference = CNAqc::chr_coordinates_GRCh38, centromere = F){
  vfrom = reference$from
  names(vfrom) = reference$chr
  if (!centromere){
    muts %>%
      mutate(
        start = start + vfrom[chr],
        end = end + vfrom[chr])
  } else if (centromere){
    muts %>%
      mutate(
        start = start + vfrom[chr],
        end = end + vfrom[chr],
        centromere = centromere + vfrom[chr])
  }
}

absolute_to_relative_coordinates_muts <- function(muts, reference = CNAqc::chr_coordinates_GRCh38){
  vfrom = reference$from
  names(vfrom) = reference$chr
  muts <- muts %>% mutate(from = from + vfrom[chr])
  return(muts) 
}

## Read data
read_ProCESS = function(spn_id,sample_id,coverage,purity){
  
  # Read CNA data
  CNA_races = readRDS(get_process_cna(spn = spn_id,sample = sample_id)) %>%
    mutate(chr = paste0('chr',chr)) %>% dplyr::rename(from=begin,to=end) %>%
    mutate(seg_id = paste0(chr,':',from,':',to)) %>% as_tibble()
  
  return(list('CNA'=CNA_races))
}

read_ASCAT = function(spn_id,sample_id,coverage,purity){
  caller = 'ascat'
  ascat_results <- get_sarek_cna_file(spn = spn_id,
                                      coverage = coverage,
                                      purity = purity,
                                      caller = caller,
                                      type = "tumour",
                                      sampleID = sample_id)
  
  CNA_ascat = read.csv(ascat_results$cnvs, sep='\t') %>%
    mutate(chr = paste0('chr',chr)) %>% dplyr::rename(from=startpos,to=endpos,major=nMajor,minor=nMinor) %>% as_tibble()
  
  purity_ploidy = read.csv(ascat_results$purityploidy, sep='\t')
  return(list('CNA'=CNA_ascat, 'purity_ploidy'= purity_ploidy))
}

read_Sequenza = function(spn_id,sample_id,coverage,purity){
  caller = 'sequenza'
  sequenza_results <- get_sarek_cna_file(spn = spn_id,
                                         coverage = coverage,
                                         purity = purity,
                                         caller = caller,
                                         type = "tumour",
                                         sampleID = sample_id)
  
  CNA_sequenza = read.csv(sequenza_results$segments, sep='\t')
  CNA_sequenza = CNA_sequenza %>% select(!c(N.BAF,sd.BAF,N.ratio,sd.ratio,CNt,LPP))
  colnames(CNA_sequenza) = c('chr', 'from', 'to','BAF','DR','major','minor')
  purity_ploidy_sequenza = read.csv(sequenza_results$confints_CP, sep='\t')
  return(list('CNA'=CNA_sequenza, 'purity_ploidy'= purity_ploidy_sequenza))
}

read_CNVkit = function(spn_id,sample_id,coverage,purity){
  caller='cnvkit'
  cnvkit_results <- get_sarek_cna_file(spn = spn_id,
                                       coverage = coverage,
                                       purity = purity,
                                       caller = caller,
                                       type = "tumour",
                                       sampleID = sample_id)
  
  CNA_cnvkit = read.csv(cnvkit_results$somatic.call,sep='\t')
  CNA_cnvkit = CNA_cnvkit %>% select(!c(gene,ci_hi,ci_lo,depth,probes,weight))
  colnames(CNA_cnvkit) = c('chr', 'from','to','logR','CN')

  return(list('CNA'=CNA_cnvkit))
}

## Create joint segmentation
create_joint_segmentation = function(CNA_ProCESS, CNA_target, caller, chromosomes){
  
  joint_segmentation = lapply(chromosomes, function(c){
    
    CNA_ProCESS_chr = CNA_ProCESS %>% 
      filter(chr==c) %>% 
      mutate(group = cumsum(
        lag(major, default = dplyr::first(major)) != major |
          lag(minor, default = dplyr::first(minor)) != minor |
          lag(ratio, default = dplyr::first(ratio)) != ratio
      )) %>%
      group_by(chr, major, minor, ratio, group, from ,to) %>%
      select(chr, from, to, major, minor, ratio) %>% 
      arrange(from, to)
    
    CNA_target_chr = CNA_target %>% filter(chr==c)
    froms = c(CNA_ProCESS_chr$from %>% unique(), CNA_target_chr$from %>% unique())
    tos = c(CNA_ProCESS_chr$to %>% unique(), CNA_target_chr$to %>% unique())
    breakpoints = c(froms, tos) %>% sort() %>% unique()
    df = data.frame()
    
    for (i in 1:(length(breakpoints)-1)){
      f = breakpoints[i]
      t = breakpoints[i+1]
      true_cna= CNA_ProCESS_chr %>% filter(from <= f, from < t, to>=t) 
      inferred_cna= CNA_target_chr %>% filter(from <= f, from < t, to>=t) 
      if (nrow(true_cna)==0){
        true_M1=NA
        true_m1=NA
        true_M2=NA
        true_m2=NA
        true_cn=NA
      }
      if (nrow(true_cna)==1){
        true_M1=true_cna %>% pull(major)
        true_m1=true_cna %>% pull(minor)
        true_M2=NA
        true_m2=NA
        true_cn=true_M1+true_m1
      }
      if (nrow(true_cna)>1){
        true_cna = true_cna %>% arrange(desc(ratio))
        true_M1=true_cna[1,] %>% pull(major)
        true_m1=true_cna[1,] %>% pull(minor)
        true_M2=true_cna[2,] %>% pull(major)
        true_m2=true_cna[2,] %>% pull(minor)
        ratio=true_cna[1,] %>% pull(ratio)
        true_cn = ratio*(true_M1+true_m1) + (1-ratio)*(true_M2+true_m2) 
      }
      if (nrow(inferred_cna)==0){
        inferred_M=NA
        inferred_m=NA
        inferred_cn=NA
      }
      if (nrow(inferred_cna)==1){
        if (caller!='cnvkit'){
          inferred_M=inferred_cna %>% pull(major)
          inferred_m=inferred_cna %>% pull(minor)
          inferred_cn=inferred_M+inferred_m
        }else{
          inferred_M=NA
          inferred_m=NA
          inferred_cn=inferred_cna %>% pull(CN)
        }
      }
      df = rbind(df, data.frame(
        'from' =f,
        'to'=t,
        'chromosome'=c,
        'TRUE_Major1'= true_M1,
        'TRUE_minor1'= true_m1,
        'TRUE_Major2'= true_M2,
        'TRUE_minor2'= true_m2,
        'TRUE_CN' = true_cn,
        'INFERRED_Major1' = inferred_M,
        'INFERRED_minor1' = inferred_m,
        'INFERRED_CN' = inferred_cn,
        'caller'=caller
      ))
    }
    df
    
  }) #%>% Reduce(rbind)
  joint_segmentation = Reduce(rbind, joint_segmentation)
  
  joint_segmentation_shifted = lapply(1:nrow(joint_segmentation), function(r){
    chromosome = joint_segmentation[r,]$chromosome
    from_chromosome = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(from)
    joint_segmentation[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome)
  })
  joint_segmentation_shifted = Reduce(rbind, joint_segmentation_shifted)
  
  if (caller %in% c('ascat','sequenza')){
    joint_segmentation_shifted_longer = joint_segmentation_shifted %>% 
      mutate(is_match = 
               case_when(
                 TRUE_Major1==INFERRED_Major1 & TRUE_minor1==INFERRED_minor1 & is.na(TRUE_Major2) ~ 'match',
                 TRUE_Major1==INFERRED_Major1 & TRUE_minor1==INFERRED_minor1 | 
                 TRUE_Major1==INFERRED_Major1 & TRUE_minor1==INFERRED_minor1 & !is.na(TRUE_Major2) ~ 'match subclone',
                 TRUE_Major1!=INFERRED_Major1 & TRUE_minor1!=INFERRED_minor1 & !is.na(TRUE_Major2) ~ 'no match subclone',
                 .default = 'no match'
               )) %>%
      mutate(type = ifelse(!is.na(TRUE_Major2) & !is.na(TRUE_Major2), 'subclonal', 'clonal')) %>% 
      pivot_longer(
        cols = c('TRUE_Major1', 'TRUE_minor1', 'TRUE_Major2', 'TRUE_minor2', 'INFERRED_Major1', 'INFERRED_minor1'),  
        names_to = "Type",  # New column for variable names
        values_to = "Value"  # New column for values
      )
    joint_segmentation_shifted = joint_segmentation_shifted %>% filter(!is.na(TRUE_CN))
    joint_segmentation_shifted_longer = joint_segmentation_shifted_longer %>% filter(!is.na(TRUE_CN))
  } else{
    joint_segmentation_shifted = joint_segmentation_shifted %>% filter(!is.na(TRUE_CN))
    joint_segmentation_shifted_longer = joint_segmentation_shifted %>% 
      mutate(is_match = 
               case_when(
                 abs(TRUE_CN - INFERRED_CN)<.1 ~ 'match',
                 .default = 'no match'
               )) %>% 
      mutate(type = ifelse(!is.na(TRUE_Major2) & !is.na(TRUE_Major2), 'subclonal', 'clonal')) %>% 
      select(from, to, chromosome, TRUE_CN, INFERRED_CN,is_match, type) %>%
      mutate(is_match = ifelse(type == 'subclonal', 'subclonal', is_match)) %>% 
      pivot_longer(
        cols = c('TRUE_CN', 'INFERRED_CN'),  
        names_to = "Type",  # New column for variable names
        values_to = "Value"  # New column for values
      )
  }
  
  return(list('joint_segmentation'=joint_segmentation_shifted,
              'joint_segmentation_long'=joint_segmentation_shifted_longer))
}

## Shift segments
shift_segments = function(CNA){
  CNA_shifted = lapply(1:nrow(CNA), function(r){
  chromosome = CNA[r,]$chr
  from_chromosome = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(from)
  
  if ("original_from" %in% colnames(CNA)){
    CNA[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome,
                       original_from = original_from + from_chromosome, original_to = original_to + from_chromosome)
  }else{
    CNA[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome)
  }
  
  })
  CNA_shifted = Reduce(rbind, CNA_shifted)
  CNA_shifted
}

## Compute correctness 
compute_correctness = function(df, caller) {
  if (caller %in% c('ascat', 'sequenza')){
    clonal = 1-( (df %>% filter(is_match == 'no match') %>% filter(type == 'clonal') %>% mutate(len=to-from) %>%
           pull(len) %>% unique() %>% sum()) / (df %>% filter(type == 'clonal') %>% mutate(len=to-from) %>%
                                                  pull(len) %>% unique() %>% sum()))
    
    all = 1-( (df %>% filter(is_match %in% c('no match', 'no match subclone')) %>% mutate(len=to-from) %>%
                    pull(len) %>% unique() %>% sum()) / (df %>% mutate(len=to-from) %>%
                                                           pull(len) %>% unique() %>% sum()))
  } else{
    clonal = 1-( (df %>% filter(is_match == 'no match') %>% filter(type == 'clonal') %>% mutate(len=to-from) %>%
                    pull(len) %>% unique() %>% sum()) / (df %>% filter(type == 'clonal') %>% mutate(len=to-from) %>%
                                                           pull(len) %>% unique() %>% sum()))
    all <- NA
  }
  
  return(list('all'=all, 'clonal'=clonal))
}

## Expected BAF and DR functions
expected_BAF = function(nA1, nB1, nA2, nB2, purity, ccf){
  if (is.na(nA2) & is.na(nB2)){
    nA2 = 0
    nB2 = 0
  }
  nom = min(nA1*ccf + nA2*(1-ccf), nB1*ccf + nB2*(1-ccf))*purity + (1-purity)
  den = ((nA1+nB1)*ccf + (nA2+nB2)*(1-ccf))*purity + 2*(1-purity)
  nom/den
}
expected_DR = function(nA1, nB1, nA2, nB2, purity, ccf, ploidy){
  if (is.na(nA2) & is.na(nB2)){
    nA2 = 0
    nB2 = 0
  }
  nom = ((nA1+nB1)*ccf + (nA2+nB2)*(1-ccf))*purity + 2*(1-purity)
  nom/(ploidy*purity + 2*(1-purity))
}

## True Ploidy
compute_true_ploidy = function(CNA_ProCESS){
  # Compute real ploidy
  CNA_ProCESS = CNA_ProCESS %>% mutate(len=to-from)
  genome_len =CNA_ProCESS$len %>% unique() %>% sum()
  CNA_ProCESS = CNA_ProCESS %>% mutate(seg_id = paste0(chr, ':', from,':',to))
  segments = CNA_ProCESS$seg_id %>% unique()
  ploidy = 0
  # purity_number = as.double(purity)
  for (s in segments){
    seg = CNA_ProCESS %>% filter(seg_id == s)
    nA1 = seg[1,]$major
    nB1 = seg[1,]$minor
    genome_fraction = seg[1,]$len/genome_len
    if (nrow(seg)>1){
      nA2 = seg[2,]$major
      nB2 = seg[2,]$minor
      ccf = seg[1,]$ratio
    }else{
      nA2 = 0
      nB2 = 0
      ccf = 1
    }
    
    # ploidy = ploidy + ((nA1+nB1)*ccf + (nA2+nB2)*(1-ccf))*purity_number*genome_fraction
    ploidy = ploidy + ((nA1+nB1)*ccf + (nA2+nB2)*(1-ccf))*genome_fraction
  }
  # ploidy = ploidy + 2*(1-purity_number)
  ploidy
}

## Theoretical vs inferred BAF and DR (ascat)
theoretical_vs_inferred_BAF_DR = function(CNA_ProCESS,CNA_ascat,purity_number,ploidy,purity_ploidy){
  already_seen_segment = c()
  CNA_ProCESS = lapply(1:nrow(CNA_ProCESS), function(r){
    chromosome = CNA_ProCESS[r,]$chr
    from_chromosome = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(from)
    CNA_ProCESS[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome)
  })
  CNA_ProCESS = Reduce(rbind, CNA_ProCESS)
  transposed_CNA_ProCESS= lapply(CNA_ProCESS$seg_id%>%unique(), function(s){
    CNA_ProCESS_s = CNA_ProCESS %>% filter(seg_id==s)
    if (nrow(CNA_ProCESS_s)>1){
      data.frame("chr"=CNA_ProCESS_s[1,]$chr,
                 "from"=CNA_ProCESS_s[1,]$from,
                 "to"=CNA_ProCESS_s[1,]$to,
                 "Major1"=CNA_ProCESS_s[1,]$major,
                 "minor1"=CNA_ProCESS_s[1,]$minor,
                 "Major2"=CNA_ProCESS_s[2,]$major,
                 "minor2"=CNA_ProCESS_s[2,]$minor,
                 "ratio"=CNA_ProCESS_s[1,]$ratio)
    }else{
      data.frame("chr"=CNA_ProCESS_s[1,]$chr,
                 "from"=CNA_ProCESS_s[1,]$from,
                 "to"=CNA_ProCESS_s[1,]$to,
                 "Major1"=CNA_ProCESS_s[1,]$major,
                 "minor1"=CNA_ProCESS_s[1,]$minor,
                 "Major2"=NA,
                 "minor2"=NA,
                 "ratio"=CNA_ProCESS_s[1,]$ratio)
    }
  })
  transposed_CNA_ProCESS = Reduce(rbind,transposed_CNA_ProCESS)
  CNA_ascat = lapply(1:nrow(CNA_ascat), function(r){
    chromosome = CNA_ascat[r,]$chr
    from_chromosome = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(from)
    CNA_ascat[r,] %>% mutate(from = from + from_chromosome, to = to + from_chromosome)
  })
  CNA_ascat = Reduce(rbind, CNA_ascat)
  
  simulated_baf_dr = lapply(1:nrow(transposed_CNA_ProCESS), function(r){
    simulated_baf=expected_BAF(transposed_CNA_ProCESS[r,]$Major1, transposed_CNA_ProCESS[r,]$minor1, transposed_CNA_ProCESS[r,]$Major2, transposed_CNA_ProCESS[r,]$minor2, purity_number, transposed_CNA_ProCESS[r,]$ratio)
    simulated_dr=expected_DR(transposed_CNA_ProCESS[r,]$Major1, transposed_CNA_ProCESS[r,]$minor1,transposed_CNA_ProCESS[r,]$Major2, transposed_CNA_ProCESS[r,]$minor2, purity_number, transposed_CNA_ProCESS[r,]$ratio, ploidy)
    cbind(transposed_CNA_ProCESS[r,], data.frame('simulated_baf'=c(simulated_baf), 'simulated_dr'=c(simulated_dr)))
  })
  simulated_baf_dr = Reduce(rbind, simulated_baf_dr)
  
  inferred_baf_dr = lapply(1:nrow(CNA_ascat), function(r){
    inferred_baf=expected_BAF(CNA_ascat[r,]$major, CNA_ascat[r,]$minor, 0, 0, purity_ploidy$AberrantCellFraction, 1)
    inferred_dr=expected_DR(CNA_ascat[r,]$major, CNA_ascat[r,]$minor, 0, 0, purity_ploidy$AberrantCellFraction, 1, purity_ploidy$Ploidy)
    cbind(CNA_ascat[r,], data.frame('inferred_baf'=c(inferred_baf), 'inferred_dr'=c(inferred_dr)))
  })
  inferred_baf_dr = Reduce(rbind,inferred_baf_dr)
  
  return(list('simulated_baf_dr' = simulated_baf_dr,'inferred_baf_dr_ascat' = inferred_baf_dr))
}

## Compute segmentation correctness
covered_genome = function(CNA_target, chromosome){
  chr_len = CNAqc:::get_reference('hg38') %>% filter(chr==chromosome) %>% pull(length)
  covered = CNA_target %>% filter(chr == chromosome) %>% mutate(l=to-from) %>% pull(l) %>% sum()
  (covered / chr_len)*100
}

breakpoint_analysis = function(CNA_ProCESS, CNA_target, chromosome, th=1e7){
  ProCESS_BP = CNA_ProCESS %>% filter(chr==chromosome) #%>% select(from, to)
  
  BP_f = ProCESS_BP$from %>% unique()
  BP_t = ProCESS_BP$to %>% unique()
  
  CNA_target_c = CNA_target %>% filter(chr==chromosome)
  BP_target_f = CNA_target_c$from %>% unique()
  BP_target_t = CNA_target_c$to %>% unique()
  
  if (nrow(CNA_target_c)==0){
    BP_df = data.frame(
      'og_coord'=c(BP_f,BP_f),
      'coord'=rep(NA, length(c(BP_f,BP_f))),
      'dist'=rep(NA,length(c(BP_f,BP_f))),
      'state'=rep('no bp',length(c(BP_f,BP_f)))
    )
  }else{
  
  BP_df = data.frame()
  
  for (bp in BP_f){
    dist = min(abs(bp - BP_target_f))
    coord = BP_target_f[which.min(abs(bp - BP_target_f))]
    if (dist < th){
      df = data.frame(
        'og_coord'=bp,
        'coord'=coord,
        'dist'=dist,
        'state'='matching'
      )
    }else{
      df = data.frame(
        'og_coord'=bp,
        'coord'=NA,
        'dist'=NA,
        'state'='missed'
      )
    }
    BP_df = rbind(BP_df, df)
  }
  
  for (bp in BP_t){
    dist = min(abs(bp - BP_target_t))
    coord = BP_target_t[which.min(abs(bp - BP_target_t))]
    if (dist < th){
      df = data.frame(
        'og_coord'=bp,
        'coord'=coord,
        'dist'=dist,
        'state'='matching'
      )
    }else{
      df = data.frame(
        'og_coord'=bp,
        'coord'=NA,
        'dist'=NA,
        'state'='missed'
      )
    }
    BP_df = rbind(BP_df, df)
  }
  
  for (bp in BP_target_f){
    if (!(bp %in% BP_df$coord)){
      df = data.frame(
        'og_coord'=NA,
        'coord'=bp,
        'dist'=NA,
        'state'='added'
      )
      BP_df = rbind(BP_df, df)
    }
  }
  
  for (bp in BP_target_t){
    if (!(bp %in% BP_df$coord)){
      df = data.frame(
        'og_coord'=NA,
        'coord'=bp,
        'dist'=NA,
        'state'='added'
      )
      BP_df = rbind(BP_df, df)
    }
  }
  }
  return(BP_df)
}
segmentation_analysis = function(CNA_ProCESS, CNA_target, chromosomes, th=1e7){
  # Mean percentage of genome covered by segments
  mean_covered_genome = lapply(chromosomes, function(c){
    covered_genome(CNA_target, c)
  }) %>% unlist() %>% mean()
  
  # Breakpoints
  BP = lapply(chromosomes, function(chr){
    #print(chr)
    bp = breakpoint_analysis(CNA_ProCESS, CNA_target, chr, th=th)
    bp$chromosome = chr
    bp
  })
  BP = Reduce(rbind, BP)
  
  # % of missed breakpoints (on total)
  missed_bp = (BP %>% filter(state == 'missed') %>% nrow()) / (BP %>% filter(state %in% c('missed','matching')) %>% nrow()) * 100
  # % of added breakpoints (on total)
  added_bp = (BP %>% filter(state == 'added') %>% nrow()) #/ (BP %>% filter(state %in% c('missed','matching')) %>% nrow()) * 100
  # average distance
  av_dist = mean(BP %>% filter(state == 'matching') %>% pull(dist))
  
  summary_df = data.frame(
    'missed_bp'=missed_bp,
    'added_bp'=added_bp,
    'av_distance'=av_dist
  )
  
  return(list('summary'=summary_df, 'breakpoints' = BP))
}
my_blank_genome = function (ref = "GRCh38", genomic_coords, chromosomes = paste0("chr", 
                                                                                 c(1:22, "X", "Y")), label_chr = -0.5, cex = 1) 
{
  reference_coordinates = CNAqc:::get_reference(ref, genomic_coords) %>% 
    filter(chr %in% chromosomes)
  low = min(reference_coordinates$from)
  upp = max(reference_coordinates$to)
  p1 = ggplot2::ggplot(reference_coordinates) + CNAqc:::my_ggplot_theme(cex = cex) + 
    ggplot2::geom_segment(ggplot2::aes(x = centromerStart, 
                                       xend = centromerEnd, y = 0, yend = Inf), size = 0, 
                          color = "black", linetype = 8)
  p1 = p1 + ggplot2::geom_rect(data = reference_coordinates, 
                               ggplot2::aes(xmin = from, xmax = from, ymin = 0, ymax = Inf), 
                               alpha = 0, colour = "grey", size=0)
  p1 = p1 + ggplot2::geom_hline(yintercept = 0, size = 1, colour = "gainsboro") + 
    ggplot2::geom_hline(yintercept = 1, size = 0.3, colour = "black", 
                        linetype = "dashed") + ggplot2::labs(x = "Chromosome", 
                                                             y = "Major/ minor allele") + ggpubr::rotate_y_text() + 
    ggplot2::scale_x_continuous(breaks = c(0, reference_coordinates$centromerStart, 
                                           upp), labels = c("", gsub(pattern = "chr", replacement = "", 
                                                                     reference_coordinates$chr), ""))+
    theme(
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
    ) 
  return(p1)
}





