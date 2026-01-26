upper <- 0.89
lower <- 0.11

get_fga_fgs = function(x) {
  final_CNA_ProCESS = x %>% 
    dplyr::rename(from = begin) %>% 
    dplyr::rename(to = end) %>% 
    mutate(CN_type = ifelse(ratio < upper & ratio > lower, 'sub-clonal', 'clonal')) %>% 
    filter(!(CN_type == 'sub-clonal' & ratio  > upper)) %>% 
    filter(!(CN_type == 'clonal' & ratio <= lower))  %>% 
    mutate(ratio = ifelse(ratio < upper & ratio > lower, ratio, 1)) %>% 
    mutate(ratio = round(ratio, digits=1))
  
  tot_genome = final_CNA_ProCESS  %>% 
    filter(!(chr %in% c('chrX', 'chrY')))  %>% 
    select(chr, from, to) %>% 
    distinct() %>%  
    mutate(len=to-from) %>% 
    pull(len) %>% 
    unique() %>% 
    sum()
  
  altered = final_CNA_ProCESS %>% 
    filter(!(chr %in% c('chrX', 'chrY'))) %>% 
    mutate(len = to-from, CN = paste(major, minor, sep=':')) %>% 
    filter(ratio < 1 | CN !='1:1') %>% 
    select(-ratio, -CN, -major, -minor) %>% 
    distinct() %>% 
    pull(len) %>%
    unique() %>% 
    sum()
  fga = (altered/tot_genome)*100
  
  subclonal = final_CNA_ProCESS %>% 
    filter(!(chr %in% c('chrX', 'chrY'))) %>% 
    mutate(len = to-from, CN = paste(major, minor, sep=':')) %>% 
    filter(ratio < 1) %>% 
    select(-ratio, -CN, -major, -minor) %>% 
    distinct() %>% 
    pull(len) %>% 
    unique() %>% 
    sum()
  
  fgs = (subclonal/tot_genome)*100
  
  return(list(
    'fga' = fga, 
    'fgs' = fgs
  ))
}




