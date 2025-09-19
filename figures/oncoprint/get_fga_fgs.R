get_fga_fgs = function(x) {
  final_CNA_ProCESS = x %>% 
    dplyr::rename(from = begin) %>% 
    rename(to = end)
  
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




