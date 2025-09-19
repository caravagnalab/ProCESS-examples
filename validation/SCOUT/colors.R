#### HEATMAP annotation coloring

SPN_colors <-c("SPN01"='steelblue', "SPN02"='seagreen', "SPN03"='goldenrod', 
               "SPN04"='coral', "SPN06"='palevioletred', "SPN07"='indianred3')

purity_colors = c("0.3"="#bfd3e6", "0.6"="#8c96c6", "0.9"="#810f7c")
coverage_colors = c("50"="#ccece6", "100"="#66c2a4", "150"="#006d2c")


col_mut_types <-c("SNV"="black", "INDEL"="grey")
col_ccf_classes<-c("Subclonal Low CCF"="#f8cec9","Subclonal High CCF"="#f09e92","Clonal"="#e65d4a")

sbs_colors = setNames(
  nm = c("SBS1",
         "SBS17b",
         "SBS18",
         "SBS5",
         "SBS88",
         "SBS10b",
         "SBS6",
         "SBS9",
         "SBS25",
         "SBS4",
         "SBS11",
         "SBS3" ,
         "SBS26"), 
  object = c('#f1696bff', 
             '#8fbd8cff', 
             '#87c7d6ff', 
             '#bac3deff', 
             '#d7bfd9ff', 
             '#a8a2a1ff', 
             '#cfadb3ff', 
             '#3c609aff', 
             '#9a4564ff', 
             '#fbcb5bff', 
             '#c2b280ff', 
             '#d47e2dff', 
             '#5f8676ff')
)
id_colors = setNames(
  nm = c('ID1', 
         'ID2', 
         "ID4",  
         "ID5",  
         "ID7",  
         "ID8",   
         "ID9",  
         "ID18"), 
  object = c('#0c8281ff', 
             '#f5a55fff', 
             '#7d287eff', 
             '#2e4f4fff', 
             '#c4ddbcff', 
             '#996869ff', 
             '#daa627ff', 
             '#bc8f8fff')
)
