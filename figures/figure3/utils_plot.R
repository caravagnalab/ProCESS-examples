library(ggplot2)
library(ggpubr)
library(dplyr)
library(colormap)
library(RColorBrewer)
library(patchwork)

my_ggplot_theme <- function () 
{
  theme_light(base_size=10) +
  theme(legend.key.size=unit(0.3, "cm"),
        panel.background=element_rect(fill="white"),
        axis.text.x=element_text(size=8),
        axis.text.y=element_text(size=8),
        axis.title=element_text(size=10),
        legend.text=element_text(size=8),
        legend.title=element_text(size=10),
        text=element_text(size=10),
        plot.title=element_text(size=12))
}


#### HEATMAP annotation coloring

SPN_colors <-c("SPN01"='steelblue', "SPN02"='seagreen', "SPN03"='goldenrod', 
               "SPN04"='coral', "SPN05"="magenta4","SPN06"='palevioletred', "SPN07"='indianred3')

purity_colors = c("0.3"="#bfd3e6", "0.6"="#8c96c6", "0.9"="#810f7c")
coverage_colors = c("50"="#ccece6", "100"="#66c2a4", "150"="#006d2c")


col_mut_types <-c("SNV"="black", "INDEL"="grey")
col_ccf_classes<-c("Subclonal Low CCF"="#f8cec9","Subclonal High CCF"="#f09e92","Clonal"="#e65d4a")


col_cna_tools <- c("ASCAT"='deepskyblue4', "Battenberg"='maroon', "Sequenza"='sienna2')
col_somatic_tools = c(
  "ProCESS" = "gray80",
  "mutect2" = "lightsteelblue",
  "mutect2 (all)" = "steelblue4",  # darker steelblue
  "strelka" = "coral",
  "strelka (all)" = "coral4",     # darker coral
  "freebayes" = "#8FBC8B",
  "freebayes (all)" = "#228B22"   # darker green (ForestGreen )
)

col_germline_tools = c(
  "ProCESS" = "gray80",
  "mutect2" = "lightsteelblue",
  "mutect2 (all)" = "steelblue4",  # darker steelblue
  "strelka" = "coral",
  "strelka (all)" = "coral4",     # darker coral
  "freebayes" = "#8FBC8B",
  "freebayes (all)" = "#228B22",   # darker green (ForestGreen )
  'haplotypecaller' = 'palevioletred'
)

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
