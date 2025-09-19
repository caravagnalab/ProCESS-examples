library(dplyr)
library(optparse)
source('../getters/process_getters.R')

option_list <- list(make_option(c("--spn_list"), type = "character", default = 'SPN03,SPN01'),
                    make_option(c("--coverages"), type = "character", default = '50,100,150'),
                    make_option(c("--purities"), type = "character", default = '0.3,0.6,0.9'))

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

spns = as.character(strsplit(opt$spn_list, ",")[[1]]) 
coverages <- as.numeric(strsplit(opt$coverages, ",")[[1]]) 
purities <- as.numeric(strsplit(opt$purities, ",")[[1]])

base = "/orfeo/cephfs/scratch/cdslab/shared/SCOUT/"
df = tibble()

for (spn in spns){
	gender = get_process_gender(spn)
        samples = c('normal_sample', get_sample_names(spn))
	for (s in samples){
		if (s == 'normal_sample'){
			tmp = tibble(spn = spn, sample = s, coverage = 50, purity = 0.9, directory = base, type = 0, sex = gender)
			df = bind_rows(df, tmp)
		} else {
			for (cov in coverages){
                        	for (pur in purities){
                                	tmp = tibble(spn = spn, sample = s, coverage = cov, purity = pur, directory = base, type = 1, sex = gender)
					df = bind_rows(df, tmp)
                        	}
                	}
		}
	}
}
write.table(x = df, file = 'samplesheet.csv', quote = F, sep = ',', row.names = F)
