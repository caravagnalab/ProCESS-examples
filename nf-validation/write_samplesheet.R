library(dplyr)
source('../getters/process_getters.R')
spns = c('SPN01', 'SPN02', 'SPN03', 'SPN04', 'SPN06', 'SPN07')
coverages = c(50)
purities = c(0.3, 0.6, 0.9)

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
