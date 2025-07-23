process SOMATIC_PROCESSOR_MUTECT2 {

    tag { "${row}_${caller}_chr${chromosome}" }

    input:
    tuple val(row), val(chromosome), val(caller)

    output:
    //path "${row.spn}_${row.coverage}_${row.purity}_chr${chromosome}_${caller}_output.txt"
    path "chr${chromosome}.rds"
    publishDir "${params.outdir}/somatic/${row.spn}/${row.coverage}x_${row.purity}p/mutect2", mode: 'copy'

    script:
    """

    #!/usr/bin/env Rscript
    library(ProCESS)
    source("${projectDir}/bin/getters/process_getters.R")
    source("${projectDir}/bin/getters/sarek_getters.R")
    source("${projectDir}/bin/somatic/utils/mutect_utils.R")
    source("${projectDir}/bin/somatic/utils/process_utils.R")

    write("$row.spn", file = "${row.spn}_${row.coverage}_${row.purity}_chr${chromosome}_${caller}_output.txt")

    data_dir = "$row.directory" #'/orfeo/scratch/cdslab/shared/SCOUT'
    spn_id = "$row.spn"
    coverage = "$row.coverage"
    purity = "$row.purity"
    chromosome <- as.numeric("$chromosome")
    
    outdir <-  file.path("${params.outdir}","/somatic/")
    
    gt_path = get_mutations(spn = spn_id, 
                            base_path = data_dir, 
                            coverage = coverage, 
                            purity = purity, 
                            type = "tumour")
    print(gt_path)
    process_mutect2_results(gt_path, spn_id, purity, coverage, chromosome, base_path = data_dir, outdir = outdir)
    """
}
