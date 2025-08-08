#!/usr/bin/env nextflow
// from https://github.com/chelauk/nf-core-mutectplatypus/blob/master/workflows/mutect_platypus.nf
nextflow.enable.dsl=2

include { SAMTOOLS_CONVERT  as  SAMTOOLS_CONVERT_T    } from './modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT  as  SAMTOOLS_CONVERT_N    } from './modules/nf-core/samtools/convert/main'
include { SEQUENZAUTILS_GCWIGGLE  } from './modules/nf-core/sequenzautils/gcwiggle/main'
include { SEQUENZAUTILS_BAM2SEQZ  } from './modules/nf-core/sequenzautils/bam2seqz/main'
include { SEQUENZAUTILS_MERGESEQZ } from './modules/local/sequenzautils/mergeseqz/main'
include { SEQUENZAUTILS_BINNING   } from './modules/local/sequenzautils/seqzbin/main'
include { SEQUENZAUTILS_RSEQZ     } from './modules/local/sequenzautils/seq_R/main'
include { BATTENBERG  } from './modules/local/battenberg/main'


include { samplesheetToList } from 'plugin/nf-schema'


workflow SEQUENZA {
    take:
    tumour_bam
    normal_bam
    fasta
    fai

    main:
           
    tumour_bam
            .combine(normal_bam, by:0)
            .map { patient, meta1, files1, files11, meta2, files2, files22 ->
                    [meta1 + [patient:patient], files2, files1, files22, files11]
            }
            .set{seq_input_matched}
    
    if (params.create_wiggle) {
        SEQUENZAUTILS_GCWIGGLE(fasta)
        wiggle = SEQUENZAUTILS_GCWIGGLE.out.wig // something wrong here
    } else {
        wiggle = Channel.fromPath(params.wiggle).collect() 
    }

    fai.map{meta, chr -> [chr]}.set{fai_fasta}
    seqz_chr = fai_fasta.splitCsv(sep: "\t")
        .map{ chr -> chr[0][0] }
        .filter( ~/^chr\d+|^chr[X,Y]|^\d+|[X,Y]/ )
        .collect()
    
    SEQUENZAUTILS_BAM2SEQZ(seq_input_matched,
                            fasta,
                            wiggle,
                            seqz_chr)

    SEQUENZAUTILS_BAM2SEQZ.out.seqz
                            .groupTuple()
                            .set{merge_seqz_input}

    SEQUENZAUTILS_MERGESEQZ(merge_seqz_input)
    SEQUENZAUTILS_BINNING(SEQUENZAUTILS_MERGESEQZ.out.concat_seqz, params.bin)
    SEQUENZAUTILS_RSEQZ(SEQUENZAUTILS_BINNING.out.seqz_bin)
}

workflow  {   
    fasta           = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()
    fai             = params.fasta ? Channel.fromPath(params.fai).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

    input_samplesheet = params.input ? Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json")) : Channel.empty()
    
    input_samplesheet.branch { meta, cram, crai ->
                        tumour: meta.status == 1
                        normal: meta.status == 0
                    }.set{seq_split}

    // process tumour
    seq_split.tumour.map{ meta, cram, crai -> 
                        meta = meta + [id:meta.sample+'_vs_normal_sample']
                        [meta, cram, crai]
                    }.set{tumour_cram}
                    
    SAMTOOLS_CONVERT_T(tumour_cram,fasta,fai)
    out_bam_T = SAMTOOLS_CONVERT_T.out.bam.join(SAMTOOLS_CONVERT_T.out.bai)
    out_bam_T.map{ meta, bam, bai -> 
            [meta.patient, [id:meta.id , status:meta.status, gender:meta.gender], bam, bai]
        }.set{tumour_bam}

    // process normal        
    seq_split.normal.map{ meta, cram, crai -> 
                        meta = meta + [id:meta.patient + '_' + meta.sample]
                        [meta, cram, crai]
                     }.set{normal_cram} 

    SAMTOOLS_CONVERT_N(normal_cram,fasta,fai)
    out_bam_N = SAMTOOLS_CONVERT_N.out.bam.join(SAMTOOLS_CONVERT_N.out.bai)
    out_bam_N.map{ meta, bam, bai -> 
                [meta.patient, [ id:meta.id, status:meta.status, gender:meta.gender], bam, bai]
            }.set{normal_bam}

    // run sequenza
    if (params.tool && params.tool.split(',').contains('sequenza')){
        SEQUENZA(tumour_bam, normal_bam, fasta, fai)
    }
    
    // run battenberg
    if (params.tool && params.tool.split(',').contains('battenberg')){
    tumour_bam
            .combine(normal_bam, by:0)
            .map { patient, meta1, files1, files11, meta2, files2, files22 ->
                    [meta1 + [patient:patient], files2, files1, files22, files11]
            }
            .set{input_join}
    BATTENBERG(input_join)
    }

}
