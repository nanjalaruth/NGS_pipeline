#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//modules
include { QC; MULTIQC } from './modules/qc.nf'
include {MAPPING; FIXMATE; MARKDUP; INDEX } from './modules/mapping.nf'
include {VARIANT_CALLING; INDEX_VARIANTS} from './modules/variant_calling.nf'
include {COMPUTE_PERBASE_COVERAGE; COMPUTE_COVERAGE_HISTOGRAM} from './modules/coverage.nf'

/*
 * pipeline input parameters
 */
log.info """\
    N G S   P I P E L I N E - P L A S M O D I U M 
    =============================================
    reference: ${params.reference}
    reads        : ${params.reads}
    outdir       : ${params.outdir}
    region       : ${params.region}
    """
    .stripIndent()

/*
 * define the `index` process that creates a binary index
 * given the transcriptome file
 */


workflow {
    input = Channel.fromFilePairs( params.reads)
    
    //Quality Control
    fastqc_ch = QC(input)
    MULTIQC(fastqc_ch.collect())
    
    //Alignment
    ref = Channel.fromPath(params.reference)
    input_ch = input.combine(ref)
    MAPPING(input_ch)
    FIXMATE(MAPPING.out)
    final_bam = MARKDUP(FIXMATE.out) 
    index_bam = INDEX(MARKDUP.out)

    //Variant Calling
    region_ch = Channel.from(params.region)
    bam_ch = final_bam.combine(index_bam, by:0)
       .combine(ref)
       .combine(region_ch)
    //bam_ch.view()
    VARIANT_CALLING(bam_ch)
    INDEX_VARIANTS(VARIANT_CALLING.out)

    //Coverage
    COMPUTE_PERBASE_COVERAGE(final_bam)
    COMPUTE_COVERAGE_HISTOGRAM(final_bam)
}

workflow.onComplete {
    log.info ( workflow.success ? "\nDone!":"Oops .. something went wrong") 
}
