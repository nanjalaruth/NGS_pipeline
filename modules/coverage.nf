process COMPUTE_PERBASE_COVERAGE {

    tag "Computing per base coverage"
    publishDir "${params.outdir}/Coverage/${dataset}", mode: 'copy', overwrite: true

    input:
        tuple val(dataset), path(final_bam)
    output:
        tuple val(dataset), path(bedgraph)
    script:
        bedgraph = "${dataset}.perbase.bedgraph"
	"""
        bedtools genomecov -ibam ${final_bam} -bg > ${bedgraph}
	"""	
}

process COMPUTE_COVERAGE_HISTOGRAM {

    tag "Coverage on a histogram"
    publishDir "${params.outdir}/Coverage/${dataset}", mode: 'copy', overwrite: true

    input:
        tuple val(dataset), path(final_bam)
    output:
        tuple val(dataset), path(histogram)
    script:
        histogram = "{dataset}.percontig.txt"
	"""
        # add a header row
        echo -e 'contig\\tdepth_of_coverage\\tnumber_of_bases_at_this_depth\\tcontig_size\\tproportion_of_bases_at_this_depth' > ${histogram}
        bedtools genomecov -ibam ${final_bam} >> ${histogram}
	"""
}
