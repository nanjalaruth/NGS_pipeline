process QC {
    tag "Quality control on ${sample_id}"
    //publishDir "${params.outdir}/fastq_preprocessing/${sample_id}", mode: 'copy', overwrite: true
    
    input:
        tuple val(sample_id), path(reads)
    output:
        path "fastqc_${sample_id}_logs"
    script:
        """
        mkdir fastqc_${sample_id}_logs
        fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
        """
}


process MULTIQC {
    publishDir "${params.outdir}/fastq_preprocessing", mode:'copy'
    tag "Compiling FASTQC reports"

    input:
        path "*"

    output:
        path "multiqc_report.html"

    script:
        """
        multiqc .
        """
}
