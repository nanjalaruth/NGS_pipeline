process VARIANT_CALLING {
    tag "Calling variants from the bam file using Octopus"
    publishDir "${params.outdir}/Variant_Calling/${dataset}", mode: 'copy', overwrite: true

    input:
        tuple val(dataset), path(final_bam), path(indexed_bam), path(ref),  val(region1), val(region2), val(region3)
    output:
        tuple val(dataset), path(vcf_file)
    script:
        vcf_file = "${dataset}.vcf"
        """
        gunzip -c ${ref} > ${ref}.fasta
        samtools faidx ${ref}.fasta
        octopus -R ${ref}.fasta -I ${final_bam} -P 1 -o ${vcf_file} --regions ${region1} ${region2} ${region3}
        """
}

process INDEX_VARIANTS {
    tag "Indexing variants"
    publishDir "${params.outdir}/Variant_Calling/${dataset}", mode: 'copy', overwrite: true

    input:
        tuple val(dataset), path(vcf_file)
    output:
        tuple val(dataset), path("${dataset}.vcf.gz.tbi")
    script:
        """
        bgzip ${vcf_file} 
        tabix ${vcf_file}.gz
        """
}
