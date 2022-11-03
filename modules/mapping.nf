process MAPPING {
    tag "Mapping ${dataset} to the reference genome"
    //publishDir "${params.outdir}/Mapping/${dataset}", mode: 'copy', overwrite: true
    
    input:
        tuple val(dataset), path(reads), path(ref)
    output:
        tuple val(dataset), path(mapping_out)
    script:
        mapping_out = "${dataset}.sam"
        if( !params.single_end)
            """
            bwa index ${ref}
            bwa mem ${ref} ${reads[0]} ${reads[1]} -R '@RG\\tID:None\\tSM:None\\tLB:None\\tPL:Illumina,g' > ${mapping_out}
            """
        else
            """
            bwa index ${ref}
            bwa mem ${ref} ${reads[0]} > ${mapping_out}  
            """
}

process FIXMATE {
    tag "Correcting mate pair position, converting to bam and sorting the bam file"
    //publishDir "${params.outdir}/Mapping/${dataset}", mode: 'copy', overwrite: false

    input:
        tuple val(dataset), path(sam_file)
    output:
        tuple val(dataset), path(sorted_bam)
    script:
        bam_file = "${dataset}.bam"
        sorted_bam = "${dataset}_sorted.bam"
        """
        samtools fixmate -rpcm ${sam_file} ${bam_file} -O BAM
        samtools sort -@ 20 -O bam -o ${sorted_bam} ${bam_file}
        """
}

process MARKDUP {
    tag "Marking duplicates on the BAM file"
    publishDir "${params.outdir}/Mapping/${dataset}", mode: 'copy', overwrite: false

    input:
        tuple val(dataset), path(sorted_bam)
    output:
        tuple val(dataset), path(final_bam)
    script:
        final_bam = "${dataset}_final.bam"
        """
        #remove all the duplicates and also print some basic stats about the result file
        samtools markdup -r -s ${sorted_bam} ${final_bam}
        """
}

process INDEX {
    tag "Marking duplicates on the BAM file"
    publishDir "${params.outdir}/Mapping/${dataset}", mode: 'copy', overwrite: false

    input:
        tuple val(dataset), path(final_bam)
    output:
        tuple val(dataset), path(indexed_final_bam)
    script:
        indexed_final_bam = "${dataset}_final.bam.bai"
        """
        samtools index ${final_bam}
        """
}
