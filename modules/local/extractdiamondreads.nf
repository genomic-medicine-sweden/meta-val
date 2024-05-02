process EXTRACTCDIAMONDREADS {

    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.8.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.8.0--h9ee0642_0':
        'biocontainers/seqkit:2.8.0--h9ee0642_0' }"

    input:
    val taxid
    tuple val (meta), path(tsv)
    tuple val (meta), path(fastq) // bowtie2/align *unmapped_{1,2}.fastq.gz

    output:
    tuple val(meta), path("*.fastq"), emit: extracted_diamond_reads
    path "versions.yml", emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    awk '\$2 == $taxid {print \$1}' $tsv > readID.txt
    if (${meta.single_end}) {
        seqkit grep -f readID.txt $fastq > ${prefix}_${taxid}.extracted_diamond_read.fastq
    } else {
        seqkit grep -f readID.txt ${fastq[0]} > ${prefix}_${taxid}.extracted_diamond_read1.fastq
        seqkit grep -f readID.txt ${fastq[1]} > ${prefix}_${taxid}.extracted_diamond_read2.fastq
    }

    rm readID.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}