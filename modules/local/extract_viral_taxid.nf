process EXTRACT_VIRAL_TAXID {

    tag "$meta.id"
    label 'process_low'

    conda "bioconda::seqkit=2.4.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqkit:2.4.0--h9ee0642_0':
        'biocontainers/seqkit:2.4.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(taxpasta_standardised_profile)

    output:
    tuple val(meta), path("*viral_taxids.tsv"), optional:true, emit: viral_taxid
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.tool}"

    """
    if grep -qi "virus" $taxpasta_standardised_profile; then
        grep -i "virus" $taxpasta_standardised_profile | cut -f 1 > ${prefix}_viral_taxids.tsv
    else
        echo "No viral taxids found." > "no_viral_taxid.txt"
    fi


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """
}
