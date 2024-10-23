/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Extract reads of taxIDs
include { EXTRACT_VIRAL_TAXID                       } from '../modules/local/extract_viral_taxid'
include { KRAKENTOOLS_EXTRACTKRAKENREADS            } from '../modules/nf-core/krakentools/extractkrakenreads/main'
include { EXTRACTCENTRIFUGEREADS                    } from '../modules/local/extractcentrifugereads'
include { EXTRACTCDIAMONDREADS                      } from '../modules/local/extractdiamondreads'
include { TAXID_READS                               } from '../subworkflows/local/taxid_reads'

// Maping subworkflow
include { BOWTIE2_BUILD as BOWTIE2_BUILD_PATHOGEN   } from '../modules/nf-core/bowtie2/build/main'
include { FASTQ_ALIGN_BOWTIE2                       } from '../subworkflows/nf-core/fastq_align_bowtie2/main'
include { LONGREAD_SCREENPATHOGEN                   } from '../subworkflows/local/longread_screenpathogen'

// Summary subworkflow
include { FASTQC                                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                          } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                    } from '../subworkflows/local/utils_nfcore_metaval_pipeline'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METAVAL {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Create input channels
    ch_input = ch_samplesheet.branch { meta, fastq_1, fastq_2, kraken2_report, kraken2_result, kraken2_taxpasta, centrifuge_report, centrifuge_result, centrifuge_taxpasta, diamond, diamond_taxpasta ->

        // Define single_end based on the conditions
        meta.single_end = ( fastq_1 && !fastq_2 )

        // reads channels
        short_reads: meta.instrument_platform != 'OXFORD_NANOPORE'
            return [ meta, fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ] ]

        long_reads: meta.instrument_platform == 'OXFORD_NANOPORE'
            return [ meta, [ fastq_1 ] ]
    }

    // channels for extracting kraken2 reads
    ch_extract_reads = ch_samplesheet.multiMap { meta, fastq_1, fastq_2, kraken2_report, kraken2_result, kraken2_taxpasta, centrifuge_report, centrifuge_result, centrifuge_taxpasta, diamond, diamond_taxpasta ->
        meta.single_end = ( fastq_1 && !fastq_2 )
        kraken2_taxpasta: [ meta + [ tool: "kraken2" ], kraken2_taxpasta ]
        kraken2_report: [ meta + [ tool: "kraken2" ], kraken2_report ]
        kraken2_result: [ meta, kraken2_result ]
        reads:[ meta, fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ] ]
        centrifuge_taxpasta: [ meta + [ tool: "centrifuge" ], centrifuge_taxpasta ]
        centrifuge_report: [ meta + [ tool: "centrifuge" ], centrifuge_report ]
        centrifuge_result: [ meta, centrifuge_result ]
        diamond_taxpasta: [ meta + [ tool: "diamond" ], diamond_taxpasta ]
        diamond_tsv: [ meta + [ tool: "diamond" ], diamond ]
    }

    /*
        SUBWORKFLOW: TAXID_READS
    */
    if ( params.perform_extract_reads ) {
        TAXID_READS (
        ch_extract_reads.reads,
        ch_extract_reads.kraken2_taxpasta,
        ch_extract_reads.kraken2_result,
        ch_extract_reads.kraken2_report,
        ch_extract_reads.centrifuge_taxpasta,
        ch_extract_reads.centrifuge_result,
        ch_extract_reads.centrifuge_report,
        ch_extract_reads.diamond_taxpasta,
        ch_extract_reads.diamond_tsv,
        )
        ch_versions            = ch_versions.mix( TAXID_READS.out.versions )
    }

    /*
        SUBWORKFLOW: Screen pathogens
    */
    ch_reference = Channel.fromPath( params.pathogens_genomes, checkIfExists: true)
        .map{ file -> [ [ id: file.baseName ], file ] }
    // Short reads
    if ( params.perform_screen_pathogens ) {
        BOWTIE2_BUILD_PATHOGEN ( ch_reference )
        ch_versions      = ch_versions.mix( BOWTIE2_BUILD_PATHOGEN.out.versions )
        FASTQ_ALIGN_BOWTIE2 (
            ch_input.short_reads,                              // ch_reads
            BOWTIE2_BUILD_PATHOGEN.out.index,                  // ch_index
            false,                                             // save unaligned
            false,                                             // sort bam
            ch_reference                                       // ch_fasta
        )
        ch_versions = ch_versions.mix( FASTQ_ALIGN_BOWTIE2.out.versions )
    // Long reads
        LONGREAD_SCREENPATHOGEN ( ch_input.long_reads, ch_reference )
        ch_versions = ch_versions.mix( LONGREAD_SCREENPATHOGEN.out.versions )
    }

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_extract_reads.reads
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  ''  + 'pipeline_software_' +  'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml',sort: true))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
