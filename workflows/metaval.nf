/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { FASTQC                          } from '../modules/nf-core/fastqc/main'
include { MULTIQC                         } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                } from 'plugin/nf-validation'
include { paramsSummaryMultiqc            } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML          } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText          } from '../subworkflows/local/utils_nfcore_metaval_pipeline'
include {EXTRACT_VIRAL_TAXID              } from '../modules/local/extract_viral_taxid'
include { KRAKENTOOLS_EXTRACTKRAKENREADS  } from '../modules/nf-core/krakentools/extractkrakenreads/main'
include { EXTRACTCENTRIFUGEREADS          } from '../modules/local/extractcentrifugereads'
include { EXTRACTCDIAMONDREADS            } from '../modules/local/extractdiamondreads'


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

    // Create channels
    ch_input = ch_samplesheet.branch { meta, fastq_1, fastq_2, kraken2_report, kraken2_result, kraken2_taxpasta, centrifuge_report, centrifuge_result, centrifuge_taxpasta, diamond, diamond_taxpasta ->

        // Define single_end based on the conditions
        meta.single_end = ( fastq_1 && !fastq_2 )

        // reads channels
        short_reads: meta.instrument_platform != 'OXFORD_NANOPORE'
            return [ meta, fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ] ]

        long_reads: meta.instrument_platform == 'OXFORD_NANOPORE'
            return [ meta, [ fastq_1 ] ]
    }

    // extract kraken2 reads
    ch_input_kraken2 = ch_samplesheet.multiMap { meta, fastq_1, fastq_2, kraken2_report, kraken2_result, kraken2_taxpasta, centrifuge_report, centrifuge_result, centrifuge_taxpasta, diamond, diamond_taxpasta ->
        meta.single_end = ( fastq_1 && !fastq_2 )
        def new_meta = meta + [ tool: "kraken2" ]
        kraken2_taxpasta: [ new_meta, kraken2_taxpasta ]
        kraken2_result: [ new_meta, kraken2_result ]
        kraken2_report: [ new_meta, kraken2_report ]
        reads:[ new_meta, fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ] ]
    }

    if ( params.extract_kraken2_reads ) {
        if ( params.taxid ) {
            KRAKENTOOLS_EXTRACTKRAKENREADS(
                params.taxid,
                ch_input_kraken2.kraken2_result,
                ch_input_kraken2.reads,
                ch_input_kraken2.kraken2_report
            )
            ch_versions            = ch_versions.mix( KRAKENTOOLS_EXTRACTKRAKENREADS.out.versions.first() )

        } else {
            kraken2_taxids = EXTRACT_VIRAL_TAXID( ch_input_kraken2.kraken2_taxpasta, ch_input_kraken2.kraken2_report )
            combined_input = kraken2_taxids.viral_taxid
                .splitText()
                .combine( ch_input_kraken2.kraken2_result, by:0 )
                .combine( ch_input_kraken2.reads, by:0 )
                .combine( ch_input_kraken2.kraken2_report, by:0 )

            ch_combined_input = combined_input.multiMap { meta,taxid,kraken2_result,reads,kraken2_report  ->
                taxid: taxid.trim()
                kraken2_result: [ meta, kraken2_result ]
                reads: [ meta, reads ]
                kraken2_report: [ meta, kraken2_report ]
            }
            KRAKENTOOLS_EXTRACTKRAKENREADS(
                ch_combined_input.taxid,
                ch_combined_input.kraken2_result,
                ch_combined_input.reads,
                ch_combined_input.kraken2_report
            )
            ch_versions            = ch_versions.mix( kraken2_taxids.versions.first(), KRAKENTOOLS_EXTRACTKRAKENREADS.out.versions.first() )
        }
    }

    // extract centrifuge reads
    ch_input_centrifuge = ch_samplesheet.multiMap { meta, fastq_1, fastq_2, kraken2_report, kraken2_result, kraken2_taxpasta, centrifuge_report, centrifuge_result, centrifuge_taxpasta, diamond, diamond_taxpasta ->
        meta.single_end = ( fastq_1 && !fastq_2 )
        def new_meta = meta + [ tool: "centrifuge" ]
        centrifuge_taxpasta: [ new_meta, centrifuge_taxpasta ]
        centrifuge_result: [ new_meta, centrifuge_result ]
        centrifuge_report: [ new_meta, centrifuge_report ]
        reads:[ new_meta, fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ] ]
    }

    if ( params.extract_centrifuge_reads ) {
        if ( params.taxid ) {
            EXTRACTCENTRIFUGEREADS(
                params.taxid,
                ch_input_centrifuge.centrifuge_result,
                ch_input_centrifuge.reads
            )
            ch_versions            = ch_versions.mix( EXTRACTCENTRIFUGEREADS.out.versions )

        } else {
            centrifuge_taxids = EXTRACT_VIRAL_TAXID( ch_input_centrifuge.centrifuge_taxpasta, ch_input_centrifuge.centrifuge_report )
            combined_input = centrifuge_taxids.viral_taxid
                .splitText()
                .combine( ch_input_centrifuge.centrifuge_result, by:0 )
                .combine( ch_input_centrifuge.reads, by:0 )

            ch_combined_input = combined_input.multiMap { meta,taxid,centrifuge_result,reads  ->
                taxid: taxid.trim()
                centrifuge_result: [ meta, centrifuge_result ]
                reads: [ meta, reads ]
            }

            EXTRACTCENTRIFUGEREADS(
                ch_combined_input.taxid,
                ch_combined_input.centrifuge_result,
                ch_combined_input.reads,
            )
            ch_versions            = ch_versions.mix( centrifuge_taxids.versions.first(), EXTRACTCENTRIFUGEREADS.out.versions )
        }
    }

    // extract diamond reads
    ch_input_diamond = ch_samplesheet.multiMap { meta, fastq_1, fastq_2, kraken2_report, kraken2_result, kraken2_taxpasta, centrifuge_report, centrifuge_result, centrifuge_taxpasta, diamond, diamond_taxpasta ->
        meta.single_end = ( fastq_1 && !fastq_2 )
        def new_meta = meta + [ tool: "diamond" ]
        diamond_taxpasta: [ new_meta, diamond_taxpasta ]
        diamond_tsv: [ new_meta, diamond ]
        reads:[ new_meta, fastq_2 ? [ fastq_1, fastq_2 ] : [ fastq_1 ] ]
    }

    if ( params.extract_diamond_reads ) {
        if ( params.taxid ) {
            EXTRACTCDIAMONDREADS(
                params.taxid,
                ch_input_diamond.diamond_tsv,
                ch_input_diamond.reads
            )
            ch_versions            = ch_versions.mix( EXTRACTCDIAMONDREADS.out.versions )

        } else {
            diamond_taxids = EXTRACT_VIRAL_TAXID( ch_input_diamond.diamond_taxpasta, ch_input_diamond.diamond_tsv )
            combined_input = diamond_taxids.viral_taxid
                .splitText()
                .combine( ch_input_diamond.diamond_tsv, by:0 )
                .combine( ch_input_diamond.reads, by:0 )

            ch_combined_input = combined_input.multiMap { meta,taxid,diamond,reads  ->
                taxid: taxid.trim()
                diamond_tsv: [ meta, diamond ]
                reads: [ meta, reads ]
            }

            EXTRACTCDIAMONDREADS(
                ch_combined_input.taxid,
                ch_combined_input.diamond_tsv,
                ch_combined_input.reads,
            )
            ch_versions            = ch_versions.mix( diamond_taxids.versions.first(), EXTRACTCDIAMONDREADS.out.versions )
        }
    }

    //
    // MODULE: Run FastQC
    //
    //FASTQC (
    //    ch_samplesheet
    //)
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    //ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_pipeline_software_mqc_versions.yml', sort: true, newLine: true)
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: false))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
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
