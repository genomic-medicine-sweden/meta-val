/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Extract reads of taxIDs
include { KRAKENTOOLS_EXTRACTKRAKENREADS            } from '../modules/nf-core/krakentools/extractkrakenreads/main'
include { EXTRACTCENTRIFUGEREADS                    } from '../modules/local/extractcentrifugereads'
include { EXTRACTCDIAMONDREADS                      } from '../modules/local/extractdiamondreads'
include { TAXID_READS                               } from '../subworkflows/local/taxid_reads'
include { RM_EMPTY_FASTQ as RM_EMPTY_KRAKEN2        } from '../modules/local/rm_empty_fastq'
include { RM_EMPTY_FASTQ as RM_EMPTY_CENTRIFUGE     } from '../modules/local/rm_empty_fastq'
include { RM_EMPTY_FASTQ as RM_EMPTY_DIAMOND        } from '../modules/local/rm_empty_fastq'

// De novo for extracted taxIDs reads
include { SPADES                                    } from '../modules/nf-core/spades/main'
include { FLYE                                      } from '../modules/nf-core/flye/main'

// Maping subworkflow
include { BOWTIE2_BUILD as BOWTIE2_BUILD_PATHOGEN   } from '../modules/nf-core/bowtie2/build/main'
include { FASTQ_ALIGN_BOWTIE2                       } from '../subworkflows/nf-core/fastq_align_bowtie2/main'
include { LONGREAD_SCREENPATHOGEN                   } from '../subworkflows/local/longread_screenpathogen'

// Summary subworkflow
include { FASTQC                                    } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                          } from 'plugin/nf-validation'
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

    // Verify whether the taxonomic IDs identified by classification are true or false positives.
    if ( params.perform_extract_reads ) {

        // SUBWORKFLOW: TAXID_READS - extract reads
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
        // Remove empty fastq files produced by extracting reads for user defined taxIDs
        if (params.extract_kraken2_reads && params.taxid) {
            RM_EMPTY_KRAKEN2(file("${params.outdir}/extracted_reads/kraken2"))
        }
        if (params.extract_centrifuge_reads && params.taxid) {
            RM_EMPTY_CENTRIFUGE(file("${params.outdir}/extracted_reads/centrifuge"))
        }
        if (params.extract_diamond_reads && params.taxid) {
            RM_EMPTY_DIAMOND(file("${params.outdir}/extracted_reads/diamond"))
        }

        // SUBWORKFLOW: DE NOVO

        // Filter out empty FASTQ files
        ch_taxid_reads_result = TAXID_READS.out.reads
            .branch {
                non_empty: it[0].single_end ? it[1].size() > 0 : it[1][0].size() > 0 || it[1][1].size() >0
                empty: true
            }
        ch_taxid_reads_result.non_empty.set { ch_taxid_reads }

        // Skip the de-novo assembly if the number of reads is lower than params.min_read_counts
        ch_taxid_reads
            .branch {
                failed: it[0].single_end ? it[1].countFastq() < params.min_read_counts : it[1][0].countFastq() < params.min_read_counts || it[1][1].countFastq() < params.min_read_counts
                passed: true
            }
            .set { ch_taxid_reads_result }
        ch_taxid_reads_result.passed.set { ch_taxid_reads_result_passed }
        //Prepare reads for de-novo assembly
        ch_taxid_reads_result_passed
            .branch { meta, reads ->
                shortreads_spades: meta.instrument_platform != 'OXFORD_NANOPORE'
                    return [ meta, reads, [], [] ]
                longreads_denovo: meta.instrument_platform == 'OXFORD_NANOPORE'
                    return [ meta, reads ]
            }
            .set { ch_denovo_input }

        // short reads de novo assembly
        if ( params.perform_shortread_denovo ) {
            SPADES( ch_denovo_input.shortreads_spades, [], [] )
            ch_versions             = ch_versions.mix( SPADES.out.versions.first() )
        }
        // long reads de novo assembly
        if ( params.perform_longread_denovo ) {
            FLYE( ch_denovo_input.longreads_denovo, params.flye_mode )
            ch_versions             = ch_versions.mix( FLYE.out.versions.first() )
        }

        // Warning message for samples that failed to run de novo assembly due to an insufficient number of reads
        ch_taxid_reads_result.failed
            .map { meta, reads -> [ meta.id ] }
            .collect()
            .subscribe {
                samples = it.join("\n")
                log.warn "The following samples skipped de novo assembly due to too few reads (<$params.min_read_counts). Run BLASTx/BLASTn directly with all reads.: \n$samples"
            }

        // BLAST

        }

    // Screen pathogens
    ch_reference = Channel.fromPath( params.pathogens_genomes, checkIfExists: true)
        .map{ file -> [ [ id: file.baseName ], file ] }
    if ( params.perform_screen_pathogens ) {
        // Map short reads to the pathogens genome
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
        // Map long reads to the pathogens genome
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
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

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
