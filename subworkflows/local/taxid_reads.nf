include { EXTRACT_VIRAL_TAXID as KRAKEN2_VIRAL_TAXID      } from '../../modules/local/extract_viral_taxid'
include { EXTRACT_VIRAL_TAXID as CENTRIFUGE_VIRAL_TAXID   } from '../../modules/local/extract_viral_taxid'
include { EXTRACT_VIRAL_TAXID as DIAMOND_VIRAL_TAXID      } from '../../modules/local/extract_viral_taxid'
include { KRAKENTOOLS_EXTRACTKRAKENREADS                  } from '../../modules/nf-core/krakentools/extractkrakenreads/main'
include { EXTRACTCENTRIFUGEREADS                          } from '../../modules/local/extractcentrifugereads'
include { EXTRACTCDIAMONDREADS                            } from '../../modules/local/extractdiamondreads'

workflow TAXID_READS {
    params.taxid

    take:
    reads                   // channel:   [mandatory] [ meta, reads ]
    kraken2_taxpasta        // channel:   [mandatory] [ meta, kraken2_taxpasta ]
    kraken2_result          // channel:   [mandatory] [ meta, kraken2_result ]
    kraken2_report          // channel:   [mandatory] [ meta, kraken2_report ]
    centrifuge_taxpasta     // channel:   [mandatory] [ meta, centrifuge_taxpasta ]
    centrifuge_result       // channel:   [mandatory] [ meta, centrifuge_result ]
    centrifuge_report       // channel:   [mandatory] [ meta, centrifuge_report ]
    diamond_taxpasta        // channel:   [mandatory] [ meta, diamond_taxpasta ]
    diamond_tsv             // channel:   [mandatory] [ meta, diamond_tsv ]

    main:
    ch_versions = Channel.empty()

    // extract kraken2 reads
    if ( params.extract_kraken2_reads ) {
        if ( params.taxid ) {
            kraken2_params_taxid = kraken2_report.map { meta, kraken2_report -> [ meta.subMap(meta.keySet() - 'tool'), kraken2_report ] }
                .combine ( Channel.of( params.taxid.split(" ") ) )
                .combine( kraken2_result, by: 0 )
                .combine( reads, by: 0)
                .multiMap { meta, kraken2_report, taxid, kraken2_result, reads  ->
                    taxid: taxid
                    kraken2_result: [ meta, kraken2_result ]
                    reads: [ meta, reads ]
                    kraken2_report: [ meta, kraken2_report ]
                    }

            KRAKENTOOLS_EXTRACTKRAKENREADS(
                kraken2_params_taxid.taxid,
                kraken2_params_taxid.kraken2_result,
                kraken2_params_taxid.reads,
                kraken2_params_taxid.kraken2_report
            )
            ch_versions            = ch_versions.mix( KRAKENTOOLS_EXTRACTKRAKENREADS.out.versions.first() )

        } else {
            kraken2_taxids = KRAKEN2_VIRAL_TAXID( kraken2_taxpasta, kraken2_report )
            kraken2_combined_input = kraken2_taxids.viral_taxid
                .map { meta, taxid -> [ meta.subMap( meta.keySet() - 'tool' ), taxid ] }
                .splitText()
                .combine( kraken2_result, by:0 )
                .combine( reads, by:0 )
                .combine( kraken2_report.map { meta, kraken2_report -> [ meta.subMap(meta.keySet() - 'tool'), kraken2_report ]}, by:0 )
                .multiMap { meta, taxid, kraken2_result, reads, kraken2_report ->
                    taxid: taxid.trim()
                    kraken2_result: [ meta, kraken2_result ]
                    reads: [ meta, reads ]
                    kraken2_report: [ meta, kraken2_report ]
                }

                KRAKENTOOLS_EXTRACTKRAKENREADS(
                kraken2_combined_input.taxid,
                kraken2_combined_input.kraken2_result,
                kraken2_combined_input.reads,
                kraken2_combined_input.kraken2_report
            )
            ch_versions            = ch_versions.mix( KRAKEN2_VIRAL_TAXID.out.versions.first(), KRAKENTOOLS_EXTRACTKRAKENREADS.out.versions.first() )
        }
    }

    // extract centrifuge reads
    if ( params.extract_centrifuge_reads ) {
        if ( params.taxid ) {
            centrifuge_params_taxid = centrifuge_result
                .combine( Channel.of( params.taxid.split(" ") ) )
                .combine( reads, by: 0 )
                .multiMap { meta, centrifuge_result, taxid, reads ->
                    taxid: taxid
                    centrifuge_result: [ meta, centrifuge_result ]
                    reads: [ meta, reads ]
                    }

            EXTRACTCENTRIFUGEREADS(
                centrifuge_params_taxid.taxid,
                centrifuge_params_taxid.centrifuge_result,
                centrifuge_params_taxid.reads
            )
            ch_versions            = ch_versions.mix( EXTRACTCENTRIFUGEREADS.out.versions )

        } else {
            centrifuge_taxids = CENTRIFUGE_VIRAL_TAXID( centrifuge_taxpasta, centrifuge_report )
            centrifuge_combined_input = centrifuge_taxids.viral_taxid
                .map { meta, taxid -> [ meta.subMap( meta.keySet() - 'tool' ), taxid ] }
                .splitText()
                .combine( centrifuge_result, by:0 )
                .combine( reads, by:0 )
                .multiMap { meta, taxid, centrifuge_result, reads ->
                    taxid: taxid.trim()
                    centrifuge_result: [ meta, centrifuge_result ]
                    reads: [ meta, reads ]
                }

            EXTRACTCENTRIFUGEREADS(
                centrifuge_combined_input.taxid,
                centrifuge_combined_input.centrifuge_result,
                centrifuge_combined_input.reads,
            )
            ch_versions            = ch_versions.mix( CENTRIFUGE_VIRAL_TAXID.out.versions.first(), EXTRACTCENTRIFUGEREADS.out.versions )
        }
    }

    // extract diamond reads
    if ( params.extract_diamond_reads ) {
        if ( params.taxid ) {
            diamond_params_taxid = diamond_tsv.map { meta, diamond_tsv -> [meta.subMap( meta.keySet() - 'tool' ), diamond_tsv ] }
                .combine( Channel.of( params.taxid.split(" ") ))
                .combine( reads, by:0)
                .multiMap { meta, diamond_tsv, taxid, reads ->
                    taxid: taxid
                    diamond_tsv: [ meta, diamond_tsv ]
                    reads: [ meta, reads ]
                    }

            EXTRACTCDIAMONDREADS(
                diamond_params_taxid.taxid,
                diamond_params_taxid.diamond_tsv,
                diamond_params_taxid.reads
            )
            ch_versions            = ch_versions.mix( EXTRACTCDIAMONDREADS.out.versions )

        } else {
            diamond_taxids = DIAMOND_VIRAL_TAXID( diamond_taxpasta, diamond_tsv )
            diamond_combined_input = diamond_taxids.viral_taxid
                .map { meta, taxid -> [ meta.subMap( meta.keySet() - 'tool' ), taxid ] }
                .splitText()
                .combine( diamond_tsv.map{ meta, diamond_tsv -> [meta.subMap( meta.keySet() - 'tool' ), diamond_tsv ] }, by:0 )
                .combine( reads, by:0 )
                .multiMap { meta, taxid, diamond, reads ->
                    taxid: taxid.trim()
                    diamond_tsv: [ meta, diamond ]
                    reads: [ meta, reads ]
                }

            EXTRACTCDIAMONDREADS(
                diamond_combined_input.taxid,
                diamond_combined_input.diamond_tsv,
                diamond_combined_input.reads,
            )
            ch_versions            = ch_versions.mix( DIAMOND_VIRAL_TAXID.out.versions.first(), EXTRACTCDIAMONDREADS.out.versions )
        }
    }

    emit:
    versions        = ch_versions          // channel: [ versions.yml ]
}
