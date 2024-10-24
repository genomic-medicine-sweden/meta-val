include { SUBSET_BAM                    } from '../../modules/local/subset_bam'
include { SAMTOOLS_SORT                 } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_IDXSTATS             } from '../../modules/nf-core/samtools/idxstats/main'


workflow TAXID_BAM {
    take:
    bam
    bai
    accession2taxid

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    input_bam =  bam.combine( bai,by: 0 )
    SAMTOOLS_IDXSTATS( input_bam )
    ch_accession = SAMTOOLS_IDXSTATS.out.idxstats
        .map { it[1] }
        .splitCsv( header: false,sep:"\t" )
        .filter { it -> it[0]!= "*" }

    ch_versions.mix( SAMTOOLS_IDXSTATS.out.versions.first() )

    // Load accession2taxid.map
    ch_accession2taxidmap = accession2taxid.splitCsv( header: false,sep:"\t" )

    ch_accession_taxid = ch_accession2taxidmap
        .join( ch_accession )
        .filter { it -> it[3] != "0" }
        .map { [ it[0], it[1] ] }
        .groupTuple( by: 1 )

    ch_samtools_view = ch_accession_taxid
        .combine(input_bam)
        //.view()
        .map {accession_list, taxid, meta, bam, bam_index ->
            def new_meta = meta.clone()
            new_meta.taxid = taxid
            return [ new_meta, bam, bam_index, accession_list ]
        }
        .multiMap {
            meta, bam, bam_index, accession_list ->
                bam:  [meta, bam, bam_index]
                accession: accession_list.flatten()
        }

    SUBSET_BAM ( ch_samtools_view.bam, ch_samtools_view.accession )
    ch_versions      = ch_versions.mix( SUBSET_BAM.out.versions.first() )

    SAMTOOLS_SORT ( SUBSET_BAM.out.bam, [[],[]] )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )
    ch_versions      = ch_versions.mix( SAMTOOLS_INDEX.out.versions.first() )

    emit:
    accession       = ch_accession
    versions        = ch_versions
    taxid_bam       = SAMTOOLS_SORT.out.bam
    taxid_bam_bai   = SAMTOOLS_INDEX.out.bai

}
