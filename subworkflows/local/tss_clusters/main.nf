//
// CTSS extraction, consensus TSS clusters, and grouping of read pairs by TSS cluster
//

include { BAM_TO_CTSS                     } from '../../../modules/local/bam_to_ctss/main'
include { POOL_CTSS                       } from '../../../modules/local/pool_ctss/main'
include { PARACLU                         } from '../../../modules/nf-core/paraclu/main'
include { ANNOTATE_TSS_CLUSTERS           } from '../../../modules/local/annotate_tss_clusters/main'
include { ASSIGN_PAIRS_TO_TSS             } from '../../../modules/local/assign_pairs_to_tss/main'
include { MERGE_COUNTS as MERGE_TSS_COUNTS } from '../../../modules/local/merge_counts/main'
include { MERGE_CAGESCAN_CLUSTERS         } from '../../../modules/local/merge_cagescan_clusters/main'
include { UCSC_BEDGRAPHTOBIGWIG           } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow TSS_CLUSTERS {

    take:
    ch_bam_bai  // channel: [ val(meta), path(bam), path(bai) ]
    ch_sizes    // channel: [ val(meta), path(chrom.sizes) ]
    min_cluster // integer: minimum total CTSS count of a paraclu cluster

    main:

    //
    // READ1 5' ends -> CTSS per sample
    //
    BAM_TO_CTSS(ch_bam_bai)

    //
    // Per-strand CTSS bigWigs (skip empty strands)
    //
    def ch_bedgraph = BAM_TO_CTSS.out.bedgraph
        .transpose()
        .filter { _meta, bedgraph -> bedgraph.size() > 0 }
        .map { meta, bedgraph -> [meta + [id: bedgraph.name - '.bedgraph'], bedgraph] }
    UCSC_BEDGRAPHTOBIGWIG(ch_bedgraph, ch_sizes.map { _meta, sizes -> sizes })

    //
    // Pool CTSS across samples -> consensus TSS clusters
    //
    POOL_CTSS(
        BAM_TO_CTSS.out.ctss
            .map { _meta, ctss -> ctss }
            .collect()
            .map { ctss -> [[id: 'all_samples'], ctss] }
    )
    PARACLU(POOL_CTSS.out.ctss, min_cluster)
    ANNOTATE_TSS_CLUSTERS(
        PARACLU.out.bed.combine(POOL_CTSS.out.ctss).map { meta, bed, _ctss_meta, ctss -> [meta, bed, ctss] }
    )

    // Single-item channels reused by every sample
    def ch_clusters_bed = ANNOTATE_TSS_CLUSTERS.out.bed.first()
    def ch_clusters_tsv = ANNOTATE_TSS_CLUSTERS.out.tsv.map { _meta, tsv -> tsv }.first()

    //
    // Keep read pairs whose READ1 starts in a TSS cluster
    //
    ASSIGN_PAIRS_TO_TSS(ch_bam_bai, ch_clusters_bed)

    MERGE_TSS_COUNTS(
        ASSIGN_PAIRS_TO_TSS.out.counts
            .map { _meta, counts -> counts }
            .collect()
            .map { counts -> [[id: 'all_samples'], counts] },
        ch_clusters_tsv,
        'tss'
    )

    //
    // FANTOM5-style CAGEscan meta-clusters across samples
    //
    MERGE_CAGESCAN_CLUSTERS(
        ASSIGN_PAIRS_TO_TSS.out.cagescan_clusters
            .map { _meta, bed -> bed }
            .collect()
            .map { beds -> [[id: 'all_samples'], beds] }
    )

    def ch_multiqc_files = BAM_TO_CTSS.out.stats
        .map { _meta, stats -> stats }
        .collectFile(name: 'ctss_stats_mqc.tsv', keepHeader: true, skip: 5, sort: true)
        .mix(
            ASSIGN_PAIRS_TO_TSS.out.stats
                .map { _meta, stats -> stats }
                .collectFile(name: 'tss_assign_stats_mqc.tsv', keepHeader: true, skip: 5, sort: true)
        )

    emit:
    ctss              = BAM_TO_CTSS.out.ctss            // channel: [ val(meta), path(ctss.bed) ]
    softclip          = BAM_TO_CTSS.out.softclip        // channel: [ val(meta), path(tsv) ]
    clusters_bed      = ch_clusters_bed                 // channel: [ val(meta), path(tss_clusters.bed) ]
    clusters_tsv      = ch_clusters_tsv                 // channel: path(tss_clusters.tsv)
    tss_bam           = ASSIGN_PAIRS_TO_TSS.out.bam     // channel: [ val(meta), path(bam), path(bai) ]
    counts            = MERGE_TSS_COUNTS.out.tsv        // channel: [ val(meta), path(tsv) ]
    cagescan_clusters = MERGE_CAGESCAN_CLUSTERS.out.bed // channel: [ val(meta), path(bed12) ]
    multiqc_files     = ch_multiqc_files                // channel: path(*_mqc.tsv)
}
