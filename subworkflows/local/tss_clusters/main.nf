//
// CTSS extraction, "broad" and "sharp" TSS cluster sets, and grouping of read pairs
// by TSS cluster (once per cluster set)
//

include { BAM_TO_CTSS                     } from '../../../modules/local/bam_to_ctss/main'
include { POOL_CTSS                       } from '../../../modules/local/pool_ctss/main'
include { PARACLU                         } from '../../../modules/nf-core/paraclu/main'
include { ANNOTATE_TSS_CLUSTERS           } from '../../../modules/local/annotate_tss_clusters/main'
include { ASSIGN_PAIRS_TO_TSS             } from '../../../modules/local/assign_pairs_to_tss/main'
include { MERGE_COUNTS as MERGE_TSS_COUNTS } from '../../../modules/local/merge_counts/main'
include { MERGE_CAGESCAN_CLUSTERS         } from '../../../modules/local/merge_cagescan_clusters/main'
include { UCSC_BEDGRAPHTOBIGWIG           } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { RECLU                           } from '../reclu'

workflow TSS_CLUSTERS {

    take:
    ch_bam_bai          // channel: [ val(meta), path(bam), path(bai) ]; meta.group = replicate group
    ch_sizes            // channel: [ val(meta), path(chrom.sizes) ]
    min_cluster         // integer: minimum total CTSS count of a broad (pooled paraclu) cluster
    run_sharp           // boolean: also build sharp (RECLU) clusters
    reclu_paraclu_min   // integer: paraclu minValue for the RECLU cluster hierarchy

    main:

    //
    // READ1 5' ends -> CTSS per sample
    //
    BAM_TO_CTSS(ch_bam_bai)

    //
    // CTSS pooled across samples
    //
    POOL_CTSS(
        BAM_TO_CTSS.out.ctss
            .map { _meta, ctss -> ctss }
            .collect()
            .map { ctss -> [[id: 'all_samples'], ctss] }
    )

    //
    // Per-strand CTSS bigWigs, per sample and pooled (skip empty strands)
    //
    def ch_bedgraph = BAM_TO_CTSS.out.bedgraph
        .mix(POOL_CTSS.out.bedgraph)
        .transpose()
        .filter { _meta, bedgraph -> bedgraph.size() > 0 }
        .map { meta, bedgraph -> [meta + [id: bedgraph.name - '.bedgraph'], bedgraph] }
    UCSC_BEDGRAPHTOBIGWIG(ch_bedgraph, ch_sizes.map { _meta, sizes -> sizes })

    //
    // Broad clusters: paraclu on the pooled CTSS
    //
    PARACLU(POOL_CTSS.out.ctss, min_cluster)
    ANNOTATE_TSS_CLUSTERS(
        PARACLU.out.bed.combine(POOL_CTSS.out.ctss).map { meta, bed, _ctss_meta, ctss -> [meta + [cluster_set: 'broad'], bed, ctss] }
    )
    def ch_sets = ANNOTATE_TSS_CLUSTERS.out.bed
        .join(ANNOTATE_TSS_CLUSTERS.out.tsv)
        .map { _meta, bed, tsv -> ['broad', bed, tsv] }

    //
    // Sharp clusters: RECLU reproducible clusters (needs replicates)
    //
    def ch_reclu_stats = channel.empty()
    if (run_sharp) {
        RECLU(BAM_TO_CTSS.out.ctss, POOL_CTSS.out.ctss, reclu_paraclu_min)
        ch_sets = ch_sets.mix(RECLU.out.clusters.map { _meta, bed, tsv -> ['sharp', bed, tsv] })
        ch_reclu_stats = RECLU.out.stats.map { _meta, stats -> stats }
    }

    //
    // Keep read pairs whose READ1 starts in a TSS cluster, for every sample x cluster set
    //
    def ch_assign = ch_bam_bai
        .combine(ch_sets)
        .multiMap { meta, bam, bai, set, bed, _tsv ->
            bam:      [meta + [cluster_set: set], bam, bai]
            clusters: [[id: "all_samples.${set}"], bed]
        }
    ASSIGN_PAIRS_TO_TSS(ch_assign.bam, ch_assign.clusters)

    def ch_counts = ASSIGN_PAIRS_TO_TSS.out.counts
        .map { meta, counts -> [meta.cluster_set, counts] }
        .groupTuple()
        .join(ch_sets.map { set, _bed, tsv -> [set, tsv] })
        .multiMap { set, counts, tsv ->
            files:    [[id: 'all_samples', cluster_set: set], counts]
            clusters: tsv
        }
    MERGE_TSS_COUNTS(ch_counts.files, ch_counts.clusters, 'tss')

    //
    // FANTOM5-style CAGEscan meta-clusters across samples, per cluster set
    //
    MERGE_CAGESCAN_CLUSTERS(
        ASSIGN_PAIRS_TO_TSS.out.cagescan_clusters
            .map { meta, bed -> [meta.cluster_set, bed] }
            .groupTuple()
            .map { set, beds -> [[id: 'all_samples', cluster_set: set], beds] }
    )

    def ch_multiqc_files = BAM_TO_CTSS.out.stats
        .map { _meta, stats -> stats }
        .collectFile(name: 'ctss_stats_mqc.tsv', keepHeader: true, skip: 5, sort: true)
        .mix(
            ASSIGN_PAIRS_TO_TSS.out.stats
                .map { _meta, stats -> stats }
                .collectFile(name: 'tss_assign_stats_mqc.tsv', keepHeader: true, skip: 5, sort: true)
        )
        .mix(ch_reclu_stats)

    emit:
    ctss              = BAM_TO_CTSS.out.ctss                            // channel: [ val(meta), path(ctss.bed) ]
    softclip          = BAM_TO_CTSS.out.softclip                        // channel: [ val(meta), path(tsv) ]
    ctss_stats        = BAM_TO_CTSS.out.stats                           // channel: [ val(meta), path(ctss_stats_mqc.tsv) ]
    clusters_tsv      = ch_sets.map { set, _bed, tsv -> [set, tsv] }    // channel: [ val(cluster_set), path(tss_clusters.tsv) ]
    tss_bam           = ASSIGN_PAIRS_TO_TSS.out.bam                     // channel: [ val(meta), path(bam), path(bai) ]; meta.cluster_set
    counts            = MERGE_TSS_COUNTS.out.tsv                        // channel: [ val(meta), path(tsv) ]
    cagescan_clusters = MERGE_CAGESCAN_CLUSTERS.out.bed                 // channel: [ val(meta), path(bed12) ]
    multiqc_files     = ch_multiqc_files                                // channel: path(*_mqc.tsv)
}
