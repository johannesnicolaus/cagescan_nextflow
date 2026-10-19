//
// CTSS extraction, "broad" and "sharp" TSS cluster sets, and grouping of read pairs
// by TSS cluster (once per cluster set)
//

include { BAM_TO_CTSS                     } from '../../../modules/local/bam_to_ctss/main'
include { CTSS_MOTIFS                     } from '../../../modules/local/ctss_motifs/main'
include { POOL_CTSS                       } from '../../../modules/local/pool_ctss/main'
include { POOL_CTSS as POOL_CTSS_BY_SL     } from '../../../modules/local/pool_ctss/main'
include { PARACLU                         } from '../../../modules/nf-core/paraclu/main'
include { ANNOTATE_TSS_CLUSTERS           } from '../../../modules/local/annotate_tss_clusters/main'
include { CAGEFIGHTR_CLUSTERS             } from '../../../modules/local/cagefightr_clusters/main'
include { COMBINE_TSS_CLUSTERS            } from '../../../modules/local/combine_tss_clusters/main'
include { ASSIGN_PAIRS_TO_TSS             } from '../../../modules/local/assign_pairs_to_tss/main'
include { MERGE_COUNTS as MERGE_TSS_COUNTS } from '../../../modules/local/merge_counts/main'
include { MERGE_CAGESCAN_CLUSTERS         } from '../../../modules/local/merge_cagescan_clusters/main'
include { UCSC_BEDGRAPHTOBIGWIG           } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { RECLU                           } from '../reclu'

workflow TSS_CLUSTERS {

    take:
    ch_bam_bai          // channel: [ val(meta), path(bam), path(bai) ]; meta.group = replicate group
    ch_sizes            // channel: [ val(meta), path(chrom.sizes) ]
    cluster_method      // string:  broad cluster caller, 'cagefightr' or 'paraclu'
    min_cluster         // integer: minimum pooled CTSS count of a broad cluster with --cluster_method paraclu
    run_sharp           // boolean: also build sharp (RECLU) clusters
    reclu_paraclu_min   // integer: paraclu minValue for the RECLU cluster hierarchy
    sl_separate         // boolean: cluster spliced-leader and other CTSS separately (broad set; needs --sl_sequence)
    has_sl              // boolean: --sl_sequence given (BAM_TO_CTSS writes SL and non-SL CTSS)
    ch_fasta_fai        // channel: [ val(meta), path(fasta), path(fai) ]

    main:

    //
    // READ1 5' ends -> CTSS per sample
    //
    BAM_TO_CTSS(ch_bam_bai)

    //
    // Sequence at TSSs (initiator) and trans-splice sites (acceptor), for MultiQC
    //
    CTSS_MOTIFS(
        has_sl
            ? BAM_TO_CTSS.out.nosl_ctss.join(BAM_TO_CTSS.out.sl_ctss, failOnDuplicate: true, failOnMismatch: true)
            : BAM_TO_CTSS.out.ctss.map { meta, ctss -> [meta, ctss, []] },
        ch_fasta_fai
    )

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
    // Broad clusters, with CAGEfightR (default: TPM, CTSS support, slice-reduce, TC expression filter) or paraclu on
    // the pooled CTSS. With a spliced leader, trans-splice sites (SL READ1s) and the other TSSs are clustered
    // separately, so that a TSS next to a trans-splice site keeps its own peak.
    //
    def ch_broad_ctss = POOL_CTSS.out.ctss
    if (sl_separate) {
        POOL_CTSS_BY_SL(
            BAM_TO_CTSS.out.sl_ctss.map { _meta, ctss -> ['sl', ctss] }
                .mix(BAM_TO_CTSS.out.nosl_ctss.map { _meta, ctss -> ['nosl', ctss] })
                .groupTuple()
                .map { kind, ctss -> [[id: "all_samples.${kind}", sl_class: kind], ctss] }
        )
        ch_broad_ctss = POOL_CTSS_BY_SL.out.ctss
    }
    def ch_sets
    if (cluster_method == 'cagefightr') {
        def ch_class_ctss = sl_separate
            ? BAM_TO_CTSS.out.sl_ctss.map { _meta, ctss -> ['sl', ctss] }
                .mix(BAM_TO_CTSS.out.nosl_ctss.map { _meta, ctss -> ['nosl', ctss] })
            : BAM_TO_CTSS.out.ctss.map { _meta, ctss -> ['all', ctss] }
        // TPM within each class: SL CTSS relative to the sample's SL READ1s, non-SL CTSS to its non-SL READ1s, so that
        // a sample's share of trans-spliced reads does not change the threshold for its TSSs (and vice versa)
        CAGEFIGHTR_CLUSTERS(
            ch_class_ctss
                .groupTuple(sort: { a, b -> a.name <=> b.name })
                .map { kind, ctss -> [[id: "all_samples.${kind}", sl_class: kind], ctss] },
            ch_fasta_fai.map { meta, _fasta, fai -> [meta, fai] },
            []
        )
        def ch_cfr = CAGEFIGHTR_CLUSTERS.out.bed
        COMBINE_TSS_CLUSTERS(
            sl_separate
                ? ch_cfr.filter { meta, _bed -> meta.sl_class == 'nosl' }
                    .combine(ch_cfr.filter { meta, _bed -> meta.sl_class == 'sl' })
                    .map { _meta, bed, _sl_meta, sl_bed -> [[id: 'all_samples', cluster_set: 'broad'], bed, sl_bed] }
                : ch_cfr.map { _meta, bed -> [[id: 'all_samples', cluster_set: 'broad'], bed, []] }
        )
        ch_sets = COMBINE_TSS_CLUSTERS.out.bed
            .join(COMBINE_TSS_CLUSTERS.out.tsv)
            .map { _meta, bed, tsv -> ['broad', bed, tsv] }
    } else {
        PARACLU(ch_broad_ctss, min_cluster)
        def ch_paraclu = PARACLU.out.bed.join(ch_broad_ctss, failOnDuplicate: true, failOnMismatch: true)
        def ch_annotate = sl_separate
            ? ch_paraclu.filter { meta, _bed, _ctss -> meta.sl_class == 'nosl' }
                .combine(ch_paraclu.filter { meta, _bed, _ctss -> meta.sl_class == 'sl' })
                .map { _meta, bed, ctss, _sl_meta, sl_bed, sl_ctss -> [[id: 'all_samples', cluster_set: 'broad'], bed, ctss, sl_bed, sl_ctss] }
            : ch_paraclu.map { meta, bed, ctss -> [meta + [cluster_set: 'broad'], bed, ctss, [], []] }
        ANNOTATE_TSS_CLUSTERS(ch_annotate)
        ch_sets = ANNOTATE_TSS_CLUSTERS.out.bed
            .join(ANNOTATE_TSS_CLUSTERS.out.tsv)
            .map { _meta, bed, tsv -> ['broad', bed, tsv] }
    }

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
        .mix(CTSS_MOTIFS.out.initiator.map { _meta, tsv -> tsv }.collectFile(name: 'tss_initiator_mqc.tsv', keepHeader: true, skip: 10, sort: true))
        .mix(CTSS_MOTIFS.out.acceptor.map { _meta, tsv -> tsv }.collectFile(name: 'sl_acceptor_mqc.tsv', keepHeader: true, skip: 10, sort: true))
        .mix(CTSS_MOTIFS.out.stats.map { _meta, tsv -> tsv }.collectFile(name: 'ctss_motif_stats_mqc.tsv', keepHeader: true, skip: 5, sort: true))
        .mix(BAM_TO_CTSS.out.sl_leaders.map { _meta, tsv -> tsv }.collectFile(name: 'sl_leaders_mqc.tsv', keepHeader: true, skip: 10, sort: true))

    emit:
    ctss_bigwig       = UCSC_BEDGRAPHTOBIGWIG.out.bigwig                  // channel: [ val(meta), path(bigWig) ]
    ctss              = BAM_TO_CTSS.out.ctss                            // channel: [ val(meta), path(ctss.bed) ]
    softclip          = BAM_TO_CTSS.out.softclip                        // channel: [ val(meta), path(tsv) ]
    ctss_stats        = BAM_TO_CTSS.out.stats                           // channel: [ val(meta), path(ctss_stats_mqc.tsv) ]
    clusters_tsv      = ch_sets.map { set, _bed, tsv -> [set, tsv] }    // channel: [ val(cluster_set), path(tss_clusters.tsv) ]
    tss_bam           = ASSIGN_PAIRS_TO_TSS.out.bam                     // channel: [ val(meta), path(bam), path(bai) ]; meta.cluster_set
    counts            = MERGE_TSS_COUNTS.out.tsv                        // channel: [ val(meta), path(tsv) ]
    cagescan_clusters = MERGE_CAGESCAN_CLUSTERS.out.bed                 // channel: [ val(meta), path(bed12) ]
    multiqc_files     = ch_multiqc_files                                // channel: path(*_mqc.tsv)
}
