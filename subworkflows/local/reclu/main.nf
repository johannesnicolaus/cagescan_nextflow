//
// "Sharp" TSS clusters: RECLU-style reproducible clusters (Ohmiya et al. 2014)
//

include { PARACLU_RAW                                    } from '../../../modules/local/paraclu_raw/main'
include { RECLU_HIERARCHY                                } from '../../../modules/local/reclu_hierarchy/main'
include { RECLU_IDR                                      } from '../../../modules/local/reclu_idr/main'
include { RECLU_COMBINE                                  } from '../../../modules/local/reclu_combine/main'
include { ANNOTATE_TSS_CLUSTERS as ANNOTATE_SHARP_CLUSTERS } from '../../../modules/local/annotate_tss_clusters/main'

workflow RECLU {

    take:
    ch_ctss        // channel: [ val(meta), path(ctss.bed) ] per sample; meta.group = replicate group
    ch_pooled_ctss // channel: [ val(meta), path(ctss.bed) ] pooled across samples
    paraclu_min    // integer: paraclu minValue for the raw cluster hierarchy

    main:

    //
    // Full paraclu hierarchy per replicate, filtered by TPM per base and scored by
    // hierarchical stability
    //
    PARACLU_RAW(ch_ctss, paraclu_min)
    RECLU_HIERARCHY(PARACLU_RAW.out.clusters.join(ch_ctss, failOnDuplicate: true, failOnMismatch: true))

    //
    // Every pair of replicates within a group: reciprocal-overlap pairing + IDR
    //
    def ch_pairs = RECLU_HIERARCHY.out.bed
        .map { meta, bed -> [meta.group, [meta.id, bed]] }
        .groupTuple()
        .flatMap { group, reps ->
            def sorted = reps.sort { rep -> rep[0] }
            def n = sorted.size()
            (0..<n).collectMany { i ->
                ((i + 1)..<n).collect { j ->
                    [[id: "${group}.${sorted[i][0]}_vs_${sorted[j][0]}", group: group], sorted[i][1], sorted[j][1]]
                }
            }
        }
    ch_pairs.count().subscribe { n ->
        if (n == 0) {
            log.warn "No replicate group has 2 or more samples: sharp (RECLU) TSS clusters need replicates and are skipped. Add a 'group' column to the samplesheet."
        }
    }
    RECLU_IDR(ch_pairs)

    //
    // Unite reproducible clusters across pairs and groups, keep the innermost
    //
    RECLU_COMBINE(
        RECLU_IDR.out.tsv
            .map { _meta, tsv -> tsv }
            .collect()
            .map { tsvs -> [[id: 'all_samples', cluster_set: 'sharp'], tsvs] }
    )
    ANNOTATE_SHARP_CLUSTERS(
        RECLU_COMBINE.out.bed.combine(ch_pooled_ctss).map { meta, bed, _ctss_meta, ctss -> [meta, bed, ctss] }
    )

    emit:
    clusters = ANNOTATE_SHARP_CLUSTERS.out.bed.join(ANNOTATE_SHARP_CLUSTERS.out.tsv) // channel: [ val(meta), path(bed), path(tsv) ]
    stats    = RECLU_COMBINE.out.stats                                                // channel: [ val(meta), path(mqc.tsv) ]
}
