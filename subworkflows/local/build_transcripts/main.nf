//
// Build isoforms from TSS-anchored read pairs, anchor them to TSS clusters,
// and write GTF / GFF3 / BED12 / FASTA plus per-sample quantification
//

include { STRINGTIE_STRINGTIE as STRINGTIE_ASSEMBLE    } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE                              } from '../../../modules/nf-core/stringtie/merge/main'
include { ANCHOR_TRANSCRIPTS                           } from '../../../modules/local/anchor_transcripts/main'
include { GFFREAD as GFFREAD_GFF3                      } from '../../../modules/nf-core/gffread/main'
include { GFFREAD as GFFREAD_FASTA                     } from '../../../modules/nf-core/gffread/main'
include { STRINGTIE_STRINGTIE as STRINGTIE_QUANT       } from '../../../modules/nf-core/stringtie/stringtie/main'
include { MERGE_COUNTS as MERGE_TRANSCRIPT_COUNTS      } from '../../../modules/local/merge_counts/main'
include { GFFCOMPARE                                   } from '../../../modules/nf-core/gffcompare/main'

workflow BUILD_TRANSCRIPTS {

    take:
    ch_tss_bam      // channel: [ val(meta), path(bam), path(bai) ] TSS-anchored pairs; meta.cluster_set
    ch_clusters_tsv // channel: [ val(cluster_set), path(tss_clusters.tsv) ]
    ch_fasta        // channel: [ val(meta), path(fasta) ]
    ch_fai          // channel: [ val(meta), path(fai) ]
    ch_gtf          // channel: [ val(meta), path(gtf) ] or [ val(meta), [] ]
    use_guide       // boolean: give the reference GTF to StringTie as a guide
    has_gtf         // boolean: a reference GTF was provided (run gffcompare)

    main:

    def ch_guide = use_guide ? ch_gtf.map { _meta, gtf -> gtf } : channel.value([])

    //
    // Per-sample isoform assembly from TSS-anchored pairs
    //
    STRINGTIE_ASSEMBLE(
        ch_tss_bam.map { meta, bam, _bai -> [meta, bam, []] },
        '',
        ch_guide
    )

    //
    // Consensus transcript set across samples, per cluster set
    //
    STRINGTIE_MERGE(
        STRINGTIE_ASSEMBLE.out.transcript_gtf
            .map { meta, gtf -> [meta.cluster_set, gtf] }
            .groupTuple()
            .map { set, gtfs -> [[id: 'merged', cluster_set: set], gtfs] },
        use_guide ? ch_gtf : channel.value([[:], []])
    )

    //
    // Anchor per-sample and merged transcripts to the TSS clusters of their set
    //
    def ch_anchor = STRINGTIE_ASSEMBLE.out.transcript_gtf
        .mix(STRINGTIE_MERGE.out.merged_gtf)
        .map { meta, gtf -> [meta.cluster_set, meta, gtf] }
        .combine(ch_clusters_tsv, by: 0)
        .multiMap { _set, meta, gtf, tsv ->
            gtf:      [meta, gtf]
            clusters: tsv
        }
    ANCHOR_TRANSCRIPTS(ch_anchor.gtf, ch_anchor.clusters)
    def ch_anchored_gtf = ANCHOR_TRANSCRIPTS.out.gtf
    def ch_merged_gtf   = ch_anchored_gtf.filter { meta, _gtf -> meta.id == 'merged' }

    //
    // GFF3 and transcript sequences
    //
    GFFREAD_GFF3(ch_anchored_gtf, [])
    GFFREAD_FASTA(ch_anchored_gtf, ch_fasta.map { _meta, fasta -> fasta })

    //
    // Quantify the merged transcript set of each cluster set in every sample
    //
    def ch_quant = ch_tss_bam
        .map { meta, bam, _bai -> [meta.cluster_set, meta, bam] }
        .combine(ch_merged_gtf.map { meta, gtf -> [meta.cluster_set, gtf] }, by: 0)
        .multiMap { _set, meta, bam, gtf ->
            bam: [meta, bam, []]
            gtf: gtf
        }
    STRINGTIE_QUANT(ch_quant.bam, 'expression-estimation', ch_quant.gtf)
    MERGE_TRANSCRIPT_COUNTS(
        STRINGTIE_QUANT.out.transcript_gtf
            .map { meta, gtf -> [meta.cluster_set, gtf] }
            .groupTuple()
            .map { set, gtfs -> [[id: 'all_samples', cluster_set: set], gtfs] },
        [],
        'transcripts'
    )

    //
    // Compare with the reference annotation (known vs novel)
    //
    def ch_gffcompare_stats = channel.empty()
    if (has_gtf) {
        GFFCOMPARE(
            ch_merged_gtf,
            ch_fasta.combine(ch_fai).map { meta, fasta, _fai_meta, fai -> [meta, fasta, fai] }.first(),
            ch_gtf
        )
        ch_gffcompare_stats = GFFCOMPARE.out.stats
    }

    def ch_multiqc_files = ANCHOR_TRANSCRIPTS.out.stats
        .map { _meta, stats -> stats }
        .collectFile(name: 'anchor_stats_mqc.tsv', keepHeader: true, skip: 5, sort: true)

    emit:
    gtf             = ch_anchored_gtf                    // channel: [ val(meta), path(gtf) ]
    gff3            = GFFREAD_GFF3.out.gffread_gff       // channel: [ val(meta), path(gff3) ]
    fasta           = GFFREAD_FASTA.out.gffread_fasta    // channel: [ val(meta), path(fasta) ]
    bed             = ANCHOR_TRANSCRIPTS.out.bed         // channel: [ val(meta), path(bed12) ]
    counts          = MERGE_TRANSCRIPT_COUNTS.out.tsv    // channel: [ val(meta), path(tsv) ]
    gffcompare      = ch_gffcompare_stats                // channel: [ val(meta), path(stats) ]
    multiqc_files   = ch_multiqc_files                   // channel: path(*_mqc.tsv)
}
