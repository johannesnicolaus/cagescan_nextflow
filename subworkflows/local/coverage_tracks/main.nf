//
// Stranded read coverage (bigWig, per million mapped read pairs) for genome browsers
//

include { SAMTOOLS_VIEW                                  } from '../../../modules/nf-core/samtools/view/main'
include { BEDTOOLS_GENOMECOV                             } from '../../../modules/nf-core/bedtools/genomecov/main'
include { UCSC_BEDGRAPHTOBIGWIG as COVERAGE_BIGWIG       } from '../../../modules/nf-core/ucsc/bedgraphtobigwig/main'

workflow COVERAGE_TRACKS {

    take:
    ch_bam_bai    // channel: [ val(meta), path(bam), path(bai) ]
    ch_ctss_stats // channel: [ val(meta), path(ctss_stats_mqc.tsv) ] (READ1s used = mapped read pairs)
    ch_sizes      // channel: [ val(meta), path(chrom.sizes) ]

    main:

    // Primary, uniquely mapped alignments only (same MAPQ filter as the CTSS)
    SAMTOOLS_VIEW(ch_bam_bai, [[:], [], []], [[:], []], [[:], []], '')

    def ch_scale = ch_ctss_stats.map { meta, stats ->
        def pairs = stats.readLines().last().split('\t')[1] as double
        [meta, pairs > 0 ? 1e6 / pairs : 1]
    }

    // CAGE READ1 is sense to the RNA: with -du the mate takes READ1's strand
    def ch_genomecov = SAMTOOLS_VIEW.out.bam
        .join(ch_scale, failOnDuplicate: true, failOnMismatch: true)
        .combine(['plus', 'minus'])
        .map { meta, bam, scale, strand -> [meta + [strand: strand], bam, scale] }
    BEDTOOLS_GENOMECOV(ch_genomecov, ch_sizes.map { _meta, sizes -> sizes }, 'bedgraph', true)

    COVERAGE_BIGWIG(
        BEDTOOLS_GENOMECOV.out.genomecov.filter { _meta, bedgraph -> bedgraph.size() > 0 },
        ch_sizes.map { _meta, sizes -> sizes }
    )

    emit:
    bigwig = COVERAGE_BIGWIG.out.bigwig // channel: [ val(meta), path(bigWig) ]; meta.strand
}
