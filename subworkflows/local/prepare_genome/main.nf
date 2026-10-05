//
// Uncompress and index the reference genome
//

include { GUNZIP as GUNZIP_FASTA } from '../../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF   } from '../../../modules/nf-core/gunzip/main'
include { SAMTOOLS_FAIDX         } from '../../../modules/nf-core/samtools/faidx/main'
include { STAR_GENOMEGENERATE    } from '../../../modules/nf-core/star/genomegenerate/main'

workflow PREPARE_GENOME {

    take:
    fasta      // string: path to genome FASTA (may be gzipped)
    gtf        // string: path to annotation GTF (optional, may be gzipped)
    star_index // string: path to a pre-built STAR index directory (optional)

    main:

    def ch_fasta = channel.value([[id: 'genome'], file(fasta, checkIfExists: true)])
    if (fasta.endsWith('.gz')) {
        ch_fasta = GUNZIP_FASTA(ch_fasta).gunzip
    }

    def ch_gtf = channel.value([[id: 'annotation'], []])
    if (gtf) {
        ch_gtf = channel.value([[id: 'annotation'], file(gtf, checkIfExists: true)])
        if (gtf.endsWith('.gz')) {
            ch_gtf = GUNZIP_GTF(ch_gtf).gunzip
        }
    }

    SAMTOOLS_FAIDX(ch_fasta.map { meta, fa -> [meta, fa, []] }, true)

    def ch_star_index
    if (star_index) {
        ch_star_index = channel.value([[id: 'star'], file(star_index, checkIfExists: true)])
    } else {
        ch_star_index = STAR_GENOMEGENERATE(ch_fasta, ch_gtf).index
    }

    emit:
    fasta      = ch_fasta                             // channel: [ val(meta), path(fasta) ]
    fai        = SAMTOOLS_FAIDX.out.fai               // channel: [ val(meta), path(fai) ]
    sizes      = SAMTOOLS_FAIDX.out.sizes             // channel: [ val(meta), path(sizes) ]
    gtf        = ch_gtf                               // channel: [ val(meta), path(gtf) ] or [ val(meta), [] ]
    star_index = ch_star_index                        // channel: [ val(meta), path(index) ]
}
