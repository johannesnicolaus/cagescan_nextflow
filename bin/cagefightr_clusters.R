#!/usr/bin/env Rscript
# Unidirectional TSS clusters with CAGEfightR (Thodberg et al. 2019), following its standard workflow:
#   CTSSs -> TPM -> pooled TPM -> CTSSs supported in >= --ctss-min-samples samples
#   -> clusterUnidirectionally (slice-reduce, pooledCutoff 0, mergeDist --merge-dist)
#   -> TCs with >= --tc-min-tpm TPM in >= --tc-min-samples samples
# The peak (thick) is the CTSS with the highest pooled TPM.
#
# usage: cagefightr_clusters.R --prefix P --fai genome.fa.fai [options] S1.ctss.bed S2.ctss.bed ...
# CTSS bed: chrom, pos0, pos0+1, ., count, strand.  Output: P.cagefightr.bed
#   chrom, start0, end, strand, peak0, total count, pooled TPM
suppressPackageStartupMessages({
    library(CAGEfightR)
    library(GenomicRanges)
    library(SummarizedExperiment)
    library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
opt <- list(prefix = NULL, fai = NULL, ctss_min_samples = 2, merge_dist = 20, tc_min_tpm = 1, tc_min_samples = 2)
files <- character()
i <- 1
while (i <= length(args)) {
    a <- args[i]
    if (startsWith(a, "--")) {
        opt[[gsub("-", "_", substring(a, 3))]] <- args[i + 1]
        i <- i + 2
    } else {
        files <- c(files, a)
        i <- i + 1
    }
}
for (k in c("ctss_min_samples", "merge_dist", "tc_min_samples")) opt[[k]] <- as.integer(opt[[k]])
opt$tc_min_tpm <- as.numeric(opt$tc_min_tpm)
stopifnot(!is.null(opt$prefix), !is.null(opt$fai), length(files) > 0)

samples <- sub("\\..*$", "", basename(files))
fai <- read.delim(opt$fai, header = FALSE, colClasses = c("character", "numeric", rep("NULL", 3)))
seqinfo <- Seqinfo(fai[[1]], fai[[2]])

# CTSS x sample count matrix
tabs <- lapply(seq_along(files), function(j) {
    x <- read.delim(files[j], header = FALSE, colClasses = c("character", "integer", "NULL", "NULL", "numeric", "character"))
    names(x) <- c("chrom", "pos", "count", "strand")
    x$sample <- j
    x
})
all <- do.call(rbind, tabs)
out_file <- paste0(opt$prefix, ".cagefightr.bed")
if (nrow(all) == 0) {
    file.create(out_file)
    quit(save = "no")
}
key <- paste(all$chrom, all$pos, all$strand, sep = "\t")
ukey <- unique(key)
parts <- do.call(rbind, strsplit(ukey, "\t", fixed = TRUE))
gr <- GRanges(parts[, 1], IRanges(as.integer(parts[, 2]) + 1L, width = 1L), strand = parts[, 3], seqinfo = seqinfo)
counts <- sparseMatrix(i = match(key, ukey), j = all$sample, x = all$count, dims = c(length(ukey), length(files)),
                       dimnames = list(NULL, samples))
o <- order(gr)
gr <- gr[o]
counts <- counts[o, , drop = FALSE]
CTSSs <- SummarizedExperiment(assays = SimpleList(counts = counts), rowRanges = GPos(gr),
                              colData = DataFrame(Name = samples, row.names = samples))

n <- length(samples)
CTSSs <- calcTPM(CTSSs, inputAssay = "counts", outputAssay = "TPM", outputColumn = "totalTags")
CTSSs <- calcPooled(CTSSs, inputAssay = "TPM")
# subsetBySupport keeps features above `unexpressed` in more than `minSamples` samples
supported <- subsetBySupport(CTSSs, inputAssay = "counts", outputColumn = "support", unexpressed = 0,
                             minSamples = min(opt$ctss_min_samples, n) - 1)
# CAGEfightR's slice-reduce, done directly: used when clusterUnidirectionally fails on very small inputs
# (its cluster statistics break when a strand has a single cluster)
slice_reduce <- function(x, merge_dist) {
    gr <- granges(rowRanges(x))
    score <- rowRanges(x)$score
    keep <- score > 0
    gr <- gr[keep]
    score <- score[keep]
    red <- reduce(gr, min.gapwidth = merge_dist, with.revmap = TRUE)
    pos <- extractList(start(gr), red$revmap)
    sc <- extractList(score, red$revmap)
    peak <- mapply(function(p, s) p[which.max(s)], as.list(pos), as.list(sc))
    red$revmap <- NULL
    red$score <- sum(sc)
    red$thick <- IRanges(peak, width = 1L)
    names(red) <- paste0(seqnames(red), ":", start(red), "-", end(red), ";", strand(red))
    red
}
TCs <- tryCatch(clusterUnidirectionally(supported, pooledCutoff = 0, mergeDist = opt$merge_dist),
                error = function(e) {
                    message("clusterUnidirectionally failed (", conditionMessage(e), "); using the same slice-reduce directly")
                    slice_reduce(supported, opt$merge_dist)
                })
TCexp <- quantifyClusters(CTSSs, clusters = TCs, inputAssay = "counts")
TCexp <- calcTPM(TCexp, inputAssay = "counts", outputAssay = "TPM", totalTags = "totalTags")
TCexp <- subsetBySupport(TCexp, inputAssay = "TPM", outputColumn = "support", unexpressed = opt$tc_min_tpm,
                         minSamples = min(opt$tc_min_samples, n) - 1)

r <- rowRanges(TCexp)
res <- data.frame(chrom = as.character(seqnames(r)), start = start(r) - 1L, end = end(r), strand = as.character(strand(r)),
                  peak = start(r$thick) - 1L, total = as.integer(round(rowSums(assay(TCexp, "counts")))),
                  pooled = signif(r$score, 6))
res <- res[order(res$chrom, res$start, res$strand), ]
write.table(res, out_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
message(sprintf("%s: %d CTSSs, %d supported, %d TCs, %d after expression filter", opt$prefix, length(CTSSs), length(supported), length(TCs), nrow(res)))
