#!/usr/bin/env Rscript
args <- commandArgs(TRUE)
if (length(args) != 3L) stop("usage: overlap_regions.R LEFT.csv RIGHT.csv OUTPUT.csv")
suppressPackageStartupMessages(library(GenomicRanges))
read_ranges <- function(path) {
  d <- read.csv(path, stringsAsFactors = FALSE)
  if (!all(c("region_id", "chrom", "start", "end") %in% names(d))) stop("区域表缺少坐标或region_id")
  if (anyNA(d$start) || anyNA(d$end) || any(d$start < 1 | d$end < d$start)) stop("需要1-based闭区间")
  GRanges(seqnames = sub("^chr", "", d$chrom), ranges = IRanges(d$start, d$end), region_id = d$region_id)
}
a <- read_ranges(args[1L]); b <- read_ranges(args[2L])
h <- findOverlaps(a, b, ignore.strand = TRUE)
i <- queryHits(h); j <- subjectHits(h)
out <- data.frame(left_region = a$region_id[i], right_region = b$region_id[j],
                  overlap_bp = pmin(end(a)[i], end(b)[j]) - pmax(start(a)[i], start(b)[j]) + 1L)
write.csv(out, args[3L], row.names = FALSE)
