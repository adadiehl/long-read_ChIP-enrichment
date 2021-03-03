#!/usr/bin/Rscript --vanilla

# Calculate enrichments for the given intersection data. Input is expected to be
# a three-column file with columns:
# 1) total reads in "reference" features
# 2) total reads in "query" features
# 3) number of reads in intersection.
#
# Tests performed will be for enrichment of "reference" features in the set of
# "query" features compared to overlaps in a random set of genomic regions.

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 4) {
   stop("Four arguments required: input_data.txt query_features.bed reference_features.bed outfile.txt.n", call.=FALSE)
}

require("rtracklayer")
require("regioneR")

# Load reads and get lengths for background feature generation.
reads = read.table(args[2], header=FALSE, stringsAsFactors=FALSE, sep="\t", col.names=c("chrom", "start", "end", "name", "score", "strand"))
reads$length = reads$end - reads$start

# Note that the command below allows overlapping regions, which seems to be a rational choice, since reads in the sequenced sample may overlap.
bg_regions = createRandomRegions(nregions=nrow(reads), length.mean = mean(reads$length), length.sd = sd(reads$length), genome = "hg38", non.overlapping = FALSE)

# Write to disk for outside comparison with bedtools...
write.table(bg_regions, col.names=FALSE, row.names=FALSE, sep="\t", quote=FALSE, file="random_seqs.bed")
cmdargs = paste("random_seqs.bed", args[3], sep=" ")
bg_int_count = as.numeric(system2("run_bedtools.sh", cmdargs, stdout=TRUE))

# Load the intersection data
int_data = read.table(args[1], header=TRUE, stringsAsFactors=FALSE, sep="\t")

# Run the enrichment test
fdata = t(matrix( c(as.numeric(int_data[3]), nrow(reads)-as.numeric(int_data[3]), bg_int_count, length(bg_regions)-bg_int_count), nrow=2, ncol=2))
f = fisher.test(fdata, alternative="g")

# Add results to table
int_data$bg_overlapping_features = bg_int_count
int_data$fisher_pval = f$p.value

# Write results
write.table(int_data, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE, file=args[4])
