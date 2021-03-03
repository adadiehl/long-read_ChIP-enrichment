#!/bin/bash

# Check bam alignments for enrichment of CTCF reads after first filtering out duplicate and low-quality reads.

BAM_IN=$1
PEAK_F=$2

WD=$(pwd)
BAM_ROOT=$(basename $BAM_IN .sorted.bam)

export PATH=$WD:$WD/scripts:$PATH

# Mark duplicates with Picard tools
echo "Checking for duplicate reads"
module load picard/2.18.0
java -jar $EBROOTPICARD/picard.jar MarkDuplicates I=$BAM_IN O=$BAM_ROOT.sorted.md.bam ASSUME_SORTED=true METRICS_FILE=.picard.metrics VALIDATION_STRINGENCY=LENIENT TMP_DIR=. 2> /dev/null

# Filter the marked bam
echo "Filtering on quality, mapping, and duplication status"
module load SAMtools/1.5
samtools view -@ 12 -b -h -F 4 -F 256 -F 1024 -F 2048 -q 30 $BAM_ROOT.sorted.md.bam > $BAM_ROOT.sorted.filtered.bam

# Convert to bed format
echo "Converting to BED format"
#module load SAMtools/1.5
bamToBed -i $BAM_ROOT.sorted.filtered.bam > $BAM_ROOT.sorted.filtered.bed

# Intersect mapped and filtered reads with given features.
echo "Checking for intersection with reference features"
module load BEDTools/2.26.0
printf "n_features\tn_reads\tn_overlapping_reads\n" > $BAM_ROOT.intersection.txt
printf "%d\t%d\t%d\n" $(zcat $PEAK_F | wc -l) $(cat $BAM_ROOT.sorted.filtered.bed | wc -l) $(bedtools intersect -a $BAM_ROOT.sorted.filtered.bed -b $PEAK_F | wc -l) >> $BAM_ROOT.intersection.txt

# Calculate enrichment p-values in R.
echo "Calculating enrichment p-value based on random intersections."
module load R
./scripts/calc_enrichments.R $BAM_ROOT.intersection.txt $BAM_ROOT.sorted.filtered.bed $PEAK_F $BAM_ROOT.enrichment.txt

echo "Done!"
