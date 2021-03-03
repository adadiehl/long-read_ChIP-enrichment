#!/bin/bash

# Runs a simple intersection with bedtools intersect and returns the number of unique features intersecting.
module load BEDTools/2.26.0
bedtools intersect -a $1 -b $2 | wc -l
