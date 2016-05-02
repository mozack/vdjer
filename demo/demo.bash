#!/bin/bash
#
# Demo script for V'DJer

# Before running this script:
# Run make to generate the vdjer executable.
# Download the vdjer references and untar, updating the VDJER_REF_DIR variable below.

VDJER_REF_DIR=<path/to/vdjer/referenes>/igh

VDJER=../vdjer

# Run V'DJer on the input BAM file for the IgH chain.
# read length = 50
# median insert size = 175
# num threads = 4
# 
# Assembled contigs appear in vdj_contigs.fa
# Reads mapped to contigs appear in vdjer.sam
$VDJER --in star.sort.bam --rl 50 --t 4 --ins 175 --chain IGH --ref-dir $VDJER_REF_DIR > vdjer.sam 2> vdjer.log