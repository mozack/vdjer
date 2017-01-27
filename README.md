# V'DJer - B Cell Receptor repertoire reconstruction from short read mRNA-Seq data

V'DJer uses customized read extraction, assembly and V(D)J rearrangement detection 
and filtering to produce contigs representing the most abundant portions of the BCR
repetoire.  V'DJer operates upon bulk short read mRNA-Seq data allowing sequence analysis
typically done using specialized long read sequencing.  The software can be run on
BCR heavy (IgH) and light (Igk/IgL) chains.

V'DJer accepts a sorted and indexed BAM file mapped to hg38 as input.  We currently
use STAR for alignment.  At the end of a run assembled contigs will appear in
```vdj_contigs.fa``` and reads aligned to those contigs are written to ```stdout```

The V'DJer output is suitable for use by downstream quantification tools such as [RSEM](http://deweylab.github.io/RSEM/).

## Getting the code and building the software

Use a recent release: https://github.com/mozack/vdjer/releases

To compile, just cd into the vdjer directory and type ```make```.  An executable named
```vdjer``` will be created.

To date, V'DJer has been tested only on linux systems.

Pre-built human indices and references can be downloaded from here: https://github.com/mozack/vdjer/releases/download/v0.10_reference/vdjer_human_references.tar.gz

The archive must be untarred and decompressed for use by V'DJer

A reference directory is included for each chain type.

## Usage:

V'DJer operates on paired end reads.  Single end reads are not currently supported.
When mapping with STAR unmapped reads must be included in the BAM file (this is not 
STAR's default behavior).  Include the unmapped reads using the following param:

```--outSAMunmapped Within```

The following runs vdjer on the input star.sort.bam file with vdj_contigs.fa
generated in the current working directory and read alignments written to vdjer.sam:

```vdjer --in star.sort.bam --t 8 --ins 175 --chain IGH --ref-dir vdjer_human_references/igh > vdjer.sam 2> vdjer.log```

The above runs on the IgH chain with read length of 50, 8 threads and median insert length of 175.

## Important Parameters
parameter | value
------ | -------
--in | input BAM
--t | num threads
--ins | median insert size
--chain | one of: IGH, IGK or IGL
--ref-dir | chain specific reference directory

## Sensitive mode:

In cases where a sample has low BCR expression levels, one may opt to run V'DJer using more sensitive settings.
Running sensitive mode on samples with abundant BCR expression levels may result in extremely high computational costs.
Example command:

```vdjer --in star.sort.bam --t 8 --ins 175 --chain IGH --ref-dir vdjer_human_references/igh --k 25 
--mq 60 --mf 2 --rs 25 --ms 2 --mcs -5.5 > vdjer.sam 2> vdjer.log ```

Decreasing mq and mf from the defaults results in less aggressive graph pruning.
Reducing rs and ms results in less aggressive coverage based filtering
Decreasing mcs allows for more exhaustive graph traversal

The values used in this example match those used when running sensitive mode in the V'DJer paper. 

## Demo
See demo.bash and quant_demo.bash under the demo directory for an example of running V'DJer.

## Manuscript
https://doi.org/10.1093/bioinformatics/btw526
