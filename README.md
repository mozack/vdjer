# V'DJer - B Cell Receptor repertoire reconstruction from short read mRNA-Seq data

V'DJer uses customized read extraction, assembly and V(D)J rearrangement detection 
and filtering to produce contigs representing the most abundant portions of the BCR
repetoire.  V'DJer operates upon short read mRNA-Seq data allowing sequence analysis
typically done using specialized long read sequencing.  The software can be run on
BCR heavy (IgH) and light (Igk/IgL) chains.

V'DJer accepts a sorted and indexed BAM file mapped to hg38 as input.  We currently
use STAR for alignment.  At the end of a run assembled contigs will appear in
```vdj_contigs.fa``` and reads aligned to those contigs will appear in ```vdjician.sam```

The V'DJer output is suitable for use by downstream quantification tools such as RSEM.

## Getting the code and building the software

Use a recent release: https://github.com/mozack/vdjer/releases

To compile, just cd into the vdjer directory and type ```make```.  An executable named
```vdjer``` will be created.

Pre-built human indices and references can be downloaded from here: TODO

The archive must be untarred and decompressed for use by V'DJer

A reference directory is included for each chain type.

## Sample usage:

The following runs vdjer on the input star.sort.bam file with vdj_contigs.fa and vdjician.sam
generated in the current working directory:

```vdjer --in star.sort.bam --rl 50 --t 8 --ins 175 --chain IGH --ref-dir vdjer_human_references/igh```

The above runs on the IgH chain with read length of 50, 8 threads and median insert length of 175.

## Important Parameters
```
--in <input BAM>
--rl <read length>
--t <threads>
--ins <median insert size>
--chain <IGH|IGK|IGL>
--ref-dir <chain specific reference directory>
```

## Demo
See demo.bash and quant_demo.bash under the demo directory for an example of running V'DJer.