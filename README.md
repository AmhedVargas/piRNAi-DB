# piRNAi-DB: A collection of *C.elegans* piRNA interference (piRNAi) targets
Piwi-interacting RNAs, or simply piRNAs, are a class of small RNAs (21 nucleotides long) that are active in the germline and aid to the regulation of gene expression. One of the mechanisms by which they can perfom this control is by targeting cytosolic mRNA transcripts. While it remains unclear which are the precise rules governing this targeting [[1](https://science.sciencemag.org/content/359/6375/587),[2](https://www.sciencedirect.com/science/article/pii/S009286741830117X)], it's quite clear that the targeting occurs by means Watson-Crick pairing and, therefore, we can assume that by altering the sequence of the piRNAs we can change their targets; indeed this have been accomplished in *C. elegans* trhough a novel method denominated piRNA interfence (piRNAi).

This repository contains custom scripts used to find 20-mer sub-sequences of *C. elegans* and *C. briggsae* transcripts that can be used as targets for piRNAi. These sequences follow the guidelines required for piRNAi targeting, those being:

- **Unique** sequences in the genome, that is sequences that appears once in the genome. Furthermore, these sequences can be subdivided into categories that define their level of **uniqueness**, or as we defined, the maximum number of base substitutions we can perfom to the sequences without targeting other part of the genome. 
- No **thymine** between the first 3 bp of the sequence. This prevents shifts to the transcription frame.
- We recommend the use no extreme levels of GC content; altough our pipeline does not filer these cases.

## Software and hardware requirements
- Unix-like enviroment (Linux preffered but the pipelines can be run in macOS too).
- ~ 100 GB of space in hard drive 
- perl
- bedtools
- awk
- gzip
- blastn-suite
- bwa

## Folder structure


## Shinny implementation

## Note on finding the uniqueness of each sequence
