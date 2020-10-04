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
C_elegans/C_briggsae-pipeline
- Stand-alone shell script that produces piRNA tracks
- Custom made perl scripts required to identify unique k-mer sequences.
Example_results
- tar.gz 5MM track
**Please note that due to space limitations we just uploaded one of the resultant outputs of the pipeline. Running the main shell scripts should recreate this and all other files used for piRNAi target identification** 

## Shinny implementation
Each resultan bed file was treated to be incorporated into the [piRNAi app](https://wormbuilder.dev/piRNAi/). The full description of how the files were prepared and the shinny app implemented can be found [here](https://github.com/AmhedVargas/piRNAi_v2). 

## Note on finding the uniqueness of each sequence
Finding words a defined length that appear once in a text is a trivial matter in computation; and so it's the case of finding sequences of n base pairs (k-mer) appearing once in a genome. However, the inclusion of mis-matches makes the number of computations required to solve this problem grow exponentially [[3]](https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7244195) rendering almost impossible to use a deterministic solution. For this reason, our pipeline used back and forth euristic and deterministic methods to come up with the list of uniquely targeting piRNAi sequences (see below).

## Rational progression of our piRNAi targetting algorithm
1. Download latest genomic files (sequence and gene annotations) of *C.elegans*/*C.briggsae* from [wormbase ftp server](ftp://ftp.wormbase.org/pub/wormbase). **Note:** When the pipeline was originally conceived the version used was WS270
2. Extract solely the CDS of protein coding genes in fasta format.
3. For each CDS, create the list of all 20 bp long sequences (overlapping 20-mers) targetting the reverse complement of the CDS.
4. Determine which 20-mers are unique in the genome (both strands) allowing 0 mismatches (0MM) via a hash method.
5. Convert the list of unique 20-mers back to fasta format while, at the same time, filtering sequences with any T between the first 3 base pairs (to exclude errouneous transcription).
6. Map the 20-mer sequences back to the genome using bwa; this mapping procedure not only allows us to recover the coordinates of 20-mers but also helps in the removal of duplicated or non-unique sequences.
7. For the curated list of unique 20-mers, we then start a series of back-and-forth blasting and mapping procedures to filter sequences with n mismatches:
- Assess uniqueness of subquences of length (20 - n) in the genome via a hash method. In other words, check for example that both 19-mers found within a 20 mer are unique (i.e. unique 20-mer that allows 1 mis-match either at the end or beginning of its sequnce). **This step is not mandatory *per se*, however, it greatly helps the following steps as it reduces the number of sequences that will be tested**
- Find any relevant and not relevant match for each 20-mer in the genome by using blast (e-value cut off equal to 1)
- Extract sequences that had solely one alignement satisfying the condition: size of alignment >= length of fragment - mismatches allowed (20 - n)
- Map them back to the genome to get coordinates and repeat the procedure once again
8. Finally, the resultant bed or fasta file can be filtered for GC content if needed.

