# GeneStructureTools

*R package for manipulation and analysis of spliced gene structures*
![](https://github.com/betsig/GeneStructureTools/blob/master/HexLogoGeneStructureTools.png)

## Installation
GeneStructureTools can be installed:
```
library(devtools)
install_github("betsig/GeneStructureTools")
```
After installation, the package can be loaded into R.
```
library(GeneStructureTools)
```

## Introduction

GeneStructureTools is an in-development package for the manipulation and analysis of transcribed gene structures.

We have provided functions for importing Whippet alternative splicing data, and the analysis of these splicing events. 
Splicing events can also be defined manually if you are using a different splicing analysis tool to Whippet.
For specific events - currently including exon skipping, intron retention, alternative splice site usage and alternative first/last exons - transcripts can be made in silico which use the two splicing modes - i.e. transcripts containing and transcripts skipping an exon. 
These transcripts do not have to be pre-annotated, and thus all potential isoforms can be compared for an event.

Current comparisons of transcripts include annotating and analysing ORF and UTR features (length, locations, difference/similarity between transcripts), and predicting nonsense mediated decay (NMD) potential.

We also have functions for re-annotation of .GTF features, such as annotating UTRs as 3' or 5', and assigning a broader biotype for genes and transcripts so more informative analysis can be performed between these classes. 

## Usage

Please check the [vignette](https://rawgit.com/betsig/GeneStructureTools/master/Vignette.html) for usage details.

## Development

GeneStructureTools is still undergoing development and documentation.
If you would like to contribute, please submit a pull request!
If you have any suggestions for features, please contact Beth Signal (b.signal@garvan.org.au).




