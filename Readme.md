# An integrative R package for analysing alternative splicing using RNAseq

## Authors
Estefania Mancini, Andrés Rabinovich, Javier Iserte, Marcelo Yanovsky, Ariel Chernomoretz

## Introduction

Alternative splicing (AS) is a common mechanism of post-transcriptional gene 
regulation in eukaryotic organisms that expands the functional and regulatory 
diversity of a single gene by generating multiple mRNA isoforms that encode 
structurally and functionally distinct proteins. The development of novel 
high-throughput sequencing methods for RNA (RNAseq) has provided a powerful 
means of studying AS under multiple conditions and in a genome-wide manner.
However, using RNAseq to study changes in AS under different experimental 
conditions is not trivial. 
In this vignette, we describe how to use ASpli, an integrative and user-friendly
R package that facilitates the analysis of changes in both annotated and novel 
AS events. This package combines statistical information from exon, intron, and 
splice junction differential usage (p-value, FDR), with information from splice 
junction reads to calculate differences in the percentage of exon inclusion 
Delta-PSI and intron retention Delta-PIR). The proposed methodology 
reliably reflect the magnitude of changes in the relative abundance of different 
annotated and novel AS events. This method can be used to analyze both simple 
and complex experimental designs involving multiple experimental conditions.

## Getting started

### Installation
library(devtools)
install_git("https://gitlab.com/chernolab/aspli2.git")
library(ASpli2)


Changes

