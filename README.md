# PHGv2Tools
Package to downstream analysis of pangenomes databases, working with Practical Haplotype Graph and its h.VCF files

## Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [References](#references)
- [License](#license)

## Introduction
Repository writen in python to perform downstream analisys of [Practical Haplotype Graph](https://phg.maizegenetics.net/). It is prepared to either work in command line to analyse the pangenome database, but also [notebooks](#notebooks) are available to be used and modified if needed by any user. 
It mainly works with [h.VCF](https://phg.maizegenetics.net/hvcf_specifications/) files, a modified version of [2.4 VCF](http://samtools.github.io/hts-specs/VCFv4.2.pdf).

It basically works with two outputs of PHG:

#### Pangenome haplotypes VCFs
With the function [Create VCF files](#https://phg.maizegenetics.net/build_and_load/#create-vcf-files) of PHG a graph formed by nodes and ranges is performed for each of the genomes. It is usefull to  [merge the VCFs](https://github.com/maize-genetics/phg_v2/blob/main/src/main/kotlin/net/maizegenetics/phgv2/cli/MergeHvcfs.kt) to do analysis of the graph as a whole. 
Functions supported now:
- [Core, accesion and unique ranges analisys](#Core-accesion-and-unique-ranges-analisys)
- [Pangenome ranges evolution](#Pangenome-ranges-evolution)
- [Plot pangenome regions/chromosomes]
- [Check haplotypes for a region]

#### Imputed haplotypes VCFs
[Imputation](https://phg.maizegenetics.net/imputation/) is a powerfull tool from PHG to achieve complete genomes even from low density sequence  or variant data. What is obtained after alligning kmers and getting the graph of the haplotype is a h.VCF file, which some data can be mined from it:
- [Check identity against pangenome]
- [Plot imputed genome]


## Installation
To use this package a conda environment is used. It is a modified version of the original phg one, adding the pygenometracks package for python plotting. For everything further needed, it will be updated.
(Here i should add the link of reference of both githubs repos, including the installation guide at phg.

#### Download & install updated package
```
git clone https://github.com/jsarriaa/PHGv2Tools.git
```
Load it
```
cd PHGv2Tools
```
Install it
```
pip install .
```
#### Install conda
Prepare
```
chmod +x CondaSetup.sh
```
Install phgtools conda
```
./CondaSetup.sh
```
Activate phgtools conda
```
conda activate phgtools
```
#### Install PHG
Ensure to have installated PHGv2. Follow the steps [here](https://phg.maizegenetics.net/installation/)

## Usage
- #### Core, accesion and unique ranges analisys
- #### Pangenome ranges evolution






## References
The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation  https://doi.org/10.1093/bioinformatics/btac410
