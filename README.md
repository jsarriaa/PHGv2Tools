# PHGv2Tools
Package to downstream analysis of pangenomes databases, working with Practical Haplotype Graph and its h.VCF files

## Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [References](#references)

## Introduction
Repository writen in python to perform downstream analisys of [Practical Haplotype Graph](https://phg.maizegenetics.net/). It is prepared to either work in command line to analyse the pangenome database, but also [notebooks](#notebooks) are available to be used and modified if needed by any user. 
It mainly works with [h.VCF](https://phg.maizegenetics.net/hvcf_specifications/) files, a modified version of [2.4 VCF](http://samtools.github.io/hts-specs/VCFv4.2.pdf).

It basically works with two outputs of PHG:

#### Pangenome haplotypes VCFs
The pangenome graph is build from the reference genome. Ranges based in annotated genes are described, stablishing by starting and ending nodes, fisic coordinates at a chromosome. Alligning new genomes against it determines the presence/absence of each range at this genome. Then, it is posible to add to the pangenome, retaining also this new haplotype for the range, enriching the variability for it. With the function [Create VCF files](#https://phg.maizegenetics.net/build_and_load/#create-vcf-files) of PHG an haplotype file for each genome is obtained, as an h.VCF. It is usefull to  [merge the VCFs](https://github.com/maize-genetics/phg_v2/blob/main/src/main/kotlin/net/maizegenetics/phgv2/cli/MergeHvcfs.kt) to do analysis of the graph as a whole. 
Functions supported now:
- [Core, accesion and unique ranges analisys](#Core-accesion-and-unique-ranges-analisys)
- [Pangenome ranges evolution](#Pangenome-ranges-evolution)
- [Plot pangenome regions/chromosomes](#Plot-pangenome-regions/chromosomes)
- [Check haplotypes for a region](#Check-haplotypes-for-a-region)

#### Imputed haplotypes VCFs
[Imputation](https://phg.maizegenetics.net/imputation/) is a powerfull tool from PHG to achieve complete genomes even from low density sequence  or variant data. What is obtained after alligning kmers and getting the graph of the haplotype is a h.VCF file, which some data can be mined from it:
- [Check identity against pangenome](#Check-identity-against-pangenome)
- [Plot imputed genome](#Plot-imputed-genome)


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
Having a merged h.VCF file of the pangenome PHG database, it is possible to determine which ranges are found in all, only one or in certain genomes of the pangenome. As PHG ranges are based in genes annotation, it is possible to perform a parallelism between them. It is stapled for pangenome data analysis to think about core, accessory and unique genes:

![Pangenome genes example. [reference](https://www.bgi.com/us/plant-and-animal-pan-genome)](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/Figure%201.%20Example%20of%20a%20pan-genome.jpg)

[img reference](https://www.bgi.com/us/plant-and-animal-pan-genome)
Extrapolating this into the PHG ranges, it is usefull to have the information of which ranges are shared or not.

```core-range-detecter``` is the function which will give as output 2 png images:
  - bar plot of ranges distributed in how many genomes are found at pangenome.
  - sector diagram plotting the % of core, accessory and unique haplotypes.
  It takes as arguments:
```
  --help
  --pangenome-hvcf // -hvcf  <Merged of all hvcf of the pangenome database>
```
![image_1_hvcf_plot](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/merged_hvcfs_19092024.h.vcf_1.png)
![image_2_hvcf_plot](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/merged_hvcfs_19092024.h.vcf_2.png)

- #### Pangenome ranges evolution

A pangenome store all variability from species that a single reference genome can not. However, as long as increasing the ammount of varieties included in a pangenome the growing slope breaks exponentially, showing a collapsed top. To plot how does this ranges storage in each species pangenome is usefull in order to optimize the ammount of data stored.

``` range-pangenome-evolution``` is the function stacking one by one the haplotype files and plotting the number of ranges. Taking as arguments:

```
--hvcf-folder // -hf <Folder path to the h.vcf files of the pangenome database>  #Built with phg create-maf-vcf
--reference-file // -ref <Reference.fasta>  #Built with phg prepare-assemblies
--range-bedfile // <reference_ranges.bed>   #Built with phg create-ranges
```
![range_evolution](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/range_evolution.png)

  
- #### Plot pangenome regions/chromosomes
Function to plot a region of a chromosome (or the whole chr) of the pangenome.
If not region is specified, it will be plotted the entire chromosome. The function is ```plot-pangenome-chromosomes``` and its agruments are:
```
--hvcf-folder // -hvcf <Folder path to the h.vcf files of the pangenome database>  #Built with phg create-maf-vcf
--reference-hvcf // -ref <Reference.h.vcf>  #Built with phg create-ref-vcf
--chromosome // -chr <chrX> i.e. [chr1], [chr22]
--reference-fasta // -fa <Reference.fasta>  #Built with phg prepare-assemblies
--region // -reg <START-END> i.e. [10000-20000]   #If not added, whole chr is plotted
```
![FULL_chr_plot](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/pangenome_FULL_chr5.png)
![Region_to_plot](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/pangenome_chr7_575000-745000.png)


- #### Check haplotypes for a region
- #### Check identity against pangenome
- #### Plot imputed genome






## References
The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation  https://doi.org/10.1093/bioinformatics/btac410
