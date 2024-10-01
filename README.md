# PHGv2Tools
Package to downstream analysis of pangenomes databases, working with Practical Haplotype Graph and its h.VCF files.
Note: This is an early release, and still under development. Not proper tests have been developed. Better annotations and code improvement will be implemented.
For anny comment, feel free to contact, any feedback is welcome:
jsarria@eead.csic.es

## Contents
- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)
- [Notebooks](#notebooks)
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

#### Quick Start
```
range-pangenome-evolution --hvcf-folder /Example_database/output/vcf_files/ --reference-file Ref.fa --range-bedfile output/ref_ranges.bed

core-range-detecter --pangenome-hvcf output/MergedLinesA_B_C.h.vcf

plot-pangenome-chromosomes --hvcf-folder output/vcf_files/ --reference-hvcf output/Ref.h.vcf.gz --chromosome chr2 --region 15000-35000 --reference-fasta Ref.fa

check-haplotype-alleles --hvcf output/MergedLinesA_B_C.h.vcf --reference-fasta Ref.fa --start 18800 --end 20100 --chromosome 2

check-imputated-haplotype --hvcf-folder output/ --hvcf-file output/LineD.h.vcf

plot-imputed-hvcf --input-hvcf output/LineD.h.vcf --pangenome-hvcf-folder output/ --reference-hvcf hvcf_files/Ref.h.vcf.gz
```


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
  --pangenome-hvcf // -hvcf  <Merged of all hvcf of the pangenome database>  #Built with phg merge-hvcfs 
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
For those who are looking for genes alleles; Ranges are built with genes annotation. Having the coordinates at the reference genome and running ```check-haplotype-alleles``` will provide a list of the haplotype ranges found in each genome and its coordinates of the ensambled respective genome.

These are the mandatory arguments:
```
--hvcf // -hf   <Merged of all hvcf of the pangenome database>  #Built with phg merge-hvcfs 
--reference-fasta // -ref   <Reference.fasta>  #Built with phg prepare-assemblies
```
And then, provide the coordinates. These are optional, but will be asked ahead if there are not as a command line argument.
```
--chromosome // -c   <INT>  i.e. <2>
--start // -s <START>   i.e. <1000>
--end // -e  <END>  i.e. <5000>
```
As output get:
```
Start:  101659993
End:  101663495
Chromosome:  chr3 

there are 10 keys of ranges containing the coordinates:
['fab0e96084dc1cd37d85e739d9142fe1', '44bd8430f89cdc117f9ff4e5686dc16a', '02d192f4d99f589c987a3d4b70688fac', '6c54df6405d0a54f48805e72edfab1bc', '15ef9d82eb9d5b3e58c8427c26800e8a', '7ce1fef0bdcce3f10d4582b6ccee55c3', '4f38def45958aff23b026a8bc092971d', 'ec5cff96aa97d55a83311afe7bf6aef8', 'a44306a57bbbb7b59741e6cbe324e9d0', '06af94fcc1afa0679b88c7053f1251f9']


Line: HOR_12184 Region: chr3H_OX460224.1:103196667-103201170    Reference range:chr3H_LR890098.1:101659493-101663995> 
Line: HOR_3474  Region: chr3H:104111380-104115859               Reference range:chr3H_LR890098.1:101659493-101663995> 
Line: HOR_14121 Region: chr3H:105286336-105290830               Reference range:chr3H_LR890098.1:101659493-101663995> 
Line: HOR_1168  Region: chr3H_OX460231.1:105381557-105386054    Reference range:chr3H_LR890098.1:101659493-101663995> 
Line: HOR_21599 Region: chr3H_OX460266.1:105308925-105313416    Reference range:chr3H_LR890098.1:101659493-101663995> 
Line: HOR_13942 Region: chr3H_OX460097.1:102757692-102762177    Reference range:chr3H_LR890098.1:101659493-101663995> 
Line: HOR_21595 Region: chr3H:103459448-103463939               Reference range:chr3H_LR890098.1:101659493-101663995> 
Line: HOR_3365  Region: chr3H_OX460308.1:103739508-103744008    Reference range:chr3H_LR890098.1:101659493-101663995> 
Line: HOR_2779  Region: chr3H_OX460210.1:101243522-101248011    Reference range:chr3H_LR890098.1:101659493-101663995> 
Line: HOR_10892 Region: chr3H_OX460112.1:105761775-105766280    Reference range:chr3H_LR890098.1:101659493-101663995> 
```
*** If instead a merged VCF from the whole pangenome you give as imput an imputed VCF file, it will return the presence (or not) of a range and to which pangenome haplotype has been associated.


- #### Plot imputed genome
Taking an imputed h.vcf from ```phg map-kmers``` and ```find-paths``` and plots it in a single image:
![Imputed_haplo_gdb136](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/1740D-268-01_S1_L001_R1_001_GDB136.h.vcf.png)

```plot-imputed-hvcf``` takes as argument:
```
--input-hvcf // -hvcf <Imputed h.vcf>  #Built with phg find-paths
--reference-hvcf // -ref <Reference.h.vcf>  #Built with phg create-ref-vcf
--pangenome-hvcf-folder // -folder  <Folder path to the h.vcf files of the pangenome database>  #Built with phg create-maf-vcf
```



- #### Check identity against pangenome
Reads all pangenome haplotypes and compair with an imputated h.vcf, to extract the percentage of the genome is used to build up the imputed haplotype. It is now only based in nº ranges itself, not based in base pairs absolute ammount. Call the function with ```check-imputated-haplotype``` and provide as arguments:
```
--hvcf-folder // -folder  <Folder path to the h.vcf files of the pangenome database>  #Built with phg create-maf-vcf
--hvcf-file // -file <Imputed h.vcf>  #Built with phg find-paths
```
Working on ploting the results, for now we get a list for the results:
```
Total ranges in ../output/ensambled_genomes/1740D-268-01_S1_L001_R1_001_GDB136.h.vcf is 65703

Checking HOR_14121
Match count for HOR_14121 is 15678 out of 65703 ranges(23.86%)

Checking HOR_2830
Match count for HOR_2830 is 15942 out of 65703 ranges(24.26%)

Checking HOR_2779
Match count for HOR_2779 is 14736 out of 65703 ranges(22.43%)

Checking HOR_12184
Match count for HOR_12184 is 22570 out of 65703 ranges(34.35%)

Checking HOR_3474
Match count for HOR_3474 is 13813 out of 65703 ranges(21.02%)

Checking HOR_21599
Match count for HOR_21599 is 14722 out of 65703 ranges(22.41%)

Checking HOR_1168
Match count for HOR_1168 is 21094 out of 65703 ranges(32.11%)

Checking HOR_13942
Match count for HOR_13942 is 23469 out of 65703 ranges(35.72%)

Checking HOR_21595
Match count for HOR_21595 is 14884 out of 65703 ranges(22.65%)

Checking HOR_10892
Match count for HOR_10892 is 17269 out of 65703 ranges(26.28%)

Checking HOR_3365
Match count for HOR_3365 is 16552 out of 65703 ranges(25.19%)
```

## Notebooks
PHGtools is ready to use as command line software. However, every script is ready to use either as a command line or a notebook. Download the files from [notebooks repository](https://github.com/jsarriaa/PHGv2Tools/tree/main/Notebooks). Find there a brief guide too.

## References
The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation  https://doi.org/10.1093/bioinformatics/btac410
