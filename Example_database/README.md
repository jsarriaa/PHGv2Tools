To guide any user about the usage of phgtools an example set is been prepared. It uses same files as [PHG testing kit]
It actually takes: 
- LineA.fa, LineB.fa, LineC.fa
- Ref.fa
- anchors.gff

*** Ensure to have the proper phgtools conda and intalled the phg software

It will give a 3 genomes pangenome using a reference genome and its gene annotation.
There is also a raw fastq (doomie) imitating a low density sequence.

This is the initial matherial to start the pangenome:
```
phg_v2_example/
├── data
│   ├── anchors.gff
│   ├── Ref-v5.fa
│   ├── LineA.fa
│   ├── LineB.fa
│   ├── LineC.fa
│   └── LineC.fastq
└── output
    └── vcf_files
```
After installing PHG, ensure to [initiate the database](https://github.com/maize-genetics/phg_v2?tab=readme-ov-file#build-and-load-data) and to [prepare properly the fastas](https://phg.maizegenetics.net/build_and_load/#prepare-assembly-fasta-files).

```
./phg initdb --db-path /path/to/dbs
./phg prepare-assemblies --keyfile /path/to/keyfile --output-dir data/updated_assemblies --threads numberThreadstoRun
```
```
phg_v2_example/
├── data
│   └── prepared_genomes
│   │   ├── Ref-v5.fa
│   │   ├── LineA.fa
│   │   ├── LineB.fa
│   │   └── LineC.fa
│   ├── anchors.gff
│   └── LineC.fastq
└── output
│    └── vcf_files
└── vcf_dbs
    └── gvcf_dataset
    └── hvcf_dataset
    └── temp
```
Then, prepare the ranges reference file, align and compress the genomes to include at the database. Follow PHG commands
```
./phg create-ranges --reference-file data/updated_assemblies/Ref.fa --gff my.gff --boundary gene --pad 500 --range-min-size 500 -o /path/to/bed/file.bed
./phg align-assemblies --gff anchors.gff --reference-file data/updated_assemblies/Ref.fa --assembly-file-list assembliesList.txt --total-threads 20 --in-parallel 4 -o /path/for/generatedFiles
./phg agc-compress --db-path /path/to/dbs --reference-file data/updated_assemblies/Ref.fa --fasta-list /my/assemblyFastaList.txt
```
Directory now should be:
```
phg_v2_example/
├── data
│   └── prepared_genomes
│   │   ├── Ref-v5.fa
│   │   ├── LineA.fa
│   │   ├── LineB.fa
│   │   └── LineC.fa
│   ├── anchors.gff
│   └── LineC.fastq
└── output
│   │   └── vcf_files
│   ├── anchorwave_gff2seq_error.log *
│   ├── anchorwave_gff2seq_output.log *
│   ├── LineA.maf *
│   ├── LineA.sam *
│   ├── LineA.svg *
│   ├── LineA_Ref.anchorspro *
│   ├── LineB.maf *
│   ├── LineB.sam *
│   ├── LineB.svg *
│   ├── LineB_Ref.anchorspro *
│   ├── LineC.maf *
│   ├── LineC.sam *
│   ├── LineC.svg *
│   ├── LineC_Ref.anchorspro *
│   ├── minimap2_LineA_error.log *
│   ├── minimap2_LineA_output.log *
│   ├── minimap2_LineB_error.log *
│   ├── minimap2_LineB_output.log *
│   ├── minimap2_LineC_error.log *
│   ├── minimap2_LineC_output.log *
│   ├── minimap2_Ref_error.log *
│   ├── minimap2_Ref_output.log *
│   ├── proali_LineA_outputAndError.log *
│   ├── proali_LineB_outputAndError.log *
│   ├── proali_LineC_outputAndError.log *
│   ├── ref.cds.fasta *
│   └── Ref.sam *
└── vcf_dbs
    ├── assemblies.agc*
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── reference *
    │   ├── Ref.bed *
    │   └── Ref.fa *
    └── temp/
```
Then, proceed creating the VCF files
```
./phg create-ref-vcf --bed /my/bed/file.bed --reference-file data/updated_assemblies/Ref.fa --reference-url https://url-for-ref --reference-name B73 --db-path /path/to/tiled/dataset folder
./phg create-maf-vcf --db-path /path/to/dbs --bed /my/bed/file.bed --reference-file data/updated_assemblies/Ref.fa --maf-dir /my/maf/files -o /path/to/vcfs
```
New files are generated.
```
phg_v2_example/
├── data
│   └── prepared_genomes
│   │   ├── Ref-v5.fa
│   │   ├── LineA.fa
│   │   ├── LineB.fa
│   │   └── LineC.fa
│   ├── anchors.gff
│   └── LineC.fastq
└── output
│   │   └── vcf_files
│   │       ├── LineA.h.vcf.gz *
│   │       ├── LineA.h.vcf.gz.csi *
│   │       ├── LineA.g.vcf.gz *
│   │       ├── LineA.g.vcf.gz.csi *
│   │       ├── LineB.h.vcf.gz *
│   │       ├── LineB.h.vcf.gz.csi *
│   │       ├── LineB.g.vcf.gz *
│   │       ├── LineB.g.vcf.gz.csi *
│   │       ├── LineC.h.vcf.gz *
│   │       ├── LineC.h.vcf.gz.csi *
│   │       ├── LineC.g.vcf.gz *
│   │       ├── LineC.g.vcf.gz.csi 
│   │       └── VCFMetrics.tsv *
│   ├── anchorwave_gff2seq_error.log 
│   ├── anchorwave_gff2seq_output.log 
│   ├── LineA.maf 
│   ├── LineA.sam 
.   .
.   .
.   .
│   └── Ref.sam 
└── vcf_dbs
    ├── assemblies.agc
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── hvcf_files *
    │   ├── Ref.h.vcf.gz *
    │   └── Ref.h.vcf.gz.csi *
    ├── reference 
    │   ├── Ref.bed 
    │   └── Ref.fa 
    └── temp/
```
Then, building the kmer index and prepare for the imputation of the test fastq file
```
./phg build-kmer-index --db-path /my/db/uri --hvcf-dir /my/hvcf/dir
./phg map-kmers --hvcf-dir /my/hvcf/dir --kmer-index /my/hvcf/dir/kmerIndex.txt --key-file /my/path/keyfile --output-dir /my/mapping/dir
./phg find-paths --path-keyfile /my/path/keyfile --hvcf-dir /my/hvcf/dir --reference-genome /my/ref/genome --path-type haploid --output-dir /my/imputed/hvcfs
```
After imputation, the new LineD.h.vcf represents the imputed haplotype.
```
phg_v2_example/
├── data
│   └── prepared_genomes
│   │   ├── Ref-v5.fa
│   │   ├── LineA.fa
│   │   ├── LineB.fa
│   │   └── LineC.fa
│   ├── anchors.gff
│   └── LineC.fastq
└── output
│   │   └── vcf_files
│   │       ├── LineA.h.vcf.gz 
│   │       ├── LineB.h.vcf.gz 
│   │       ├── LineC.h.vcf.gz 
│   │       .
│   │       .
│   │       ├── kmerIndex.txt *
│   │       ├── kmerIndexStatistics.txt *
│   │       ├── LineD_readMapping.txt *
│   │       ├── LineD_imputed.h.vcf *
│   │       └── VCFMetrics.tsv 
│   ├── anchorwave_gff2seq_error.log 
│   ├── anchorwave_gff2seq_output.log 
│   ├── LineA.maf 
│   ├── LineA.sam 
.   .
.   .
.   .
│   └── Ref.sam 
└── vcf_dbs
    ├── assemblies.agc
    ├── gvcf_dataset/
    ├── hvcf_dataset/
    ├── hvcf_files *
    │   ├── Ref.h.vcf.gz *
    │   └── Ref.h.vcf.gz.csi *
    ├── reference 
    │   ├── Ref.bed 
    │   └── Ref.fa 
    └── temp/
```
Last step, to merge all hvcfs of Pangenome (LineA, LineB, LineC) in a single file using:
```
phg merge-hvcfs --input-dir my/hvcf/directory --output-file output/merged_hvcfs.h.vcf --id-format CHECKSUM --reference-file --range-bedfile
```

This is then the files to have in mind for phgtools analysis:

```
├── data
│   └── prepared_genomes
└── output
│   │   └── vcf_files
│   │       ├── LineA.h.vcf.gz ***
│   │       ├── LineB.h.vcf.gz ***
│   │       ├── LineC.h.vcf.gz ***
│   │       ├── LineD_imputed.h.vcf ***
│   │       └── MergedLinesA_B_C.h.vcf ***
└── vcf_dbs
    └── hvcf_files 
            └──── Ref.h.vcf.gz ***

```
## PHG tools utilities:

#### Pangenome ranges evolution
```
range-pangenome-evolution --hvcf-folder /Example_database/output/vcf_files/ --reference-file Ref.fa --range-bedfile output/ref_ranges.bed
```
The testing files are small, and new ranges are not included after including new genomes
[RangesAmplificationSlope](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/RangesAmplificationSlope.png)


#### Core, accesory and unique ranges
```
core-range-detecter --pangenome-hvcf output/MergedLinesA_B_C.h.vcf
```
These pangenome has either ranges in all haplotypes or only in one, as shown in plots:

[MergedLinesA_B_C.2.png](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/MergedLinesA_B_C.h.vcf_2.png)
[MergedLinesA_B_C.2.png](https://github.com/jsarriaa/PHGv2Tools/blob/main/Misc/Images/MergedLinesA_B_C.h.vcf_1.png)


