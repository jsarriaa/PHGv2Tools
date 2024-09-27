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
│   ├── ref_ranges.bed
│   ├── assemblies.agc
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
└── gvcf_dataset
└── hvcf_dataset
└── temp
```
Then, proceed creating the VCF files
```
./phg create-ref-vcf --bed /my/bed/file.bed --reference-file data/updated_assemblies/Ref.fa --reference-url https://url-for-ref --reference-name B73 --db-path /path/to/tiled/dataset folder
./phg create-maf-vcf --db-path /path/to/dbs --bed /my/bed/file.bed --reference-file data/updated_assemblies/Ref.fa --maf-dir /my/maf/files -o /path/to/vcfs
```
