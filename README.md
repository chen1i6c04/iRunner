# iRunner
A simple bacterial genome assembly pipeline for Illumina paired-end reads from NCBI SRA run

## Requirements
* Python3.6+
* [sra-tools](https://github.com/ncbi/sra-tools)
  * fastq-dump
  * sra-stat
* [fastp](https://github.com/OpenGene/fastp)
* [shovill](https://github.com/tseemann/shovill)

## Quick Start
```
irunner.py --outdir [path/of/output] --threads [number] [SRA accession]
```
## Usage
```
arguments:
  accession   NCBI SRA accession
  --outdir    Output folder
  --tmpdir    Directory of temp folder default: '/tmp'
  --gsize     Estimated genome size(MB) eg. 3.2M, If blank will AUTODETECT. default: ''
  --threads   Number of threads to use. default: 8
```
