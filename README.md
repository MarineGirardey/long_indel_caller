# LONG INDEL DETECTION.

**Institut Curie - Long indel detection application**


### Introduction
This application has been developed to detect long INDELs in tumor's DNA in cancer patient samples after a panel
 sequencing. The application take in input the BAM file to analyse, associated FastQ to compute the Variant Allele
Frenquency (VAF) and reference genome.

### Requirements
To run the analysis, packages have to be installed on your machine. Run this following command:

```bash 
pip install blat
pip install samtools
pip install bedtools
pip install pysam
pip install pyfastx
pip install os
pip install pandas
pip install collections
pip install heapq
pip install operator
```

### How to use LID ? - step by step

1) Run the script ```appli_debut_launcher.py```

2) Choose an empty output folder where all output files will be created.

3) Give input files : 1 BAM, 2 FASTQ, 1 reference FASTA with explicit names.

4) Define your parameters. All parameters are needed.

### Example
```bash
$ python3 appli_debut_launcher.py -fp your_chr_position -o output_folder -b your_bam_file.bam -r refrence_genome.fa -fq1 your_fastq.R1.fastq -fq2 your_fastq.R2.fastq
```

### Quick Help
```bash
$ python3 appli_debut_launcher.py - h

usage: appli_debut_launcher.py [-h] [-sp span] [-mq map_qua] [-bsp min_block_size] [-mbs span] [-nbp nb_pos] -fp filter_position -o output -b BAM -r ref -fq1 fastq1 -fq2 fastq2

Detect long INDELs in a BAM file an return a informations about the INDEL and the Variant Allele Frequency (VAF).

optional arguments:
  -h, --help            show this help message and exit
  -sp span, --span span
                        Position span to filter BAM (default: 200)
  -mq map_qua, --map_quality map_qua
                        Minimum mapping quality for the filtering part (default: 10)
  -bsp min_block_size, --blocks_span min_block_size
                        Block span to map the break point (default: 15)
  -mbs span, --min_len_blocks_size span
                        Min blocks size length acceptable (default: 5)
  -nbp nb_pos, --nb_pos_analyse nb_pos
                        Number of position to analyse (default: 1)
  -fp filter_position, --filter_pos filter_position
                        Position to filter the BAM file
  -o output, --output_folder output
                        PPath to the output folder
  -b BAM, --bam BAM     Path to the bam file
  -r ref, --reference ref
                        Path to the reference genome
  -fq1 fastq1, --fastq1 fastq1
                        Path to the fastq sample
  -fq2 fastq2, --fastq2 fastq2
                        Path to the fastq sample
```

### Run unit Tests
```bash
python -m pytest unit_test_lic.py
```
