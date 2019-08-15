# read_counter_HTSeq

Counts paired-end fragments from SAM/BAM files for defined input regions.

This repository contains a single script [read_counter_HTSeq.py](read_counter_HTSeq.py) to count mapped reads on given input sequences.

## Dependencies

- Python 2.7.12
- Python packages:
  - argparse
  - HTSeq (https://htseq.readthedocs.io/en/release_0.11.1/)
  - collections
  

## Installation

Just download this script: [read_counter_HTSeq.py](read_counter_HTSeq.py)

## Usage

```
usage: read_counter_HTSeq.py [-h] -i [INPUT_FILES [INPUT_FILES ...]] -b
                             BED_FILE [-l [LABELS [LABELS ...]]] [-c]
                             [-ft {sam,bam}] [-r {name,pos}] -o OUTPUT_FILE

 
Counts paired-end fragments from SAM/BAM files for each region (e.g. genes, transcripts or exons) in a .bed file.

optional arguments:
  -h, --help            show this help message and exit
  -i [INPUT_FILES [INPUT_FILES ...]], --input_files [INPUT_FILES [INPUT_FILES ...]]
                        input SAM/BAM file(s).
  -b BED_FILE, --bed_file BED_FILE
                        Bed file.
  -l [LABELS [LABELS ...]], --labels [LABELS [LABELS ...]]
                        labels for each sam in same number and order as
                        input_files
  -c, --use_chrom_name  In case off mapping against transcriptome, use the
                        chrom-name as feature for counting.
  -ft {sam,bam}, --file_type {sam,bam}
                        Read-alignment file format.
  -r {name,pos}, --order {name,pos}
                        Sorting order of <alignment_file> (default: name).
                        Paired-end sequencing data must be sorted either by
                        position or by read name, and the sorting order must
                        be specified. Ignored for single-end data.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        output file in bed format.

Jonas Ibn-Salem <ibnsalem@molgen.mpg.de> 29.09.12


```
