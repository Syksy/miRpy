# miRpy - misc miRNA Python tools

## Installation

Clone with ```git``` and install using ```pip```:
```commandline
cd /my/path/
git clone git@github.com:Syksy/miRpy.git
cd miRpy
pip install -e .
mirpy
```

Key dependency is ```bamnostic``` (version >=1.1) for processing BAM-files.

## Example usage

### Download and process latest miRBase annotations 

Default download location is from 
```commandline
https://mirbase.org/download/hsa.gff3
```

to a local destination:

```commandline
mirpy download /my/path/hsa.gff3
head hsa.gff3 -n 20
```

Notice that this also includes miRNA precursors (e.g. miRNA_primary_transcript).

Extract mature miRNA sequences with from the tab-separated gff3 
file with a matching criteria from the 2nd column (has to match 'miRNA' or the
row is omitted):

```commandline
mirpy subset --in hsa.gff3 --out hsa_mature.gff3 --col 2 --criteria miRNA
head hsa_mature.gff -n 20
```

### Example dataset

```miRpy``` delivers with 3 samples with the first 1k alignments (both *.bam and *.bai):

```commandline
mirpy examples --out ./
ls sample*.bam
```

Check the headers of the BAMs for sanity checking (checking a head of the BAMs): 
```commandline
mirpy head sample*.bam
```

### Map aligned reads to the miRBase annotations

Test mapping aligned files to the miRBase mature MRI database

```commandline
mirpy map sample*.bam --gff hsa_mature.gff3 --out samples_map.tsv --metric approx --shift 4 --multi fractional
```

## Acknowledgements

<TODO>
