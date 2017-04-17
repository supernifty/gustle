
## gustle

Given a set of fasta sequences and a set of genomes, 
finds exact matches for the fasta sequence in each genome.

### Usage

List all available organisms for gMLST:
```
gustle list
```

Specify organisms to download or update (downloads set of sequences and cgST definitions):
```
gustle get organism1 [organism2...]
```

Index specified organisms:
```
gustle index organism1 [organism2...]
```

Perform gMLST genotyping with the specified genome FASTA file against the specified organism:
```
gustle genotype [--readlength rl --mismatches mm --cgst cgst_file] alleles.fq.gz genome1.fa [genome2.fa...]
```

e.g.
```
./gustle genotype data/test_query.fa.gz data/test_genome.fa data/test_genome_2.fa
```

### Output Format
The output is TSV, one line per genome.
Each line contains the following columns:
* Genome filename
* cgST
* Each gene name. Each allele in the gene that are present in the genome will be listed here. If no alleles for that gene are present, "N" is shown.

The gene names are the sorted list of genes in alleles.fq.gz, unless a CGST file is specified, in which case the list of genes is taken from that file.

### Implementation

* Hash the first 16 bases of each sequence
* Iterate over each genome
* Look for matching hash
* Check that rest of sequence matches

### Installation

* Download the binaries from the releases tab on github.

### Authors

* Peter Georgeson
* Torsten Seemann
* Bernie Pope
