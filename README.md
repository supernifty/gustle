
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
Genome filename: name of input file
Gene: allele gene name
Count: number of exactly matching alleles
Alleles: comma separated list of exactly matching alleles

### Implementation

* hash the first 16 bases of each sequence
* iterate over each genome
* look for matching hash
* check that rest of sequence matches

### Installation

* download the binaries from the releases tab on github.

### Authors

* Peter Georgeson
* Torsten Seemann
* Bernie Pope
