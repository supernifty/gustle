
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

Index specified organisms (writes alleles.fq.gz.gus):
```
gustle index --output alleles.gus alleles.fq.gz 
```

Perform gMLST genotyping with the specified genome FASTA file against the specified organism:
```
gustle genotype [--verbose] --index alleles.gus genome1.fa [genome2.fa...]
```

e.g.
```
./gustle index --cgst data/test.cgst --output data/test_query.gus data/test_query.fa.gz 
./gustle genotype --index data/test_query.gus data/test_cgst.fa
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

### TODO
* equality check should support multiple allele match
* command line argument for dealing with missing alleles N or 0?
* populate cgst
* list, get, index subcommands
* tests
* tests
