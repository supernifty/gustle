
## fastaquery

Given a set of fasta sequences and a set of genomes, 
finds exact matches for the fasta sequence in each genome.

### Usage

```
fastaquery [--readlength rl --mismatches mm] sequences.fa.gz genome.fa genome2.fa ...
```

e.g.
```
./fastaquery data/test_query.fa.gz data/test_genome.fa data/test_genome_2.fa
```


### Implementation

* hash the first 16 bases of each sequence
* iterate over each genome
* look for matching hash
* check that rest of sequence matches

### Installation

* download the binaries from the releases tab on github.
