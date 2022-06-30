
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

* e.g. use data/test*
```
./gustle index --cgst data/test.cgst --output data/test_query.gus data/test_query.fa.gz 
./gustle genotype --index data/test_query.gus data/test_cgst.fa
```

* e.g. use data/second-test*
```
cd ~/gustle/data
../gustle index \
   --cgst second-test.cgst \
   --output second-test_query.gus \
   second-test_query.fa
../gustle genotype \
   --index second-test_query.gus \
   second-test_cgst.fa \
   second-test_genome.fa \
   second-test_genome_2.fa > second-test.summary.tsv
cat second-test.summary.tsv
```

| Filename | cgST | geneA | geneB | geneC | geneD |
| --- | --- | --- | --- | --- |--- |
| *second-test_cgst.fa* | 8675309 (4/8) | 0 | 0 | 1 | 16 |
| *second-test_genome.fa* | 8675309 (3/14) | 0,16,32 | 0 | N | 16 |
| *second-test_genome_2.fa* | 8675309 (3/14) | 16,32,0 | 0 | N | 16 |

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
* For the latest on GitHub with new Golang without root for $USER in Linux:
    1. Set a Go version
    `GO_VERSION=1.18.3`
    1. Dowload the go tarball
    `curl -LO https://go.dev/dl/go"${GO_VERSION}".linux-amd64.tar.gz`
    1. Make paths for go to be unpacked into
    `mkdir -p "${HOME}/.go/${GO_VERSION}/go-pkgs"`
    1. Uncompress the tarball into a hidden Home path
    `tar -C "${HOME}/.go/${GO_VERSION}" -xzf go"${GO_VERSION}".linux-amd64.tar.gz`
    1. Add 4 new lines into bashrc to make this version available
    ```
    echo '# Golang paths' >> ~/.bashrc
    echo "export GOROOT="$HOME"/.go/"${GO_VERSION}"/go" >> ~/.bashrc
    echo "export GOPATH="$HOME"/.go/"${GO_VERSION}"/go-pkgs" >> ~/.bashrc
    echo 'export PATH="$PATH":"$GOROOT"/bin:"$GOPATH"/bin' >> ~/.bashrc
    ```
    1. Re-source the bashrc to make `go` available
    `source ~/.bashrc`
    1. Fetch an updated copy of the Gustle repository
    `cd $HOME && git clone git@github.com:supernifty/gustle.git`
    1. Compile the `gustle` binary
    `cd gustle && go build main/gustle.go`
    1. Optionally make a symlink to a common area of binaries
    `ln -sv ~/gustle/gustle ~/bin/gustle`


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
