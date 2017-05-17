// gustle performs a genome wide MLST
package main

import (
	"flag"
	"fmt"
	"github.com/supernifty/gustle/cgmlst"
	"os"
)

const (
	VERSION           = "0.1"
	DEFAULT_SEED_SIZE = 16
)

func showUsage() {
	fmt.Fprintf(os.Stderr, "gustle version %v\nUsage: gustle [list|get|index|genotype|version] [subcommand arguments]\n", VERSION)
	flag.PrintDefaults()
	os.Exit(1)
}

// version outputs the current software version
func version() {
	fmt.Fprintf(os.Stderr, "gustle version %v\n", VERSION)
}

// list shows the available cgst organisms
func list(args []string) {
	fmt.Fprintf(os.Stderr, "list: TODO\n")
}

// get downloads the specified organisms
func get(args []string) {
	fmt.Fprintf(os.Stderr, "get: TODO\n")
}

// indexes an organism
func index(args []string) {
	indexCommand := flag.NewFlagSet("index", flag.ExitOnError)
	var seedSize int
	var verbose bool
	indexCommand.IntVar(&seedSize, "readlength", 16, "minimum read length of queries (16)")
	indexCommand.BoolVar(&verbose, "verbose", false, "include additional logging")
	indexCommand.Parse(args)
	if !indexCommand.Parsed() || indexCommand.NArg() < 1 {
		fmt.Fprintf(os.Stderr, "gustle version %v\nUsage: gustle index queries.fa.gz\n", VERSION)
		os.Exit(1)
	}
	sequencesFilename := indexCommand.Arg(0)
	var sequenceNames []string
	var geneNames map[string]bool = make(map[string]bool)
	// queryKmers := genotype.IndexSequences(sequencesFilename, &sequenceNames, geneNames, seedSize, verbose)
	genotype.IndexSequences(sequencesFilename, &sequenceNames, geneNames, seedSize, verbose)
}

// genotype finds matching alleles on a specified genome
func doGenotype(args []string) {
	genotypeCommand := flag.NewFlagSet("genotype", flag.ExitOnError)
	var seedSize int
	var mismatches int
	var cgst string
	var verbose bool
	genotypeCommand.IntVar(&seedSize, "readlength", 16, "minimum read length of queries (16)")
	genotypeCommand.IntVar(&mismatches, "mismatches", 0, "mismatches to include (0)")
	genotypeCommand.StringVar(&cgst, "cgst", "", "cgST file")
	genotypeCommand.BoolVar(&verbose, "verbose", false, "include additional logging")
	genotypeCommand.Parse(args)
	if !genotypeCommand.Parsed() || genotypeCommand.NArg() < 2 {
		fmt.Fprintf(os.Stderr, "gustle version %v\nUsage: gustle genotype [-kmer kmer -mismatches mismatches -cgst cgst_file -verbose] queries.fa.gz genome.fa [genome.fa...]\n", VERSION)
		os.Exit(1)
	}
	sequencesFilename := genotypeCommand.Arg(0)

	genotype.FindAlleles(seedSize, mismatches, cgst, sequencesFilename, genotypeCommand.Args()[1:], verbose)
}

func main() {
	// parse command line arguments
	if len(os.Args) == 1 {
		showUsage()
	}

	switch os.Args[1] {
	case "list":
		list(os.Args[2:])
	case "get":
		get(os.Args[2:])
	case "index":
		index(os.Args[2:])
	case "genotype":
		doGenotype(os.Args[2:])
	case "version":
		version()
	default:
		showUsage()
	}
}
