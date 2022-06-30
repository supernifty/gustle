// gustle performs a genome wide MLST
package main

import (
	"flag"
	"fmt"
	"github.com/supernifty/gustle/cgmlst"
	"os"
)

const (
	VERSION           = "0.2.1"
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
// gustle index --output x --cgst cgst inputs
func index(args []string) {
	indexCommand := flag.NewFlagSet("index", flag.ExitOnError)
	var seedSize int
	var verbose bool
	var cgst string
	var output string
	indexCommand.IntVar(&seedSize, "readlength", 16, "minimum read length of queries (16)")
	indexCommand.BoolVar(&verbose, "verbose", false, "include additional logging")
	indexCommand.StringVar(&cgst, "cgst", "", "cgST input file")
	indexCommand.StringVar(&output, "output", "", "write index to this file")
	indexCommand.Parse(args)
	if !indexCommand.Parsed() || cgst == "" || indexCommand.NArg() < 1 {
		fmt.Fprintf(os.Stderr, "gustle version %v\nUsage: gustle index --output output_file --cgst cgst_file queries.fa.gz [queries.fa.gz...]\n", VERSION)
		os.Exit(1)
	}
	db := genotype.IndexSequences(indexCommand.Args(), cgst, seedSize, verbose)
	genotype.SaveIndex(output, db, verbose)
}

// genotype finds matching alleles on a specified genome
func doGenotype(args []string) {
	genotypeCommand := flag.NewFlagSet("genotype", flag.ExitOnError)
	var mismatches int
	var verbose bool
	var index string
	genotypeCommand.IntVar(&mismatches, "mismatches", 0, "mismatches to include (0)")
	genotypeCommand.BoolVar(&verbose, "verbose", false, "include additional logging")
	genotypeCommand.StringVar(&index, "index", "", "index file")
	genotypeCommand.Parse(args)
	if !genotypeCommand.Parsed() || index == "" || genotypeCommand.NArg() < 1 {
		fmt.Fprintf(os.Stderr, "gustle version %v\nUsage: gustle genotype [-mismatches mismatches -verbose] --index index_file genome.fa [genome.fa...]\n", VERSION)
		os.Exit(1)
	}

	var db genotype.QueryList = genotype.LoadIndex(index, verbose)
	genotype.FindAlleles(db, mismatches, genotypeCommand.Args(), verbose)
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
