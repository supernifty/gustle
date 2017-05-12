// gustle performs a genome wide MLST
package main

import (
	"flag"
	"fmt"
	"log"
	"os"
	"time"
)

const (
	VERSION           = "0.1"
	DEFAULT_SEED_SIZE = 16
	NO_CGST           = "NA"
)

func checkResult(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

func logm(level string, msg string) {
	fmt.Fprintf(os.Stderr, "%s: %s: %s\n", time.Now().String(), level, msg)
}

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
	fmt.Fprintf(os.Stderr, "index: TODO\n")
}

// genotype finds matching alleles on a specified genome
func genotype(args []string) {
	genotypeCommand := flag.NewFlagSet("genotype", flag.ExitOnError)
	var seedSize int
	var mismatches int
	var cgst string
	genotypeCommand.IntVar(&seedSize, "readlength", 16, "minimum read length of queries (16)")
	genotypeCommand.IntVar(&mismatches, "mismatches", 0, "mismatches to include (0)")
	genotypeCommand.StringVar(&cgst, "cgst", "", "cgST file")
	genotypeCommand.Parse(args)
	if !genotypeCommand.Parsed() || genotypeCommand.NArg() < 2 {
		fmt.Fprintf(os.Stderr, "gustle version %v\nUsage: gustle genotype [-kmer kmer -mismatches mismatches -cgst cgst_file] queries.fa.gz genome.fa [genome.fa...]\n", VERSION)
		os.Exit(1)
	}
	sequencesFilename := genotypeCommand.Arg(0)

	findAlleles(seedSize, mismatches, cgst, sequencesFilename, genotypeCommand.Args()[1:])
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
		genotype(os.Args[2:])
	case "version":
		version()
	default:
		showUsage()
	}
}
