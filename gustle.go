// gustle performs a genome wide MLST
package main

import (
	"bytes"
	"flag"
	"fmt"
	"hash/fnv"
	"log"
	"os"
	"strings"
	"time"
)

const (
	VERSION           = "0.1"
	DEFAULT_SEED_SIZE = 16
	NO_CGST           = "NA"
)

type QueryPos struct {
	name    int    // array index to name of query
	pos     int    // position of kmer in query
	content string // entire query sequence
}

// list of matching alleles
type AlleleResult map[string]bool

// maps gene to list of alleles
type GenomeAlleleResult map[string]AlleleResult

func checkResult(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

func logm(level string, msg string) {
	fmt.Fprintf(os.Stderr, "%s: %s: %s\n", time.Now().String(), level, msg)
}

// stringToHash gives a 32-bit hash of a string
func stringToHash(s string) uint32 {
	h := fnv.New32a()
	h.Write([]byte(s))
	return h.Sum32()
}

// addSequence adds a sequence to a kmer hash
func addSequence(content string, sequence int, hash map[uint32][]QueryPos, seedSize int) {
	// for an exact match we only have to hash the start of the query
	if len(content) < seedSize {
		logm("WARN", fmt.Sprintf("sequence %v is length %v, shorter than seed size %v", sequence, len(content), seedSize))
		return
	}
	position := 0
	kmer := content[position : position+seedSize]
	kmerHash := stringToHash(kmer)
	entry := QueryPos{name: sequence, pos: position, content: content}
	query, ok := hash[kmerHash]
	if !ok {
		query = make([]QueryPos, 0)
	}
	query = append(query, entry)
	hash[kmerHash] = query
}

func searchSequence(content string, hash map[uint32][]QueryPos, sequenceName string, genomeFilename string, sequenceNames []string, seedSize int, result GenomeAlleleResult) {
	// populates result with a map from gene name to list of found alleles
	for position := 0; position <= len(content)-seedSize; position += 1 {
		kmer := content[position : position+seedSize]
		kmerHash := stringToHash(kmer)
		query, ok := hash[kmerHash]
		if ok {
			for _, candidate := range query {
				// check for a match
				if position+len(candidate.content) <= len(content) && content[position:position+len(candidate.content)] == candidate.content {
					// it's a match, add to result
					geneAllele := strings.Split(sequenceNames[candidate.name], "_")
					alleles, ok := result[geneAllele[0]]
					if !ok {
						alleles = make(AlleleResult, 0)
					}
					alleles[geneAllele[1]] = true
					result[geneAllele[0]] = alleles
				}
			}
		}
	}
}

func joinAlleles(alleles AlleleResult) string {
	keys := make([]string, 0, len(alleles))
	for k := range alleles {
		keys = append(keys, k)
	}
	return strings.Join(keys, ",")
}

func writeResults(filename string, result GenomeAlleleResult, geneNames map[string]bool, cgstGeneNames []string, cgstIds map[uint32]string) {
	genomeAlleles := make([]string, 0, len(cgstGeneNames))
	for _, gene := range cgstGeneNames {
		alleles, ok := result[gene]
		if ok {
			genomeAlleles = append(genomeAlleles, joinAlleles(alleles))
		} else {
			genomeAlleles = append(genomeAlleles, "N")
		}
	}
	alleles := strings.Join(genomeAlleles, "\t")
	// see if it matches a cgstId
	cgstId, ok := cgstIds[stringToHash(alleles)]
	if !ok {
		cgstId = NO_CGST
	}
	fmt.Fprintf(os.Stdout, "%s\t%s\t%s\n", filename, cgstId, alleles) // filename, cgST, alleles
}

// reverse complements a single byte
func reverse(in byte) byte {
	result, ok := REVERSE_MAP[in]
	if !ok {
		log.Fatal("failed to reverse complement")
	}
	return result
}

// reverseComplement reverses and complements a string
func reverseComplement(in *bytes.Buffer) []byte {
	var result []byte = make([]byte, in.Len(), in.Len())
	for pos := in.Len() - 1; pos >= 0; pos-- {
		current, ok := in.ReadByte()
		checkResult(ok)
		result[pos] = reverse(current)
	}
	return result
}

var REVERSE_MAP = map[byte]byte{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

func showUsage() {
	fmt.Fprintf(os.Stderr, "gustle version %v\nUsage: gustle [list|get|index|genotype|version] [subcommand arguments]\n", VERSION)
	flag.PrintDefaults()
	os.Exit(1)
}

// version outputs the current software version
func version() {
	fmt.Fprintf(os.Stderr, "gustle version %v\n", VERSION)
}

func list(args []string) {
	fmt.Fprintf(os.Stderr, "list: TODO\n")
}

func get(args []string) {
	fmt.Fprintf(os.Stderr, "get: TODO\n")
}

func index(args []string) {
	fmt.Fprintf(os.Stderr, "index: TODO\n")
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
