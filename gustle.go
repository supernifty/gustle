// gustle performs a genome wide MLST
package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"flag"
	"fmt"
	"hash/fnv"
	"log"
	"os"
	"sort"
	"strings"
	"time"
)

const (
	VERSION           = "0.1"
	DEFAULT_SEED_SIZE = 16
)

type QueryPos struct {
	name    int    // array index to name of query
	pos     int    // position of kmer in query
	content string // entire query sequence
}

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

func searchSequence(content string, hash map[uint32][]QueryPos, sequenceName string, genomeFilename string, sequenceNames []string, seedSize int, result map[string][]string) {
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
						alleles = make([]string, 0)
					}
					alleles = append(alleles, geneAllele[1])
					result[geneAllele[0]] = alleles
				}
			}
		}
	}
}

func writeResults(filename string, result map[string][]string, geneNames map[string]bool, cgstGeneNames []string) {
	genomeAlleles := make([]string, 0, len(cgstGeneNames))
	for _, gene := range cgstGeneNames {
		alleles, ok := result[gene]
		if ok {
			alleleList := strings.Join(alleles, ",")
			genomeAlleles = append(genomeAlleles, alleleList)
		} else {
			genomeAlleles = append(genomeAlleles, "N")
		}
	}
	fmt.Fprintf(os.Stdout, "%s\t%s\t%s\n", filename, "TODO", strings.Join(genomeAlleles, "\t")) // filename, cgST, alleles
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
	logm("INFO", fmt.Sprintf("processing with seed size %v, up to %v mismatches, cgst: %s: %s", seedSize, mismatches, cgst, sequencesFilename))

	// open fa.gz file
	file, err := os.Open(sequencesFilename)
	checkResult(err)
	defer file.Close()

	gr, err := gzip.NewReader(file)
	checkResult(err)
	defer gr.Close()

	scanner := bufio.NewScanner(gr)
	lines := 0
	sequenceCount := 0

	var sequenceNames []string
	var geneNames map[string]bool = make(map[string]bool)
	var queryKmers map[uint32][]QueryPos = make(map[uint32][]QueryPos)
	var content *bytes.Buffer

	// index sequences
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if content != nil {
				addSequence(content.String(), sequenceCount-1, queryKmers, seedSize)
			}
			sequenceNames = append(sequenceNames, line[1:])
			sequenceCount++
			content = new(bytes.Buffer)
			geneName := strings.Split(line[1:], "_")[0]
			_, ok := geneNames[geneName]
			if !ok {
				geneNames[geneName] = true
			}
		} else {
			(*content).WriteString(line)
		}
		lines++
		if lines%1000000 == 0 {
			logm("INFO", fmt.Sprintf("processing: %s: %v lines %v sequences. %v kmers", sequencesFilename, lines, len(sequenceNames), len(queryKmers)))
		}
	}
	addSequence(content.String(), sequenceCount-1, queryKmers, seedSize)

	var cgstGeneNames = make([]string, 0, len(geneNames))
	if cgst == "" {
		for gene := range geneNames {
			cgstGeneNames = append(cgstGeneNames, gene)
		}
		sort.Strings(cgstGeneNames)
	} else {
		// read cgst details
		file, err := os.Open(cgst)
		checkResult(err)
		defer file.Close()
		r := bufio.NewReader(file)
		scanner := bufio.NewScanner(r)
		first := true
		for scanner.Scan() {
			if first {
				first = false
				fields := strings.Split(scanner.Text(), "\t")
				for _, gene := range fields[1:] {
					cgstGeneNames = append(cgstGeneNames, gene)
				}
			} else {
				// TODO
			}
		}
	}

	logm("INFO", fmt.Sprintf("done processing: %s: %v lines %v sequences.", sequencesFilename, lines, len(sequenceNames)))

	// genotype each genome - list matching sequences
	fmt.Fprintf(os.Stdout, "Filename\tcgST\t%v\n", strings.Join(cgstGeneNames, "\t"))
	for _, genomeFilename := range genotypeCommand.Args()[1:] {
		logm("INFO", fmt.Sprintf("processing genome: %s", genomeFilename))
		var result map[string][]string = make(map[string][]string)
		file, err := os.Open(genomeFilename)
		checkResult(err)
		defer file.Close()

		lines = 0
		sequenceCount = 0
		content = new(bytes.Buffer)
		r := bufio.NewReader(file)
		scanner := bufio.NewScanner(r)
		var sequenceName string
		for scanner.Scan() {
			line := scanner.Text()
			if strings.HasPrefix(line, ">") {
				if content != nil {
					searchSequence(content.String(), queryKmers, sequenceName, genomeFilename, sequenceNames, seedSize, result)
					searchSequence(string(reverseComplement(content)), queryKmers, sequenceName, genomeFilename, sequenceNames, seedSize, result)
				}
				sequenceName = line[1:]
				sequenceCount++
				content = new(bytes.Buffer)
			} else {
				(*content).WriteString(line)
			}
			lines++
		}
		searchSequence(content.String(), queryKmers, sequenceName, genomeFilename, sequenceNames, seedSize, result)
		searchSequence(string(reverseComplement(content)), queryKmers, sequenceName, genomeFilename, sequenceNames, seedSize, result)
		logm("INFO", fmt.Sprintf("done genome: %s. %v lines. %v sequences.", genomeFilename, lines, sequenceCount))
		writeResults(genomeFilename, result, geneNames, cgstGeneNames)
	}
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
