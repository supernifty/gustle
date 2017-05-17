package genotype

//
// this package implements
//

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"hash/fnv"
	"log"
	"os"
	"sort"
	"strings"
	"time"
)

const (
	NO_CGST = "NA"
)

// QueryIndex maps hashes of kmers to all possible source locations
type QueryIndex map[uint32][]QueryPos

// QueryPos provides an index into a sequence name, and a position
type QueryPos struct {
	name    int    // array index to name of query
	pos     int    // position of kmer in query
	content string // entire query sequence
}

// Allele is names a short sequence
type GeneName string

// AlleleResult is a list of allele names
type AlleleResult map[string]bool

// GenomeAlleleResult maps gene to list of allele names
type GenomeAlleleResult map[GeneName]AlleleResult

// stringToHash gives a 32-bit hash of a string
func stringToHash(s string) uint32 {
	h := fnv.New32a()
	h.Write([]byte(s))
	return h.Sum32()
}

func logm(level string, msg string, verbose bool) {
	if verbose || level != "DEBUG" {
		fmt.Fprintf(os.Stderr, "%s: %s: %s\n", time.Now().String(), level, msg)
	}
}

func checkResult(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

// joinAlleles takes a list of alleles and returns a comma separated list of them
func joinAlleles(alleles AlleleResult) string {
	keys := make([]string, 0, len(alleles))
	for k := range alleles {
		keys = append(keys, k)
	}
	return strings.Join(keys, ",")
}

// reverse complements a single nucleotide
func reverse(in byte) byte {
	result, ok := REVERSE_MAP[in]
	if !ok {
		log.Fatal("failed to reverse complement")
	}
	return result
}

var REVERSE_MAP = map[byte]byte{'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

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

// addSearchableSequence adds an allele sequence to the hash
func addSearchableSequence(content string, sequence int, db QueryIndex, seedSize int) {
	// for an exact match we only have to hash the start of the query
	if len(content) < seedSize {
		logm("WARN", fmt.Sprintf("sequence %v is length %v, shorter than seed size %v", sequence, len(content), seedSize), false)
		return
	}
	position := 0
	kmer := content[position : position+seedSize]
	kmerHash := stringToHash(kmer)
	entry := QueryPos{name: sequence, pos: position, content: content}
	query, ok := db[kmerHash]
	if !ok {
		query = make([]QueryPos, 0)
	}
	query = append(query, entry)
	db[kmerHash] = query
}

// searchSequence iterates over content and populates result with genes and matching
func searchSequence(content string, db QueryIndex, sequenceNames []string, seedSize int, result GenomeAlleleResult, verbose bool) {
	// populates result with a map from gene name to list of found alleles
	for position := 0; position <= len(content)-seedSize; position += 1 {
		kmer := content[position : position+seedSize] // kmer at curreent position in content
		kmerHash := stringToHash(kmer)
		query, ok := db[kmerHash]
		if ok {
			for _, candidate := range query {
				// check for a match
				if position-candidate.pos >= 0 && position-candidate.pos+len(candidate.content) <= len(content) && content[position-candidate.pos:position-candidate.pos+len(candidate.content)] == candidate.content {
					// it's a match, split the sequence name into gene and allele
					geneAllele := strings.Split(sequenceNames[candidate.name], "_")
					alleles, ok := result[GeneName(geneAllele[0])]
					if !ok {
						alleles = make(AlleleResult, 0)
					}
					alleles[geneAllele[1]] = true
					logm("DEBUG", fmt.Sprintf("%s found at %v", sequenceNames[candidate.name], position-candidate.pos), verbose)
					result[GeneName(geneAllele[0])] = alleles
				}
			}
		}
	}
}

func IndexSequences(sequencesFilename string, sequenceNames *[]string, geneNames map[string]bool, seedSize int, verbose bool) QueryIndex {
	logm("INFO", fmt.Sprintf("processing with seed size %v: %s", seedSize, sequencesFilename), verbose)
	var queryKmers QueryIndex = make(QueryIndex)
	var content *bytes.Buffer

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

	// index sequences
	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if content != nil {
				addSearchableSequence(content.String(), sequenceCount-1, queryKmers, seedSize)
			}
			*sequenceNames = append(*sequenceNames, line[1:])
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
			logm("INFO", fmt.Sprintf("processing: %s: %v lines %v sequences. %v kmers", sequencesFilename, lines, len(*sequenceNames), len(queryKmers)), false)
		}
	}
	addSearchableSequence(content.String(), sequenceCount-1, queryKmers, seedSize)

	logm("INFO", fmt.Sprintf("done processing: %s: %v sequences.", sequencesFilename, len(*sequenceNames)), verbose)
	return queryKmers
}

// FindAlleles finds alleles that match a genome and generates cgst information
func FindAlleles(seedSize int, mismatches int, cgst string, sequencesFilename string, genomes []string, verbose bool) {
	var sequenceNames []string
	var geneNames map[string]bool = make(map[string]bool)
	queryKmers := IndexSequences(sequencesFilename, &sequenceNames, geneNames, seedSize, verbose)

	var cgstGeneNames = make([]string, 0, len(geneNames))
	var cgstIds map[uint32]string = make(map[uint32]string)
	if cgst == "" {
		// no cgst file provided
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
		lines := 0
		for scanner.Scan() {
			if lines == 0 {
				// first line contains the gene names
				fields := strings.Split(scanner.Text(), "\t")
				for _, gene := range fields[1:] {
					cgstGeneNames = append(cgstGeneNames, gene)
				}
			} else {
				fields := strings.SplitN(scanner.Text(), "\t", 2)
				// subsequent lines are of the form id,allele_id1,...
				cgstIds[stringToHash(fields[1])] = fields[0]
			}
			lines++
		}
		logm("INFO", fmt.Sprintf("done processing: %s: %v lines %v sequences.", cgst, lines, len(sequenceNames)), verbose)
	}

	// genotype each genome - list matching sequences
	fmt.Fprintf(os.Stdout, "Filename\tcgST\t%v\n", strings.Join(cgstGeneNames, "\t"))
	var lines int
	var sequenceCount int
	var content *bytes.Buffer

	for _, genomeFilename := range genomes {
		logm("INFO", fmt.Sprintf("processing genome: %s", genomeFilename), verbose)
		var result GenomeAlleleResult = make(GenomeAlleleResult)
		file, err := os.Open(genomeFilename)
		checkResult(err)
		defer file.Close()

		lines = 0
		sequenceCount = 0
		content = new(bytes.Buffer)
		r := bufio.NewReader(file)
		scanner := bufio.NewScanner(r)
		for scanner.Scan() {
			line := scanner.Text()
			if strings.HasPrefix(line, ">") {
				if content != nil {
					searchSequence(content.String(), queryKmers, sequenceNames, seedSize, result, verbose)
					searchSequence(string(reverseComplement(content)), queryKmers, sequenceNames, seedSize, result, verbose)
				}
				sequenceCount++
				content = new(bytes.Buffer)
			} else {
				(*content).WriteString(line)
			}
			lines++
		}
		searchSequence(content.String(), queryKmers, sequenceNames, seedSize, result, verbose)
		searchSequence(string(reverseComplement(content)), queryKmers, sequenceNames, seedSize, result, verbose)
		logm("INFO", fmt.Sprintf("done genome: %s. %v lines. %v sequences.", genomeFilename, lines, sequenceCount), verbose)
		writeResults(genomeFilename, result, geneNames, cgstGeneNames, cgstIds)
	}
}

func writeResults(filename string, result GenomeAlleleResult, geneNames map[string]bool, cgstGeneNames []string, cgstIds map[uint32]string) {
	genomeAlleles := make([]string, 0, len(cgstGeneNames))
	for _, gene := range cgstGeneNames {
		alleles, ok := result[GeneName(gene)]
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
