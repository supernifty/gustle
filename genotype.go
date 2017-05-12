package main

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
)

type QueryPos struct {
	name    int    // array index to name of query
	pos     int    // position of kmer in query
	content string // entire query sequence
}

type QueryIndex map[uint32][]QueryPos

// AlleleResult is a list of matching alleles
type AlleleResult map[string]bool

// maps gene to list of matching alleles
type GenomeAlleleResult map[string]AlleleResult

// stringToHash gives a 32-bit hash of a string
func stringToHash(s string) uint32 {
	h := fnv.New32a()
	h.Write([]byte(s))
	return h.Sum32()
}

func joinAlleles(alleles AlleleResult) string {
	keys := make([]string, 0, len(alleles))
	for k := range alleles {
		keys = append(keys, k)
	}
	return strings.Join(keys, ",")
}

// reverse complements a single byte
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

func indexSequences(sequencesFilename string, sequenceNames *[]string, geneNames map[string]bool, seedSize int) QueryIndex {
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
				addSequence(content.String(), sequenceCount-1, queryKmers, seedSize)
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
			logm("INFO", fmt.Sprintf("processing: %s: %v lines %v sequences. %v kmers", sequencesFilename, lines, len(*sequenceNames), len(queryKmers)))
		}
	}
	addSequence(content.String(), sequenceCount-1, queryKmers, seedSize)

	return queryKmers
}

// findAlleles finds alleles that match a genome and generates cgst information
func findAlleles(seedSize int, mismatches int, cgst string, sequencesFilename string, genomes []string) {
	logm("INFO", fmt.Sprintf("processing with seed size %v, up to %v mismatches, cgst: %s: %s", seedSize, mismatches, cgst, sequencesFilename))
	var sequenceNames []string
	var geneNames map[string]bool = make(map[string]bool)
	queryKmers := indexSequences(sequencesFilename, &sequenceNames, geneNames, seedSize)

	logm("INFO", fmt.Sprintf("done processing: %s: %v sequences.", sequencesFilename, len(sequenceNames)))

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
		logm("INFO", fmt.Sprintf("done processing: %s: %v lines %v sequences.", cgst, lines, len(sequenceNames)))
	}

	// genotype each genome - list matching sequences
	fmt.Fprintf(os.Stdout, "Filename\tcgST\t%v\n", strings.Join(cgstGeneNames, "\t"))
	var lines int
	var sequenceCount int
	var content *bytes.Buffer

	for _, genomeFilename := range genomes {
		logm("INFO", fmt.Sprintf("processing genome: %s", genomeFilename))
		var result GenomeAlleleResult = make(GenomeAlleleResult)
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
		writeResults(genomeFilename, result, geneNames, cgstGeneNames, cgstIds)
	}
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
