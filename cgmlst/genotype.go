package genotype

//
// this package implements
//

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/gob"
	"fmt"
	"hash/fnv"
	"log"
	"os"
	"strings"
	"time"
)

const (
	NO_CGST = "NA"
)

// QueryList contains indexed queries and the names of all sequences
type QueryList struct {
	Index QueryIndex
	Names []string // list of query fasta names
	SeedSize int
	Cgst CGST
}

// QueryIndex maps hashes of kmers to all possible source locations
type QueryIndex map[uint32][]QueryPos

// QueryPos provides an index into a sequence name, and a position
type QueryPos struct {
	Name    int    // array index to name of query
	Pos     int    // position of kmer in query
	Content string // entire query sequence
}

// GeneName names a short sequence
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
func addSearchableSequence(content string, sequence int, db QueryList) {
	// for an exact match we only have to hash the start of the query
	if len(content) < db.SeedSize {
		logm("WARN", fmt.Sprintf("sequence %v is length %v, shorter than seed size %v", sequence, len(content), db.SeedSize), false)
		return
	}
	position := 0
	kmer := content[position : position+db.SeedSize]
	kmerHash := stringToHash(kmer)
	entry := QueryPos{Name: sequence, Pos: position, Content: content}
	query, ok := db.Index[kmerHash]
	if !ok {
		query = make([]QueryPos, 0)
	}
	query = append(query, entry)
	db.Index[kmerHash] = query
}

// searchSequence iterates over content and populates result with genes and matching
func searchSequence(content string, db QueryList, result GenomeAlleleResult, verbose bool, reverseComplement bool) {
	// populates result with a map from gene name to list of found alleles
	for position := 0; position <= len(content)-db.SeedSize; position += 1 {
		kmer := content[position : position+db.SeedSize] // kmer at curreent position in content
		kmerHash := stringToHash(kmer)
		query, ok := db.Index[kmerHash]
		if ok {
			// logm("DEBUG", fmt.Sprintf("found %v potential locations for %s", len(query), content), verbose)
			for _, candidate := range query {
				// check for a match
				if position-candidate.Pos >= 0 && position-candidate.Pos+len(candidate.Content) <= len(content) && content[position-candidate.Pos:position-candidate.Pos+len(candidate.Content)] == candidate.Content {
					// it's a match, split the sequence name into gene and allele
					geneAllele := strings.Split(db.Names[candidate.Name], "_")
					alleles, ok := result[GeneName(geneAllele[0])]
					if !ok {
						alleles = make(AlleleResult, 0)
					}
					alleles[geneAllele[1]] = true
					if reverseComplement {
						logm("DEBUG", fmt.Sprintf("%s found at reverse complement -%v (%v)", db.Names[candidate.Name], position-candidate.Pos, len(content)-len(candidate.Content)-position+candidate.Pos), verbose)
					} else {
						logm("DEBUG", fmt.Sprintf("%s found at %v", db.Names[candidate.Name], position-candidate.Pos), verbose)
					}
					result[GeneName(geneAllele[0])] = alleles
				} else {
					// logm("DEBUG", fmt.Sprintf("didn't match %s", candidate.Content), verbose)
				}
			}
		} else {
			// logm("DEBUG", fmt.Sprintf("didn't find hash for %s", content), verbose)
		}
	}
}

// IndexSequences generates an index from a gzipped list of alleles for genotyping
func IndexSequences(sequences []string, cgstFilename string, seedSize int, verbose bool) QueryList {
	logm("INFO", fmt.Sprintf("processing with cgst file '%s', seed size %v: %v sequence file(s)", cgstFilename, seedSize, len(sequences)), verbose)
	var queryList QueryList
	queryList.SeedSize = seedSize
	queryList.Index = make(QueryIndex)
	queryList.Cgst = CreateCGST(cgstFilename, verbose)
	var content *bytes.Buffer

	lines := 0
	sequenceCount := 0

	for _, sequenceFilename := range sequences {
		logm("INFO", fmt.Sprintf("processing '%s'...", sequenceFilename), verbose)
		// open sequences (either .fa or .fa.gz) file for reading
		file, err := os.Open(sequenceFilename)
		checkResult(err)

		var scanner *bufio.Scanner
		if strings.HasSuffix(sequenceFilename, ".gz") {
			gr, err := gzip.NewReader(file)
			checkResult(err)
			scanner = bufio.NewScanner(gr)
		} else {
			scanner = bufio.NewScanner(file)
		}


		// index sequences
		for scanner.Scan() {
			line := scanner.Text()
			if strings.HasPrefix(line, ">") {
				if content != nil {
					addSearchableSequence(content.String(), sequenceCount-1, queryList)
				}
				queryList.Names = append(queryList.Names, line[1:])
				sequenceCount++
				content = new(bytes.Buffer)
			} else {
				(*content).WriteString(line)
			}
			lines++
			if lines%1000000 == 0 {
				logm("INFO", fmt.Sprintf("processing %s: %v lines %v sequences. %v kmers", sequenceFilename, lines, len(queryList.Names), len(queryList.Index)), false)
			}
		}
		addSearchableSequence(content.String(), sequenceCount-1, queryList)
		logm("INFO", fmt.Sprintf("processing '%s': done", sequenceFilename), verbose)

		file.Close()
	}

	logm("INFO", fmt.Sprintf("processing %v file(s): done. %v sequences", len(sequences), len(queryList.Names)), verbose)
	return queryList
}

// SaveIndex writes an indexed collection of indexes for use with the genotype command
func SaveIndex(target string, source QueryList, verbose bool) {
	logm("INFO", fmt.Sprintf("saving index to %s...", target), verbose)
	file, err := os.Create(target)
	checkResult(err)
	defer file.Close()

	gr := gzip.NewWriter(file)
	defer gr.Close()

	encoder := gob.NewEncoder(gr)

	err = encoder.Encode(source.Names)
	checkResult(err)
	logm("INFO", fmt.Sprintf("%v sequence names saved", len(source.Names)), verbose)

	err = encoder.Encode(source.SeedSize)
	checkResult(err)

	err = encoder.Encode(source.Cgst)
	checkResult(err)

	// save the index, but go has a size limit
	indexSize := len(source.Index)
	err = encoder.Encode(indexSize)
	checkResult(err)
	logm("INFO", fmt.Sprintf("%v queries to save...", indexSize), verbose)

	count := 0
	for key, value := range source.Index {
		err = encoder.Encode(key)
		checkResult(err)
		err = encoder.Encode(value)
		checkResult(err)
		count++
		if count%10000 == 0 {
			logm("INFO", fmt.Sprintf("processing: saved %v items", count), false)
		}
	}

	logm("INFO", fmt.Sprintf("saving index to %s: done", target), verbose)
}

func LoadIndex(source string, verbose bool) QueryList {
	logm("INFO", fmt.Sprintf("loading index from %s...", source), verbose)

	var result QueryList

	// open fa.gz file
	file, err := os.Open(source)
	checkResult(err)
	defer file.Close()

	gr, err := gzip.NewReader(file)
	checkResult(err)
	defer gr.Close()

	decoder := gob.NewDecoder(gr)

	err = decoder.Decode(&result.Names)
	checkResult(err)
	logm("INFO", fmt.Sprintf("%v sequence names restored", len(result.Names)), verbose)

	err = decoder.Decode(&result.SeedSize)
	checkResult(err)

	err = decoder.Decode(&result.Cgst)
	checkResult(err)

	var indexSize int
	err = decoder.Decode(&indexSize)
	checkResult(err)

	result.Index = make(QueryIndex)

	count := 0
	for i := 0; i < indexSize; i++ {
		var key uint32
		var val []QueryPos
		err = decoder.Decode(&key)
		checkResult(err)
		err = decoder.Decode(&val)
		checkResult(err)
		result.Index[key] = val
		count++
		if count%1000 == 0 {
			logm("INFO", fmt.Sprintf("processing: loaded %v items", count), false)
		}
		// logm("DEBUG", fmt.Sprintf("last key: %v, values: %v", key, len(val)), verbose)
	}

	logm("INFO", fmt.Sprintf("loading index from %s - loaded %v: done", source, len(result.Index)), verbose)

	return result
}

// FindAlleles finds alleles that match a genome and generates cgst information
func FindAlleles(db QueryList, mismatches int, genomes []string, verbose bool) {
	logm("INFO", "find alleles...", verbose)

	// genotype each genome - list matching sequences
	fmt.Fprintf(os.Stdout, "Filename\tcgST\t%v\n", strings.Join(db.Cgst.GeneNames, "\t"))
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
					searchSequence(content.String(), db, result, verbose, false)
					searchSequence(string(reverseComplement(content)), db, result, verbose, true)
				}
				sequenceCount++
				content = new(bytes.Buffer)
			} else {
				(*content).WriteString(line)
			}
			lines++
		}
		searchSequence(content.String(), db, result, verbose, false)
		searchSequence(string(reverseComplement(content)), db, result, verbose, true)
		logm("INFO", fmt.Sprintf("done genome: %s. %v lines. %v sequences.", genomeFilename, lines, sequenceCount), verbose)
		writeResults(genomeFilename, result, db.Cgst.GeneNames, db.Cgst.Ids)
	}
	logm("INFO", "find alleles: done", verbose)
}

func writeResults(filename string, result GenomeAlleleResult, cgstGeneNames []string, cgstIds map[uint32]string) {
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
	// TODO need more than just a perfect match
	cgstId, ok := cgstIds[stringToHash(alleles)]
	if !ok {
		cgstId = NO_CGST
	}
	fmt.Fprintf(os.Stdout, "%s\t%s\t%s\n", filename, cgstId, alleles) // filename, cgST, alleles
}
