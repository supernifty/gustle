package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"fmt"
	"hash/fnv"
	"log"
	"os"
	"strings"
	"time"
)

const (
	SEED_SIZE = 16
)

func check_result(e error) {
	if e != nil {
		log.Fatal(e)
	}
}

type QueryPos struct {
	name    int
	pos     int
	content string
}

func string_to_hash(s string) uint32 {
	h := fnv.New32a()
	h.Write([]byte(s))
	return h.Sum32()
}

func addSequence(content string, sequence int, hash map[uint32][]QueryPos) {
	// for an exact match we only have to hash the start of the query
	position := 0
	kmer := content[position : position+SEED_SIZE]
	kmer_hash := string_to_hash(kmer)
	entry := QueryPos{name: sequence, pos: position, content: content}
	query, ok := hash[kmer_hash]
	if !ok {
		query = make([]QueryPos, 0)
	}
	query = append(query, entry)
	hash[kmer_hash] = query
}

func searchSequence(content string, hash map[uint32][]QueryPos, sequence_name string, genome_filename string, sequence_names []string) {
	for position := 0; position <= len(content)-SEED_SIZE; position += 1 {
		kmer := content[position : position+SEED_SIZE]
		kmer_hash := string_to_hash(kmer)
		query, ok := hash[kmer_hash]
		if ok {
			for _, candidate := range query {
				if position+len(candidate.content) <= len(content) && content[position:position+len(candidate.content)] == candidate.content {
					fmt.Fprintf(os.Stdout, "%s\t%s\t%s\t%v\n", genome_filename, sequence_name, sequence_names[candidate.name], position)
				}
			}
		}
	}
}

func main() {
	sequences_filename := os.Args[1]
	fmt.Fprintf(os.Stderr, "%s: processing: %s\n", time.Now().String(), sequences_filename)

	file, err := os.Open(sequences_filename)
	check_result(err)
	defer file.Close()

	gr, err := gzip.NewReader(file)
	check_result(err)
	defer gr.Close()

	scanner := bufio.NewScanner(gr)
	lines := 0
	sequence_count := 0

	var sequence_names []string
	var query_kmers map[uint32][]QueryPos = make(map[uint32][]QueryPos)
	var content *bytes.Buffer

	for scanner.Scan() {
		line := scanner.Text()
		if strings.HasPrefix(line, ">") {
			if content != nil {
				addSequence(content.String(), sequence_count-1, query_kmers)
			}
			sequence_names = append(sequence_names, line[1:])
			sequence_count++
			content = new(bytes.Buffer)
		} else {
			(*content).WriteString(line)
		}
		lines++
		if lines%1000000 == 0 {
			fmt.Fprintf(os.Stderr, "%s: processing: %s: %v lines %v sequences. %v kmers\n", time.Now().String(), sequences_filename, lines, len(sequence_names), len(query_kmers))
			// debug.FreeOSMemory()
		}
	}
	addSequence(content.String(), sequence_count-1, query_kmers)
	fmt.Fprintf(os.Stderr, "%s: done processing: %s: %v lines %v sequences.\n", time.Now().String(), sequences_filename, lines, len(sequence_names))

	// process each genome - list matching sequences
	fmt.Fprintf(os.Stdout, "Genome\tGenomeContig\tQueryContig\tGenomePosition\n")
	for _, genome_filename := range os.Args[2:] {
		fmt.Fprintf(os.Stderr, "%s: processing genome: %s\n", time.Now().String(), genome_filename)
		file, err := os.Open(genome_filename)
		check_result(err)
		defer file.Close()

		lines = 0
		sequence_count = 0
		content = new(bytes.Buffer)
		r := bufio.NewReader(file)
		scanner := bufio.NewScanner(r)
		var sequence_name string
		for scanner.Scan() {
			line := scanner.Text()
			if strings.HasPrefix(line, ">") {
				if content != nil {
					searchSequence(content.String(), query_kmers, sequence_name, genome_filename, sequence_names)
				}
				sequence_name = line[1:]
				sequence_count++
				content = new(bytes.Buffer)
			} else {
				(*content).WriteString(line)
			}
			lines++
		}
		searchSequence(content.String(), query_kmers, sequence_name, genome_filename, sequence_names)
		fmt.Fprintf(os.Stderr, "%s: done genome: %s. %v lines. %v sequences.\n", time.Now().String(), genome_filename, lines, sequence_count)
	}
}
