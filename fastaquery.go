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
	"strings"
	"time"
)

const (
	DEFAULT_SEED_SIZE = 16
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

func addSequence(content string, sequence int, hash map[uint32][]QueryPos, seed_size int) {
	// for an exact match we only have to hash the start of the query
	if len(content) < seed_size {
		fmt.Fprintf(os.Stderr, "%s: WARN: sequence %v is length %v, shorter than seed size %v\n", time.Now().String(), sequence, len(content), seed_size)
		return
	}
	position := 0
	kmer := content[position : position+seed_size]
	kmer_hash := string_to_hash(kmer)
	entry := QueryPos{name: sequence, pos: position, content: content}
	query, ok := hash[kmer_hash]
	if !ok {
		query = make([]QueryPos, 0)
	}
	query = append(query, entry)
	hash[kmer_hash] = query
}

func searchSequence(content string, hash map[uint32][]QueryPos, sequence_name string, genome_filename string, sequence_names []string, seed_size int) {
	for position := 0; position <= len(content)-seed_size; position += 1 {
		kmer := content[position : position+seed_size]
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
	var seed_size int
	var mismatches int
	flag.IntVar(&seed_size, "readlength", 16, "minimum read length of queries (16)")
	flag.IntVar(&mismatches, "mismatches", 0, "mismatches to include (0)")
	flag.Parse()
	if !flag.Parsed() || flag.NArg() < 2 {
		fmt.Fprintf(os.Stderr, "Usage: fastaquery [-readlength min_readlength] [-mismatches mismatches] queries.fq.gz genome.fa [genome.fa...]\n")
		flag.PrintDefaults()
		os.Exit(1)
	}
	sequences_filename := flag.Arg(0)
	fmt.Fprintf(os.Stderr, "%s: processing with seed size %v and up to %v mismatches: %s\n", time.Now().String(), seed_size, mismatches, sequences_filename)

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
				addSequence(content.String(), sequence_count-1, query_kmers, seed_size)
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
	addSequence(content.String(), sequence_count-1, query_kmers, seed_size)
	fmt.Fprintf(os.Stderr, "%s: done processing: %s: %v lines %v sequences.\n", time.Now().String(), sequences_filename, lines, len(sequence_names))

	// process each genome - list matching sequences
	fmt.Fprintf(os.Stdout, "Genome\tGenomeContig\tQueryContig\tGenomePosition\n")
	for _, genome_filename := range flag.Args()[1:] {
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
					searchSequence(content.String(), query_kmers, sequence_name, genome_filename, sequence_names, seed_size)
				}
				sequence_name = line[1:]
				sequence_count++
				content = new(bytes.Buffer)
			} else {
				(*content).WriteString(line)
			}
			lines++
		}
		searchSequence(content.String(), query_kmers, sequence_name, genome_filename, sequence_names, seed_size)
		fmt.Fprintf(os.Stderr, "%s: done genome: %s. %v lines. %v sequences.\n", time.Now().String(), genome_filename, lines, sequence_count)
	}
}
