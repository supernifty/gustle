package main

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"flag"
	"fmt"
	"os"
	"sort"
	"strings"
)

type QueryIndex map[uint32][]QueryPos

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

// genotype finds alleles that match a genome and generates cgst information
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

	for _, genomeFilename := range genotypeCommand.Args()[1:] {
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
