package genotype

import (
	"bufio"
	"fmt"
	"os"
	"strings"
)

// CGST maintains a list of gene names and a mapping from a list of alleles to an id
type CGST struct {
	GeneNames []string
	Ids map[string]string // alleles and cgstid
}

// CreateCGST creates a new instance of CGST
func CreateCGST(filename string, verbose bool) CGST {
	var result CGST
	result.GeneNames = make([]string, 0)
	result.Ids = make(map[string]string)
	file, err := os.Open(filename)
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
				result.GeneNames = append(result.GeneNames, gene)
			}
		} else {
			fields := strings.SplitN(scanner.Text(), "\t", 2)
			// subsequent lines are of the form id<tab>allele_id1<tab>...
			result.Ids[fields[1]] = fields[0]
		}
		lines++
	}
	logm("INFO", fmt.Sprintf("processing '%s': done reading %v lines", filename, lines), verbose)
	return result
}

// matchesCGST returns number of matching alleles (e.g. 0,16 0) against this cgst entry (e.g. 0 0)
func matchesCGST(alleles []string, cgst []string) int {
	// look at each allele
	matches := 0
	for idx, allele := range alleles {
		// look at each found allele
		for _, possible := range strings.Split(allele, ",") {
			if possible == cgst[idx] {
				matches += 1
				break
			}
		}
	}
	return matches
}

// FindFirst returns the first matching CGST entry 
func (cgst CGST) FindFirst(alleles []string) (string, bool) {
	// check multiple options
  for key, value := range cgst.Ids {
		if matchesCGST(alleles, strings.Split(key, "\t")) == len(alleles) {
			return value, true
		}
	}
	return "", false
}

// FindBest returns the best matching cgst entry, as well as the number of matches
func (cgst CGST) FindBest(alleles []string) (string, int) {
	// check multiple options
	best := ""
	matches := -1
  for key, value := range cgst.Ids {
		candidate := matchesCGST(alleles, strings.Split(key, "\t"))
		if candidate > matches {
			matches = candidate
			best = value
		}
	}
	return best, matches
}
