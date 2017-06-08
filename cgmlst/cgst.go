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
	Ids map[uint32]string
}

func CreateCGST(filename string, verbose bool) CGST {
	var result CGST
	result.GeneNames = make([]string, 0)
	result.Ids = make(map[uint32]string)
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
			// subsequent lines are of the form id,allele_id1,...
			result.Ids[stringToHash(fields[1])] = fields[0]
		}
		lines++
	}
	logm("INFO", fmt.Sprintf("processing '%s': done reading %v lines", filename, lines), verbose)
	return result
}
