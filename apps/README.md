# parasail applications

Though the parasail library is the focus of this work, we do provide some example applications which may be of immediate use.

## parasail_aligner

The aligner tool will take a FASTA-formatted query file as input and align all of the query sequences against a set of database sequences found in a second FASTA-formatted database file.  Alternatively, if only one file is supplied, all of the sequences in the file will be compared against themselves.  The alignment routine to use defaults to one of the Smith-Waterman routines, but any of the parasail routines may be selected using the appropriate command-line parameter.

### Command-Line Interface

```bash
usage: parasail_aligner [-a funcname] [-c cutoff] [-x] [-e gap_extend] [-o gap_open] [-m matrix] -f file [-q query_file] [-g output_file] 

Defaults:
   funcname: sw_stats_striped_16
     cutoff: 7, must be >= 1, exact match length cutoff
         -x: if present, don't use suffix array filter
 gap_extend: 1, must be >= 0
   gap_open: 10, must be >= 0
     matrix: blosum62
       file: no default, must be in FASTA format
 query_file: no default, optional, must be in FASTA format
output_file: parasail.csv
```

### Using the Enhanced Suffix Array Filter

One feature of this tool is its ability to filter out sequence pairs based on an exact-match cutoff.  Using the cutoff paramter (-c), the filter will keep only those pairs of sequences which contain an exact match of length greater than or equal to the cutoff.  The assumption is that any pair of sequences which are highly similar should also contain an exact-matching k-mer of at length at least c, our cutoff.  This is similar to the seed and extend model of sequence alignment found in other tools, however, our filter allows for arbitrarily long exact-matches (aka seeds) and once a match is found the entire alignment is performed rather than extending the seed.  The filter is turned on by default but can be disabled using the -x command-line parameter.

### Output

The output is always comma-separated values (CSV).  However, the exact output will depend on how the tool is used.  At minimum, there will be three values per line in the file.

index1, index2, score[, matches, similarities, length]

If a query was used, then index1 is the query index and index2 is the database index.

If a query was not used, then both index1 and index2 refer to the same single FASTA file that was supplied.

If a statistics-calculating function is used, for example 'sw_stats_striped_16', then the number of exact matches, similarities, and alignment length are also computed and returned.

If a score-only function is used, for example 'sw_striped_16', then only the score is returned.

## parasail_all

The parasail_all tool will take a FASTA-formatted set of sequences and all of the sequences in the file will be compared against themselves.  The alignment routine defaults to a statistics-calculating semi-global (Needleman-Wunsch without end gap penalties), but any of the statistics-calculating parasail routines may be selected using the appropriate command-line parameter.  The reason statistics must be calculated is that the output depends on them.  This application is used in a metagenomics workflow, creating a homology graph as output which is later processed by a community detection application.  The 'edges' in the graph consist of any highly similar pair of sequences such that their alignment meets certain criteria.

 * AOL = percent alignment length -- the alignment must cover at least XX percent of the longer sequence.
 * SIM = percent exact matches -- the alignment must contain at least XX percent exact character matches.
 * OS = percent optimal score -- the calculated score must be XX percent of the longer sequence's self score.

### Command-Line Interface

```bash
usage: parasail_all [-a funcname] [-c cutoff] [-e gap_extend] [-o gap_open] [-m matrix] [-l AOL] [-s SIM] [-i OS] -f file [-g output] 

Defaults:
  funcname: sg_stats_striped_16
    cutoff: 7, must be >= 1, exact match length cutoff
gap_extend: 1, must be >= 0
  gap_open: 10, must be >= 0
    matrix: blosum62
       AOL: 80, must be 0 <= AOL <= 100, percent alignment length
       SIM: 40, must be 0 <= SIM <= 100, percent exact matches
        OS: 30, must be 0 <= OS <= 100, percent optimal score over self score
      file: no default, must be in FASTA format
    output: edges.csv
```

### Output

The output is always comma-separated values (CSV).  An 'edge' is only output if it meets the criteria described above.

index1, index2, length/max_lengh, matches/length, score/self_score

## parasail_query

The parasail_query application is nearly identical to parasail_all except it takes an optional query file as input.  The output then contains only pairs of sequences which connect both input FASTA files; it will not report pairs of sequences where both sequences come from the same file.  This is useful if you have analyzed two files and now want to analyze the combined file without duplicating work.
