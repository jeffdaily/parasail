## Comparison with other aligners

**Special thanks to Martin Šošić for the original version of this
document in the OpAl software package. See
https://github.com/Martinsos/opal/blob/master/aligner_comparison.md for
the original. **

Here are described results of speed comparison of parasail with other
aligners: SSW, OpAl, and SWIPE.  SSW and SWIPE do only Smith-Waterman
alignment so we compared only for SW.

Aligners were tested by quering sequences against UniProtKB/Swiss-Prot
database (contains 547964 sequences).  Database can be obtained from
www.uniprot.org/downloads -> UniProtKB/Swiss-Prot.  

Specific sequence can also be obtained from www.uniprot.org by searching
it by name (Search tab).

All aligners were tested with following parameters:
* number of threads = 1
* gap opening = 3
* gap extension = 1
* score matrix = BLOSUM50

Only scores were calculated (not alignments). Time spent to read
sequences and database was not measured. In the case of SSW, output was
programmatically suppressed to avoid its effect on alignment timing.

How aligners were called:
* SSW: `./ssw_test -p uniprot_sprot.fasta <query_file>`
* OpAl: `./opal_aligner -s <query_file> uniprot_sprot.fasta`
* SWIPE: `./swipe -a 1 -p 1 -G 3 -E 1 -M BLOSUM50 -b 0 -i <query_file> -d uniprot_sprot` NOTE: database had to be preprocessed for SWIPE using _makeblastdb_
* parasail(16-bit): `parasail_aligner -x -t 1 -a sw_striped_profile_16 -o 3 -e 1 -m blosum50 -f uniprot_sprot -q <query_file>`
* parasail(8-bit): `parasail_aligner -x -t 1 -a sw_striped_profile_8 -o 3 -e 1 -m blosum50 -f uniprot_sprot -q <query_file>`

Following table shows how much time took for different sequences to be
aligned against UniProtKB/Swiss-Prot database.

All times are in seconds. Test were performed on a MacBook Pro i5 CPU @ 2.53GHz
with 8GB RAM (SSE41 support). The compiler was Apple LLVM version 6.0
(clang-600.0.57. The times are an average of three runs.

|                               |O74807  |P19930  |Q3ZAI3  |P18080|
|-------------------------------|--------|--------|--------|------|
| **query length**              |110     |195     |390     |513   |
| **SSW**                       |14.9    |22.7    |44.4    |54.4  |
| **OpAl(SSE4.1)**              |15.2    |21.5    |35.9    |44.6  |
| **SWIPE**                     |7.68    |13.3    |24.7    |32.0  |
| **_parasail(SSE4.1) 16-bit_** |10.7    |15.2    |24.8    |30.5  |
| **_parasail(SSE4.1) 8-bit_**  |9.9     |13.7    |20.0    |23.7  |
