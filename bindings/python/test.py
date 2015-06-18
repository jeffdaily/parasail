import parasail

result = parasail.sw_table("asdf","asdfasdf",10,1,parasail.blosum62)
print result
print result.saturated
print result.score
print result.matches
print result.similar
print result.length
print result.score_table

result = parasail.sw_stats_table("asdf","asdfasdf",10,1,parasail.blosum62)
print result
print result.saturated
print result.score
print result.matches
print result.similar
print result.length
print result.score_table
print result.matches_table
print result.similar_table
print result.length_table

result = parasail.sw_rowcol("asdf","asdfasdf",10,1,parasail.blosum62)
print result
print result.saturated
print result.score
print result.matches
print result.similar
print result.length
print result.score_row
print result.score_col

result = parasail.sw_stats_rowcol("asdf","asdfasdf",10,1,parasail.blosum62)
print result
print result.saturated
print result.score
print result.matches
print result.similar
print result.length
print result.score_row
print result.score_col
print result.matches_row
print result.matches_col
print result.similar_row
print result.similar_col
print result.length_row
print result.length_col

result = parasail.sw("asdf","asdf",10,1,parasail.blosum62)
print result
print result.saturated
print result.score
print result.matches
print result.similar
print result.length
print result.score_table

result = parasail.sw_stats("asdf","asdf",10,1,parasail.blosum62)
print result
print result.saturated
print result.score
print result.matches
print result.similar
print result.length

result = parasail.sw_scan_32("asdf","asdf",10,1,parasail.blosum62)
print result
print result.saturated
print result.score
print result.matches
print result.similar
print result.length

result = parasail.sw_scan_16("asdf","asdf",10,1,parasail.blosum62)
print result
print result.saturated
print result.score
print result.matches
print result.similar
print result.length

result = parasail.sw_stats_striped_16("asdf","asdf",10,1,parasail.blosum62)
print result
print result.saturated
print result.score
print result.matches
print result.similar
print result.length

print parasail.blosum62.name
print parasail.blosum62.size
print parasail.blosum62.matrix
