import parasail

result = parasail.sw("asdf","asdf",10,1,parasail.blosum62)
print result
print result.saturated
print result.score
print result.matches
print result.similar
print result.length

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
