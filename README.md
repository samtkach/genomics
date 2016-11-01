## Error-Tolerant Short Read Alignment for Genome Sequence Data

usage: `python search_bwt.py [--no-indels] [test|<reference file name>] [<read file name>]`
  * `--no-indels`:  do not allow insertions or deletions in matching
  * `test`:  run default test
  * `<reference file name>`:  reference sequence to match against
  * `<read file name>`:  shotgun reads to align
  
See included sequence files for examples.

This is a proof of concept. It is largely an implementation of string matching using the Burrows-Wheeler transform with a few extra things: FM-indexing, suffix tree search, and some heuristics that "prune" branches at search time. 
