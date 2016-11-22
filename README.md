hCED: heuristic for Cyclic Edit Distance computation
===

<b>Description</b>: Given two sequences x and y, hCED finds (an approximation of) the cyclic edit 
distance between x and y.

<b>Installation</b>: To compile hCED, please follow the instructions given in file INSTALL.
```
 hCED <options>
 Mandatory arguments:
  -i, --input-file            <str>     (Multi)FASTA input filename.
  -o, --output-file           <str>     Output filename for the rotated sequences.

 Optional for Stage 1 (Algorithm saCSC):
  -q, --q-length              <int>     The q-gram length. Default: 5.
  -l, --block-length          <int>     The length of each block. Default: sqrt(seq_len).
 Optional for Stage 2 (Refinement):
  -P, --refine-blocks         <dbl>     Refine the alignment of saCSC by
                                        checking P blocks of the ends. Default: 1.
  -m, --score-match           <int>     Score of match for refinement. Default: 1.
  -r, --score-mismatch        <int>     Score of mismatch for refinement. Default: -1.
  -f, --score-insertion       <int>     Score of insertion for refinement. Default: -1.
  -g, --score-deletion        <int>     Score of deletion for refinement. Default: -1.
  -O, --gap-open              <int>     Score of gap opening for refinement. Default: NOT USED.
  -E, --gap-extend            <int>     Score of gap extension for refinement. Default: NOT USED.
 Optional for Stage 3 (Edit distance model):
  -e, --edit-distance         <str>     Choose 'Y' to calculate edit distance and 'N' to output
                                        rotation only. Default: Y.
  -S, --cost-substitution     <int>     Cost of substitution. Default: 1.
  -I, --cost-insertion        <int>     Cost of insertion. Default: 1.
  -D, --cost-deletion         <int>     Cost of deletion. Default: 1.
```
<b>License</b>: GNU GPLv3 License; Copyright (C) 2016 Lorraine A. K. Ayad and Solon P. Pissis
