# `kelp` - contig tiling to reference sequence simulator

`kelp` generates "small" reference sequences given contig tilings.

## Heuristic used in `kelp` to generate reference sequences

Given reference seqences represented as sequences of contigs.
Kelp heuristically finds compatible sequences for said contigs by:
1. Iterating over references inserting (k-1) overlaps to yield 
terminating (k-1) bases for contigs
2. Iterating over contigs, nucleotides are then inserted for
each contig so that no k-mer appears more than once.

Simulated references are "small" since contigs need only valid
overlaps with no repeated k-mers. Compatible contigs can thus (usually)
be found with less than `(k-1) + k + (k-1)` bases.