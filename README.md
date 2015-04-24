SLIDESORT-BPR
==========================
Copyright &copy; 2013 Kana Shimizu. All rights reserved.

Slidesort-BPR (Break Point Reads) detector is a reference-free method for
detecting clusters of breakpoints reads from the chromosomal rearrangements.

Input format
------------

Multi-fasta format is acceptable.
Break Point Reads Detector accepts both DNA sequence and protein sequence.

```
> seq 1
ATGCTAGCTGATACATCTAGCTCGTACGTACGTCAGTCGTAGT
CTGACTGACTAGCTAGCTAGCATCGTACGTACGTCGTAGCTAC
> seq 2
ATGCTAGCTGATACATCTAGCTCGTACGTACGTCAGTCGTAGT
CTGACTGACTAGCTAGCTAGCATCGTACGTACGTCGTAGCTAC
```

Execution
------------

Initially, the user is required to set the path by executing this line:

```
export LD_LIBRARY_PATH=. 
```

Typical usage is as follows:

```
./bpr-detector -d <distance> -IC <input_control_fasta_file> -IT <input_target_fasta_file> -sl <size_of_splitted_sequence> -b <value of TAU> -o <output_file> -M I (-M:default set to 'I')
```


Sequences can include unknown characters.
In default setting,  Slidesort-BPR excludes sequences with unknown characters such as `N`, `X`.
To include sequences with unknown characters, use `-u` option.

Options
-------

Basic  options:
```
-d  distance threshold
-IT input target file name
-IC input control file name
-o  output filename
-sl size of splitted sequence
-M  method of caluclating M-value I: calculate from the mean degree of target, T: calculate from the table(default=I)
-MD Depth(when you set -M T, It is necessary to specify this value)
-ME Error Rate(when you set -M T, It is necessary to specify this value)
-b  parameter of TAU
-t  distance type  E: edit-distance  H: hamming-distance (default=E)
-v  search with both original seq and reverse complement seq.
```
Advanced options:
```
-c  type of input string.
    DNA: DNA seq, PROTEIN: protein seq, INT: integer seq (default=DNA)
-g  gap extention cost (default=1, must be positive value. it is better to use larger value to avoid slow-down of the search.)
-G  gap open cost (default=0, must be positive value)
-k  size of sorting key
-u  do not exclude sequences with unknown character. ex) n, N, Z, etc...
-V  output a same pair twice if dist(A,B)<=d and dist(A,B')<=d. B' is reverse complement of B. (use with -v)
-mt number of threads (multi-threading mode)
-mp number of tasks (multi-threading mode)
-mr ratio of tasks(<number of tasks> = <number of threads> * <this value>) (multi-threading mode)
```
