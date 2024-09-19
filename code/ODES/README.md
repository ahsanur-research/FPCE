Downloaded from http://dense.sourceforge.net/ \
The corresponding paper is : http://bioinformatics.oxfordjournals.org/content/26/21/2788.abstract

1) if desired, uncomment "# define CONSISTENCY_CHECK"
   at line 133 of utils.c for the ReadGraph() routine
   (slow for large graphs), and edit "#define NUM_THREADS"
   at line 6 of odes.h to reflect your architecture
2) make
3) time ./denseSG [graph file] [min density] [min nodes] mxcsvstat \
   or, \
   time ./denseSG [graph file] [min density] [min nodes] fscsvstat


test graphs are in the "../../data/odes/test-inputs/" directory

If you want to experiment with excluding some edges, uncomment 
"#define EXCLUDED_EDGES" in odes.h and define some excluded edges
in denseSG.c

e.g. run: to get order >= 4 maximal (indicated by "mx" of the argument "mxcsvstat") dense subgraphs which have density >= 0.7
```
time ./denseSG ../../data/odes/test-inputs/sample1.txt 0.7 4 mxcsvstat > result-sample1-0.7-4-mxcsvstat.txt
```
to get dense subgraphs of order = 4 (fixed-ordered indicated by "fs" of the argument "fscsvstat") which have density >= 0.7
```
time ./denseSG ../../data/odes/test-inputs/sample1.txt 0.7 4 fscsvstat > result-sample1-0.7-4-fscsvstat.txt
```
Input file format (almost like .grh format used by the PCE paper, see "sample1.txt" file 
which is a 7 node graph whose density = 0.52): \
First line contains a single number: number of nodes \
Each subsequent line contains node_label\tspace-separated-list-of-nodes-adjacent-to-it \
e.g. input file
```
7
0	1 2
1	0 2 6
2	0 1 6
3	4 5 6
4	3 5 6
5	3 4 6
6	1 2 3 4 5
```

Unlike .grh format (used by PCE and QCE), ODES requires that both (a,b) and (b,a) edges are to be specified; 
for e.g. in the above graph both (0,1) and (1,0) were mentioned as an edge.

Output file is in csv format as indicated by "csvstat" option (tab-separated) which contains the following information for each dense subgraph:
* number of nodes
* density
* space separated list of nodes in the dense subgraph

E.g. of an output file (that we get from the above input file):

For mxcsvstat option, to get order >= 4 maximal dense subgraphs which have density >= 0.7
```
4	0.833333	0 1 2 6 
5	0.700000	1 3 4 5 6 
5	0.700000	2 3 4 5 6
```

For fscsvstat option, to get dense subgraphs of order = 4 which have density >= 0.7
```
4	0.833333	0 1 2 6 
4	1.000000	3 4 5 6
```

