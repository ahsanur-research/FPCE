# FPCE

**SYSTEM REQUIREMENTS:**
* Operating System: Linux (we are using Ubuntu version 20.04)
* Java (we are using openjdk version 11.0.22) – to run ClusterONE
* GCC (we are using GCC 9.4.0) - to compile our codes

**COMPILING:**
* To compile PCE:
cd code/PCE/
make

* To compile FPCE:
cd code/FPCE/
make

* To compile ODES:
cd code/ODES/
make

* To compile FastQC:
cd code/SIGMOD24-MQCE-main
Then follow the instructions in https://github.com/KaiqiangYu/SIGMOD24-MQCE to compile FastQC code, which has been slightly modified to make it work on our inputs 


**SUPPORTED GRAPH FILE FORMATS:**

1. EDGELIST format: 
Each line describes edge whose two endpoints/node-IDs are separated by a space for e.g.
0 1
0 3
1 2
2 3

2. GRH FORMAT:
It assumes that node-IDs of the graph are between 0 to n-1 (where n is the total number of nodes). This format describes an adjacency list where the i-th line (1st line is line-0, 2nd line is line-1, etc.) lists those neighbors of node i that have IDs > i. See details in: https://research.nii.ac.jp/~uno/code/pce.html

For e.g the aforementioned EDGELIST formatted graph file will look like below when converted to GRH format:
1 3
2
3

Typically, graphs are available in networkrepository.org in EDGELIST format. To be able to apply PCE/FPCE on a graph, we need to convert it to GRH format, which can be done via transgrh.pl (see an e.g. given below). Also, to be able to apply transgrh.pl on an EDGELIST formatted file, the node-IDs must be between 0 to n-1; otherwise, we need to convert the EDGELIST file to another EDGELIST file using transnum.pl (see details in: https://research.nii.ac.jp/~uno/code/pce.html).


**INSTRUCTIONS TO REPLICATE OUR EXPERIMENTS:**

__1. REAL and SYNTHETIC GRAPH ANALYSIS__

* Run: cd real-synth-graph-analysis

* To generate synthetic graphs, run jupyter notebook (in bash run: jupyter-notebook), and then run generate_synth_graphs.ipynb there. However, running it will replace all the synthetic graphs generated+used for our paper.


* To run FPCE on a graph:
** preprocess (remove self-loop and multi-edges) an edgelist formatted graph, for e.g.:
python ../Graph_check.py --i real-graphs/bio-grid-human.txt --o real-graphs/bio-grid-human.pel

For your convenience, we have preprocessed all real graphs and put them inside real-synth-graph-analysis/real-graphs directory.

** Convert it to a grh file:
../transgrh.pl < real-graphs/bio-grid-human.pel > real-graphs/bio-grid-human.grh

** Finally to compute all maximal (10,0.9)-pseudo-cliques from this graph using FPCE as well as to report its running time & memory consumption (RSS), run:
/usr/bin/time -o log/time-fpce-biogrid-10-0.9.txt -f "RSS=%M TIME=%S+%U" ../code/FPCE/fpce M -l 10 real-graphs/bio-grid-human.grh 0.9 outputs/biogrid-10-0.9.out > log/log-fpce-biogrid-10-0.9.txt

The output file (here biogrid-10-0.9.out) simply contains a list of maximal (10,0.9)-pseudo-cliques, one per line. Each line is a list of nodes in a pseudo-clique separated by spaces.

To run FPCE without generating the output file, simply omit it from the command line:
/usr/bin/time -o log/time-fpce-biogrid-10-0.9.txt -f "RSS=%M TIME=%S+%U" ../code/FPCE/fpce M -l 10 real-graphs/bio-grid-human.grh 0.9 > log/log-fpce-biogrid-10-0.9.txt

To run FPCE without Edge Bound code, run the command line:
/usr/bin/time -o time-fpce-biogrid-10-0.9.txt -f "RSS=%M TIME=%S+%U" ../code/FPCE/fpce_wo_eb M -l 10 real-graphs/bio-grid-human.grh 0.9 > log/log-fpce-biogrid-10-0.9.txt


* To run PCE on a graph:
Preprocess and prepare a grh file, like before. Then run PCE in the following way:
/usr/bin/time -o log/time-pce-biogrid-10-0.9.txt -f "RSS=%M TIME=%S+%U" ../code/PCE/pce M -l 10 real-graphs/bio-grid-human.grh 0.9 > log/log-pce-biogrid-10-0.9.txt

* To run ODES on a graph:
** Convert the previously generated grh file to ODES compatible format:
../code/FPCE/grh2odes real-graphs/bio-grid-human.grh > real-graphs/bio-grid-human.odes

** Run ODES on that graph file:
/usr/bin/time -o log/time-odes-biogrid-10-0.9.txt -f "RSS=%M TIME=%S+%U" ../code/ODES/denseSG real-graphs/bio-grid-human.odes 0.9 10 mx > outputs/biogrid-10-0.9-odes.out

__2. E. COLI DATA ANALYSIS:__

* Run: cd e-coli-analysis
* Copy the contents of https://github.com/tanghaibao/goatools into a folder called goatools-main
* run the commands mentioned in process_ecoli_data.sh
* run jupyter notebook: jupyter-notebook
* run evaluate_ecoli_pseudocliques.ipynb in jupyter notebook

__3. Other tasks:__

* To prepare tables describing information of graphs (Tables 2,3,4 in our paper) you can use the following steps to gather info of all graphs
cd real-synth-graph-analysis
Run collect_graph_info.ipynb using jupyter-notebook

* To prepare list of commands to run PCE/FPCE/ODES:
cd real-synth-graph-analysis
PATH=$PATH:../code/PCE:../code/FPCE:../code/ODES
Then run command_generator.ipynb in jupyter-notebook

