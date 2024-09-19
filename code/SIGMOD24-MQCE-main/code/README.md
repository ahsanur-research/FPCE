# SIGMOD24-MQCE
Our code, data and additional materials are avaliable here.
## Index  
```shell
.
|- README.md
|- SIGMOD_MQCE_report.pdf
|- Code
|  | - dataset
|  |   |- ca-GrQC
|  |   |- opsahl
|  |   |- condmat
|  |   |- Enron
|  | - FastQC
|  |   |....
```


# DCFastQC
An efficient algorithm for enumerating maximal $\gamma$-quasi-cliques (MQCs).


## Source code info
Programming Language: `C++`
 
Compiler Info: `gcc (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0 ` 

Packages/Libraries Needed: `CMake 2.8`, `makefile`

## Datasets info
All datasets used in our experiments are from the [KONECT](http://konect.cc/networks/ "KONECT") website. We provided four small datasets, namely Ca-GrQc, Opsahl, CondMat, and Enron, used in our paper.



## Setup
### Option 1 (without CMake)
```shell
mkdir build
cd build
g++ ../main.cpp -O3 -o FastQC
```
### Option 2 (using CMake)
```shell
mkdir build
cd build
cmake ..
make
```

## Usage
  ./FastQC {OPTIONS}

    DCFastQC, an efficient algorithm for enumerating MQCs

  OPTIONS:

      -h, --help                          Display this help menu
      -f[dataset], --file=[dataset]       Path to dataset
      -g[para gamma], --g=[para gamma]    The parameter gamma
      -u[para theta] --u=[para theta]     The threshold of number of vertices within a found MQC
      -q[quiet], --q=[quiet]              -q 0 (by default): without output; -q 1: output QCs to file

## Scripts for our experiments
### Compile the code under fold './code/FastQC/'
```shell
mkdir build
cd build
cmake ..
make
```
### Run `DCFastQC` on `Ca-GrQc`:
```shell
./FastQC -f "../../dataset/ca-GrQc.txt" -g 0.9 -u 10
```

### Run `DCFastQC` on `Opsahl`:
```shell
./FastQC -f "../../dataset/opsahl.txt" -g 0.9 -u 20
```

### Run `DCFastQC` on `CondMat`:
```shell
./FastQC -f "../../dataset/condmat.txt" -g 0.9 -u 10
```

### Run `DCFastQC` on `Enron`:
```shell
./FastQC -f "../../dataset/Enron.txt" -g 0.9 -u 23
```

## Running Example

```shell
> mkdir build
> cd build
> cmake ..
> make
>
> ./FastQC -f "../../dataset/ca-GrQc.txt" -g 0.9 -u 10
> # of returned QCs: 1725
> Running Time: 0.015155s
> ./FastQC -f "../../dataset/ca-GrQc.txt" -g 0.9 -u 10 -q 1
> # of returned QCs: 1725
> Running Time: 0.020845s
>
> cat ./output
10 0 22 7 12 13 15 17 2 5 352
10 1 24 8 14 16 3 26 23 4 38
10 6 85 10 11 18 19 25 332 341 353
11 39 65 66 63 40 41 42 64 43 62 78
12 44 79 58 71 72 73 74 75 76 46 54 81
12 45 77 68 48 52 69 70 57 61 67 47 328
12 49 65 66 78 56 59 62 63 64 50 55 81
12 51 74 75 76 79 80 71 72 73 53 60 81
```

### Format of output files
DCFastQC would output all found QCs to a file './code/FastQC/build/output'. For example, in the output file, we have 

    10 0 22 7 12 13 15 17 2 5 352
    10 1 24 8 14 16 3 26 23 4 38
    10 6 85 10 11 18 19 25 332 341 353
    11 39 65 66 63 40 41 42 64 43 62 78
    12 44 79 58 71 72 73 74 75 76 46 54 81
    12 45 77 68 48 52 69 70 57 61 67 47 328
    12 49 65 66 78 56 59 62 63 64 50 55 81
    12 51 74 75 76 79 80 71 72 73 53 60 81
where each line corresponds to a QC. To illusrate, consider the first line "10 0 22 7 12 13 15 17 2 5 352", where the first integer '10' denotes the size of the QC and the rest denotes the vertices that induce a QC, i.e., G[\{0,22,7,12,13,15,17,2,5,352\}].

## Post-processing
Recall that a post-processing procedure is required for filtering out those non-maximal QCs from the output. We put the implementation under './code/maximal_check/'.

### Compile under folder './code/maximal_check/'
```shell
make
```

### Usage
```shell
./maximal graph_file output_file
```

### Runing example
```shell
> ./maximal "../FastQC/build/output" "./output"
num_of_sets: 1725 #maximal cliques: 1665
Maximality checking time: 0.021
```

### Remark
When adopting this implementation for other problems, note that it only accepts the introduced format of output files, e.g.,

    10 0 22 7 12 13 15 17 2 5 352
    10 1 24 8 14 16 3 26 23 4 38 