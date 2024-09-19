/*
 odes.h                                                  sudowrestler@gmail.com
 
 A breadth-first algorithm to find dense subgraphs with density >= 1/2
 Copyright (C) 2010 James Long
 
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define ODES_H
#define NUM_THREADS   1
//#define EXCLUDED_EDGES
/*
 EXCLUDED_EDGES is experimental:
 This exact algorithm can be combined with a faster heuristic routine by passing
 it a NULL terminated list of edges found by the heuristic routine. This list 
 will be excluded from the list of edges that initialize the algorithm.
 
 The format of each list entry is "label<whitespace>label". 
 
 Using this feature runs the risk of not finding dense subgraphs that overlap 
 the set of excluded edges.
*/

#include <float.h>
#include <limits.h>
#include <math.h>
#include <pthread.h>
#include <search.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>

typedef struct microarrayDataSet_jlong
{
  int numCols, numRows;
  char **labels; /* each label is unique */
  double *data;  /* really a data table  */
} microarrayDataSet;

// Node data struct which takes a string as a label
typedef struct graphNode_jlong
{
  char *label;   /* from microarrayDataSet */
  char **edges;  /* array of graphNode labels that the node connects to */
  double *edgeWeights;
  int numEdges;
} graphNode;

// graph/adj. list data struct. which uses string-labeled nodes (graphNode)
typedef struct gnGraph_jlong
{
  graphNode **gnList;	//adj. list of nodes
  int numNodes;
} gnGraph;

typedef struct edgeSupportVector_jlong
{
  char *edgeName;
  double *supportVector;
} edgeSupportVector;

typedef struct edgeSupportMatrix_jlong
{
  edgeSupportVector **esvList; /* NULL terminated */
  int vecLen;                  /* length of supportVectors */
} edgeSupportMatrix;

// Node data struct which takes an int as a label
/* structures for more efficient searching */
typedef struct gintNode_jlong
{
  int label;
  int *edges;
  double *edgeWeights;
  int numEdges;
} gintNode;

// graph/adj. list data struct. which uses int-labeled nodes (graphNode)
typedef struct gintGraph_jlong
{
  gintNode **nodeList;
  int numNodes;
} gintGraph;

//singly-linked-list of gintGraph
typedef struct subGraphChain_jlong subGraphChain;
struct subGraphChain_jlong
{
  gintGraph *subG;
  int *lmbs;	//stores a sorted list of boundary nodes (nodes which are connected to some node in the subgraph subG) for the last subgraph
  subGraphChain *next;
};

typedef struct thread_args_jlong
{
  int idx;
  int minNodes;
  int numThreads;
  double desiredDensity;
} thread_args;

#ifdef EXCLUDED_EDGES
gnGraph ** ODES(gnGraph *bigGraph, 
                double density, 
                int minNodes,
                char ** excludedNodes);
int compIntA (const void *p1, const void *p2);
int compIntA2(const void *p1, const void *p2);
#else
gnGraph ** ODES(gnGraph *bigGraph, 
                double density, 
                int minNodes);
#endif

/* utility functions defined in util.c */
int compGN      (const void *p1, const void *p2);
int compInt     (const void *p1, const void *p2);
int compIntGintN(const void *p1, const void *p2);
int compStrLM   (const void *p1, const void *p2);
double Density(gnGraph *gp);
void FreeGraph(gnGraph *gp);
void PrintGraph(gnGraph *gp);
gnGraph * ReadGraph(char *filename);

