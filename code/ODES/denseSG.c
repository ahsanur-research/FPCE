/*
 denseSG.c                                               sudowrestler@gmail.com
 
 A driver for ODES, an algorithm to find dense subgraphs with density >= 1/2
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

#ifndef ODES_H
  #include "odes.h"
  #define ODES_H
#endif

//#include "declaration.h"

int returnLen = 0; //indicates the no. of clusters found (set in ODES function in ODES_master.c)
int DONTEXTEND = 0;	//Author-A: =1 => report only fixed size clusters, where fixed size == min-nodes param
					// =2 => report all maximal clusters with: min-nodes-param <= size <= max-nodes-param
int CSVPRINT = 0;	//Author-A: = 1 => print clusters in a csv formatted file
int STATPRINT = 0; //Author-A:  =1 => print densities & sizes of the clusters as well, when CSVPRINT == 1

int maxNodes = 0; //Author-A: used only when DONTEXTEND == 2 and <max-nodes> param is specified

int NUM_SORT = 0; //Author-A: =1 if num param is mentioned which indicates that the vertex labels are numbers; currntly, it can only be used with csv option


void PrintUsage(void)
{
  printf("usage: denseSG <graph file> <min density> <min nodes> <option> [max-nodes]\n");
  printf("       <graph file>  - a graph description file, see ReadGraph() for format.\n");
  printf("       <min density> - min density for a subgraph\n");
  printf("       <min nodes>   - min number of nodes for a subgraph\n");
  printf("       <option>  - use this to control clustering: fs* = fixed size clusters of size min-nodes\n\tmx = maximal clusters of size >= min-nodes\nto control printing: *csv* = print output in csv format (default is adj. list format), only if you want fixed sized (of size <min nodes>) clusters\n\t*num* = consider the vertex-labels to be numbers and sort them numerically in the output\n\t*stat = print statistics about each cluster along with the cluster nodes (only use with csv option)\nE.g. option: fscsvstat");
  printf("       [max-nodes]  - use this only if you want maximal clusters of size <= max-nodes, use only with mx option, i.e. option must starts with mx\n");
}

int main(int argc, char *argv[])
{
  int i, minNodes;
  double density;
  gnGraph *g, **denseSubGraphSet;
  
#ifdef EXCLUDED_EDGES /* experimental, set in odes.h */
  char * exclude[]={"10 1","0 10","1 0",'\0'}; /* must be NULL terminated */
#endif
  
  if(argc < 5 || argc > 6)
  {
    PrintUsage();
    return 1;
  }
  
//  if(argc >= 5)
//  {
  if(strstr(argv[4], "fs"))	//==fixed-size
  {
	  /// Author-A: set this when you want fixed sized (of minNodes) clusters
	  DONTEXTEND = 1;
  }
  else if(strstr(argv[4], "mx"))	//==maximal
  {
	  if(argc == 6)
		  DONTEXTEND = 2;
	  else
		  DONTEXTEND = 0;
  }
  if(strstr(argv[4], "csv"))
  {
	  CSVPRINT = 1;
  }
  if(strstr(argv[4], "num"))
  {
	  NUM_SORT = 1;
  }
  if(strstr(argv[4], "stat"))
  {
	  //CSVPRINT = 1;
	  STATPRINT = 1;
  }
//	  printf("DONTEXTEND = %d, CSVPRINT = %d, STATPRINT = %d\n", DONTEXTEND, CSVPRINT, STATPRINT);
//	  return 1;
//  }
  if(argc == 6)
  {
	  if(!strstr(argv[4], "mx"))
	  {
		  fprintf(stderr, "option must start with mx, exiting...");
		  return 1;
	  }

	  maxNodes = (int) strtol(argv[5], NULL, 10);
  }

  density = strtod(argv[2], NULL);
  if(density<0.5 || density>1.0)
  {
    fprintf(stderr, "denseSG: density out of range, must be 0.5 <= den <= 1.0, returning...\n");
    return 1;
  }
  
  minNodes = (int) strtol(argv[3], NULL, 10);
  if(minNodes < 2)
  {
    fprintf(stderr, "denseSG: minNodes must be > 1, and should be > 3, exiting...\n");
    return 1;
  }

  if(argc == 6 && maxNodes < minNodes)
  {
	  fprintf(stderr, "#max-nodes must be >= #min-nodes...");
	  return 1;
  }

//  printf("===============================================================\n");
  
  g = ReadGraph(argv[1]); //PrintGraph(g); continue;
// printf("===============================================================\n");

#ifdef EXCLUDED_EDGES
  denseSubGraphSet = ODES(g, density, minNodes, exclude);
#else
  denseSubGraphSet = ODES(g, density, minNodes);
#endif

  if(denseSubGraphSet!=NULL)
  {
//    printf("\nDense subgraphs of: %s\n", argv[1]);
//    printf("===================\n");
//printf("#Size\tDensity\tNodes\n");
    i = 0;
    if(DONTEXTEND != 0){
		while(denseSubGraphSet[i] != NULL)
		{
			if(DONTEXTEND == 1)
			{
				if(denseSubGraphSet[i]->numNodes == minNodes)
					PrintGraph(denseSubGraphSet[i++]);
			}
			else if(DONTEXTEND == 2)
			{
				if(denseSubGraphSet[i]->numNodes >= minNodes && denseSubGraphSet[i]->numNodes <= maxNodes)
					PrintGraph(denseSubGraphSet[i++]);
			}
			/*else
			{
				PrintGraph(denseSubGraphSet[i++]);
			}*/
		}
	}
	else{
		printf("Found %d clusters", returnLen);	
	}

    FreeGraph(g);
    i = 0;
    while(denseSubGraphSet[i]!=NULL)
    {
      FreeGraph(denseSubGraphSet[i]);
      denseSubGraphSet[i] = NULL;
      i++;
    }
    
//    printf("Done\n\n");
  }
  else
  {
    fprintf(stderr, "denseSG: denseSubGraphSet is NULL for %s, returning...\n", argv[1]);
    return 0;
  }
  
  
  return 0;
}
