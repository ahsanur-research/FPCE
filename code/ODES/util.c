/*
 util.c                                                  sudowrestler@gmail.com
 
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

#ifndef ODES_H
  #include "odes.h"
  #define ODES_H
#endif

//#include "declaration.h"

#define BUFSZ 65536

extern int DONTEXTEND;	//Author-A: =1 => report only fixed size clusters
extern int CSVPRINT;	//Author-A: = 1 => print clusters in a csv formatted file
extern int STATPRINT; //Author-A:  =1 => print densities & sizes of the clusters as well, when CSVPRINT == 1
extern int NUM_SORT;

/* qsort & bsearch compare functions, all are case insensitive */
int compGN(const void *p1, const void *p2)
{
  return strcasecmp((*(graphNode **)p1)->label, (*(graphNode **)p2)->label);
}

int compInt(const void *p1, const void *p2)
{
  return *(int *)p1 - *(int *)p2;
}

#ifdef EXCLUDED_EDGES /* set in odes.h */
int compIntA(const void *p1, const void *p2)
{
  if((*(int **)p1)[0] == (*(int **)p2)[0])
    return (*(int **)p1)[1] - (*(int **)p2)[1];
  else
    return (*(int **)p1)[0] - (*(int **)p2)[0];
}

int compIntA2(const void *p1, const void *p2)
{
  if(((int *)p1)[0] == (*(int **)p2)[0])
    return ((int *)p1)[1] - (*(int **)p2)[1];
  else
    return ((int *)p1)[0] - (*(int **)p2)[0];
}
#endif

int compIntGintN(const void *p1, const void *p2)
{
  return *(int *)p1 - (*(gintNode **)p2)->label;
}

int compStrLM(const void *p1, const void *p2)
{
  return strcasecmp((char *)p1, *(char **)p2);
}

double Density(gnGraph *gp)
{
  int i, j;
  double totalWt;
  
  if(gp->numNodes > 1)
  {
	totalWt = 0;
    for(i=0; i<gp->numNodes; i++)
    {
    	for(j=0;j<gp->gnList[i]->numEdges;j++)
    	{
    		totalWt += gp->gnList[i]->edgeWeights[j];
    	}
    }

    totalWt /= 2; /* each edge was counted twice */
  
    return ((double)totalWt/(double)(((gp->numNodes)*(gp->numNodes-1))/2.0));
  }
  else
    return (0.0);
}

/* free everything that gp points to, and gp itself */
void FreeGraph(gnGraph *gp)
{
  int i, j;
  graphNode *p;

  for(i=0; i<gp->numNodes; i++)
  {
    p = gp->gnList[i];
    
    if(p->label      !=NULL) free(p->label);
    if(p->edgeWeights!=NULL) free(p->edgeWeights);
    
    if(p->edges != NULL)
    {
      for(j=0; j<p->numEdges; j++)
        if(p->edges[j]!=NULL)
          free(p->edges[j]);
        
      free(p->edges);
    }
    
    free(p);
  }
  
  free(gp->gnList);
  free(gp);
  return;
}


/*Author-A: adding radix sort code*/
int getMax(int arr[], int n) 
{ 
    int mx = arr[0]; 
    for (int i = 1; i < n; i++) 
        if (arr[i] > mx) 
            mx = arr[i]; 
    return mx; 
} 
  
// A function to do counting sort of arr[] according to 
// the digit represented by exp. 
void countSort(int arr[], int n, int exp) 
{ 
    int output[n]; // output array 
    int i, count[10] = {0}; 
  
    // Store count of occurrences in count[] 
    for (i = 0; i < n; i++) 
        count[ (arr[i]/exp)%10 ]++; 
  
    // Change count[i] so that count[i] now contains actual 
    //  position of this digit in output[] 
    for (i = 1; i < 10; i++) 
        count[i] += count[i - 1]; 
  
    // Build the output array 
    for (i = n - 1; i >= 0; i--) 
    { 
        output[count[ (arr[i]/exp)%10 ] - 1] = arr[i]; 
        count[ (arr[i]/exp)%10 ]--; 
    } 
  
    // Copy the output array to arr[], so that arr[] now 
    // contains sorted numbers according to current digit 
    for (i = 0; i < n; i++) 
        arr[i] = output[i]; 
} 
  
// The main function to that sorts arr[] of size n using  
// Radix Sort 
void radixsort(int arr[], int n) 
{ 
    // Find the maximum number to know number of digits 
    int m = getMax(arr, n); 
  
    // Do counting sort for every digit. Note that instead 
    // of passing digit number, exp is passed. exp is 10^i 
    // where i is current digit number 
    for (int exp = 1; m/exp > 0; exp *= 10) 
        countSort(arr, n, exp); 
} 

/* Author-A added radix-sort code above */

/* mainly for debugging */
void PrintGraph(gnGraph *gp)
{
	if(CSVPRINT)
	{
	  int i;
	  graphNode *p;
	  if(STATPRINT)
		  printf("%d\t%1.6lf\t", gp->numNodes, Density(gp));

		if(!NUM_SORT){
			for(i=0; i<gp->numNodes; i++)
			{
			p = gp->gnList[i];
			printf("%s ", p->label);
			}
			printf("\n");
	  }
	  else{
	  	int *a = (int *)malloc(gp->numNodes * sizeof(int));
			for(i=0; i<gp->numNodes; i++)
			{
			  p = gp->gnList[i];
			  a[i] = atoi(p->label);
			}
	  	radixsort(a,gp->numNodes);//sort numerically
			for(i=0; i<gp->numNodes; i++)
			{
  			printf("%d ", a[i]);
			}
			printf("\n");
	  }
	}
	else
	{
	  int i, j, len;
	  graphNode *p;

	  printf("\n    Density = %1.5lf\n", Density(gp));
	  printf("    Node        Edges                   Weights\n");
	  printf("    ===========================================\n");

	  for(i=0; i<gp->numNodes; i++)
	  {
		p = gp->gnList[i];

		printf("%3d) %s\t\t", i, p->label);

		len = 0;
		if(p->edges!=NULL)
		{
		  for(j=0; j<p->numEdges; j++)
		  {
			printf("%s ", p->edges[j]);
			fflush(stdout);
			len += 1 + strlen(p->edges[j]);
		  }

		  if(len < 8)
			printf("\t\t\t");
		  else if(len < 16)
			printf("\t\t");
		  else
			printf("\t");

		  if(p->edgeWeights!=NULL)
			for(j=0; j<p->numEdges; j++)
			  printf("%3.1lf ", p->edgeWeights[j]);

		}
		printf("\n");
	  }

	  printf("    ===========================================\n\n");
	}
  fflush(stdout);
}

/*
 Format for infile:
 Number of nodes is on first line; subsequent lines are a series of tokens,
 where the first token on a line is the node label, and subsequent tokens are
 the edges for that node, i.e. the labels of the nodes this node connects to.
 Edge weights are appended to the label with a '__', i.e. label__0.7
 Comments are preceeded by a '#', thus no label may contain a '#' sign
*/
//# define CONSISTENCY_CHECK
gnGraph * ReadGraph(char *filename)
{
  int i, firstTok;
  char *b, buf[BUFSZ], *p, *w;
  double edgeWeight=0.0;
  gnGraph *g;
  graphNode *gn;
  FILE *infile;
  
#ifdef CONSISTENCY_CHECK
  int j, k, m, found;
#endif
  
  infile = fopen(filename, "r");
  if(infile==NULL)
  {
    fprintf(stderr, "ReadGraph: unable to open %s, returning...\n", filename);
    return NULL;
  }
    
  g = (gnGraph *) malloc(sizeof(gnGraph));
  if(g==NULL)
  {
    fprintf(stderr, "ReadGraph: malloc error, returning...\n");
    return NULL;
  }
  
  g->numNodes = 0;
  while(fgets(buf, BUFSZ, infile))
    if((g->numNodes = (int) strtol(buf, NULL, 10)))
      break;
  
  if(g->numNodes)
  {
    g->gnList = (graphNode **) malloc(g->numNodes * sizeof(graphNode *));
    if(g->gnList==NULL)
    {
      fprintf(stderr, "ReadGraph: malloc error, returning...\n");
      return NULL;
    }
  
    for(i=0; i<g->numNodes; i++)
    {
      gn = g->gnList[i] = (graphNode *) malloc(sizeof(graphNode));
      if(gn==NULL)
      {
        fprintf(stderr, "ReadGraph: malloc error, returning...\n");
        return NULL;
      }
      
      gn->edges       = NULL;
      gn->edgeWeights = NULL;
    }
    
    i = 0;
    while(fgets(buf, BUFSZ, infile))
    {
      if(buf[0]=='\n') continue;
      if((p = strstr(buf, "#")))  p[0] = '\0';
      if((p = strstr(buf, "\n"))) p[0] = '\0';
      
      b = &buf[0];
      p = &buf[0]; 
      firstTok = 1;
      while(p!=NULL && b!=NULL)
      {
        p = strsep(&b, " ,\t");
        if(p[0]=='\0') continue;
        
        gn = g->gnList[i];
        if(firstTok) /* the node label */
        {
          gn->label = (char *) malloc(strlen(p)+1);
          if(gn->label==NULL)
          {
            fprintf(stderr, "ReadGraph: malloc error, returning...\n");
            return NULL;
          }
          strcpy(gn->label, p);
          
          firstTok = 0;
        }
        else /* the edges */
        {
	  if((w = strstr(p, "__"))!=NULL) 
	  {
	    w[0] = '\0';
	    w += 2;
	    edgeWeight = strtod(w, NULL);
//printf("label %s ew = %lf\n", p, edgeWeight);
	  }

          if(gn->edges==NULL)
          {
            gn->edges       = (char **)  malloc(sizeof(char *));
            gn->edgeWeights = (double *) malloc(sizeof(double));
            if(gn->edges==NULL || gn->edgeWeights==NULL)
            {
              fprintf(stderr, "ReadGraph: malloc error, returning...\n");
              return NULL;
            }
            gn->numEdges = 1;
            
            gn->edges[0] = (char *) malloc(strlen(p)+1);
            if(gn->edges[0]==NULL)
            {
              fprintf(stderr, "ReadGraph: malloc error, returning...\n");
              return NULL;
            }
            strcpy(gn->edges[0], p);
	    
			if(w!=NULL) gn->edgeWeights[0] = edgeWeight;
			else        gn->edgeWeights[0] = 1.0;
          }
          else
          {
            gn->numEdges++;
            gn->edges       = (char **)  realloc(gn->edges,       gn->numEdges*sizeof(char *));
            gn->edgeWeights = (double *) realloc(gn->edgeWeights, gn->numEdges*sizeof(double));
            if(gn->edges==NULL || gn->edgeWeights==NULL)
            {
              fprintf(stderr, "ReadGraph: realloc error, returning...\n");
              return NULL;
            }
            
            gn->edges[gn->numEdges-1] = (char *) malloc(strlen(p)+1);
            if(gn->edges[0]==NULL)
            {
              fprintf(stderr, "ReadGraph: malloc error, returning...\n");
              return NULL;
            }
            strcpy(gn->edges[gn->numEdges-1], p);
	    
			if(w!=NULL) gn->edgeWeights[gn->numEdges-1] = edgeWeight;
			else        gn->edgeWeights[gn->numEdges-1] = 1.0;
          }
        }
      }
      
      i++;
    }
  }
  else
  {
    fprintf(stderr, "ReadGraph: Error, numNodes is zero, returning...\n");
    return NULL;
  }
  
#ifdef CONSISTENCY_CHECK
  /*
   check that every claimed edge & weight is claimed by both vertices
   very slow for large graphs
  */
  for(i=0; i<g->numNodes; i++)              /* start node */
    for(j=0; j<g->gnList[i]->numEdges; j++) /* edges      */
      for(k=0; k<g->numNodes; k++)          /* end node   */
      {
        if(!strcasecmp(g->gnList[i]->edges[j], g->gnList[k]->label))
        {
          found = 0;
          for(m=0; m<g->gnList[k]->numEdges; m++)
          {
            if(!strcasecmp(g->gnList[k]->edges[m], g->gnList[i]->label))
            {
              found = 1;
              if(g->gnList[i]->edgeWeights[j] != g->gnList[k]->edgeWeights[m])
              {
                fprintf(stderr, "ReadGraph: error in file %s, vertex %s claims edgeWeight %lf to %s\n", 
                                 filename, g->gnList[i]->label, 
                                 g->gnList[i]->edgeWeights[j], 
                                 g->gnList[i]->edges[j]);
                return NULL;
              }
              
              break;
            }
          }
          if(!found)
          {
            fprintf(stderr, "ReadGraph: error in file %s, vertex %s claims edge to %s\n", 
                    filename, g->gnList[i]->label, g->gnList[i]->edges[j]);
            return NULL;
          }
        }
      }
#endif

  return g;
}
