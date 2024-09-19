/*
 ODES_slave.c                                            sudowrestler@gmail.com
 
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

//#define DEBUG_SLAVE

void FreeSubGraphChain(subGraphChain *sgc);
void GIntGraphs2GnGraphs(int idx);
int ProcessSubGraph(int idx, gintGraph *subG, int *limbs, int minNodes,
                    double desiredDensity);

extern int **hits[NUM_THREADS], lenL[NUM_THREADS], product[NUM_THREADS], 
             rLen[NUM_THREADS];
extern subGraphChain *L[NUM_THREADS], *masterCurr, *sgcCurr[NUM_THREADS];

/*
 Each thread will consume items in its buffer, placing productions in the next
 buffer mod NUM_THREADS.
 
 Algorithm:
 1) watch linked list for subGraphChain unit
 2) save pointer p to subGraphChain unit, and grab the subG
 3) process subG, if extendable, place at end of linked list for next thread
 4) if subG==NULL, and no subG of current order was extended, exit
 5) otherwise, set location of linked list to watch at p->next, and goto 1)
 
 Don't use any optimization when compiling, as sgcThis = sgcLast->next
 gets optimized away.
*/

void ThreadFunc(thread_args *ta)
{
  int i, idx, minNodes, numThreads;
  double density;
  subGraphChain *sgcLast, *sgcPt, *sgcThis;
  
#ifdef DEBUG_SLAVE
  int count=0;
  char buf[16];
  FILE *fp;
  
  printf("howdy from %d\n", ta->idx); fflush(NULL);
  sprintf(buf, "%d_out", ta->idx);
  buf[6] = '\0';
  fp = fopen(buf, "w");
#endif
  
  density    = ta->desiredDensity;
  idx        = ta->idx;
  minNodes   = ta->minNodes;
  numThreads = ta->numThreads;
  
  hits[idx]    = NULL; /* subGs of a set order that were extended */
  lenL[idx]    = 0;    /* length of hits */
  product[idx] = 0;    /* was at least one subG extended? */
  rLen[idx]    = 0;    /* length of allocated hits buffer */
  
  /* sgcCurr is where to place extended subG */
  if((idx+1)%numThreads == 0)	//true for last/(numThreds-1)-th thread
    sgcCurr[idx] = masterCurr; 
  else
    sgcCurr[idx] = L[(idx+1)%numThreads];

  sgcLast = L[idx]; /* watch at sgcLast->next for next one */
  
  /*
   grow any dense subgraphs, i.e. search for every v that satisfies
   den(subG + v) >= density >= 1/2 
 */
  while(1)
  {
    sgcThis = sgcLast->next;
    
    if(sgcThis != NULL)
    {
#ifdef DEBUG_SLAVE
      if(sgcThis == sgcLast) fprintf(fp, "%d)***********SAME ONE!!!\n", idx);
#endif
      FreeSubGraphChain(sgcLast);
      sgcLast = sgcThis;
      
#ifdef DEBUG_SLAVE
      fprintf(fp, "%d)***********got one, count = %d\n", idx, count); fflush(NULL);
#endif

      if(sgcThis->subG==NULL) /* all subGs of current order have been processed */
      {
      
#ifdef DEBUG_SLAVE
        fprintf(fp, "%d) **********enqueue NULL frame\n", idx); fflush(NULL);
#endif

        /* pass the NULL frame along */
        sgcPt = (subGraphChain *) malloc(sizeof(subGraphChain));
        if(sgcPt==NULL)
        {
          fprintf(stderr, "ODES: slave malloc error, thread exiting...\n");
          fflush(stderr);
          pthread_exit(0);
        }
        sgcPt->subG = NULL;
        sgcPt->lmbs = NULL;
        sgcPt->next = NULL;
        
        sgcCurr[idx]->next = sgcPt;
        sgcCurr[idx]       = sgcPt;
        
        if(!product[idx]) /* exit */
        {
        
#ifdef DEBUG_SLAVE
          fprintf(fp, "%d) **********will terminate\n", idx); fflush(NULL);
#endif

          GIntGraphs2GnGraphs(idx);
          
          pthread_exit(0);
        }
        else         /* look for more work */
        {
        
#ifdef DEBUG_SLAVE
          fprintf(fp, "%d) **********continue\n", idx); fflush(NULL);
#endif

          for(i=0; i<lenL[idx]; i++)
            free(hits[idx][i]);
        
          lenL[idx]    = 0;
          hits[idx]    = NULL;
          product[idx] = 0;
          rLen[idx]    = 0;
         
          continue;
        }
      }// if(sgcThis->subG==NULL) /* all subGs of current order have been processed */
      
      if(ProcessSubGraph(idx, sgcThis->subG, sgcThis->lmbs, minNodes, density))
        pthread_exit(0);
      
#ifdef DEBUG_SLAVE
      count++;
#endif

    }//if not-NULL
  }//while
}
