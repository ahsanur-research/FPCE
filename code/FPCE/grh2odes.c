/* PCE: Pseudo Clique Enumerater */
/* ver 1.0 1/Sep/2005 Takeaki Uno   e-mail:uno@nii.jp,
    homepage:   http://research.nii.ac.jp/~uno/index.html  */
/* This program is available for only academic use, basically.
   Anyone can modify this program, but he/she has to write down
    the change of the modification on the top of the source code.
   Neither contact nor appointment to Takeaki Uno is needed.
   If one wants to re-distribute this code, do not forget to
    refer the newest code, and show the link to homepage of
    Takeaki Uno, to notify the news about the code for the users.
   For the commercial use, please make a contact to Takeaki Uno. */

#ifndef _pce_c_
#define _pce_c_

#define WEIGHT_DOUBLE

#include "alist.c"
#include "sgraph.c"
#include "problem.c"
#include "core.c"
#include "subgraph.c"
#include <math.h>
#include <stdlib.h>

// int PCE_connected=1;
void PCE_error (){
  ERROR_MES = "command explanation";
  print_err ("convert2odes_graph: CMqsS [options] input-filename threshhold [output-filename]\n\
%%:show progress, _:no message, +:write solutions in append mode\n\
C:enumerate pseudo cliques, M:enumerate maximal pseudo cliques\n\
Q,f:print density preceding/following to pseudo cliques\n\
[options]\n\
-K [num]:output [num] most frequent itemsets\n\
-l [num]:output pseudo cliques with size at least [num]\n\
-u [num]:output pseudo cliques with size at most [num]\n\
-U [real-num]:upper bound for density\n\
-# [num]:stop after outputting [num] solutions\n\
-, [char]:give the separator of the numbers in the output\n\
-Q [filename]:replace the output numbers according to the permutation table given by [filename]\n\
threshold is given by real-number from 0 to 1\n\
# the 1st letter of input-filename cannot be '-'.\n\
# if the output file name is '-', the solutions will be output to standard output.\n");
  EXIT;
}
// c:enumerate non-connected pseudo cliques\n");
// -w [filename]:read weights of edges from [filename]\n");

/***********************************************************************/
/*  read parameters given by command line  */
/***********************************************************************/
void PCE_read_param (int argc, char *argv[], PROBLEM *PP){
  ITEMSET *II = &PP->II;
/*  int c=1;
  if ( argc < c+3 ){ PCE_error (); return; }

  if ( !strchr (argv[c], '_') ){ II->flag |= SHOW_MESSAGE; PP->SG.flag |= SHOW_MESSAGE; }
  if ( strchr (argv[c], '%') ) II->flag |= SHOW_PROGRESS;
  if ( strchr (argv[c], '+') ) II->flag |= ITEMSET_APPEND;
  if ( strchr (argv[c], 'f') ) II->flag |= ITEMSET_FREQ;
  if ( strchr (argv[c], 'Q') ) II->flag |= ITEMSET_PRE_FREQ;
  if ( strchr (argv[c], 'C') ) PP->problem |= PROBLEM_FREQSET;
//  if ( strchr( argv[c], 'c' ) ) PP->problem |= FREQSET;
  else if ( strchr (argv[c], 'M') ) PP->problem |= PROBLEM_MAXIMAL;
  else error ("M or C has to be specified", EXIT);
  c++;

  //AuthorA new start
  II->lb = 3; //default value of lower-bound of the order of quasi-clique we wanna find
  //II->ub = INT_MAX;//default value of upper-bound of the order of quasi-clique we wanna find
  //AuthorA new end


  while ( argv[c][0] == '-' ){
    switch (argv[c][1]){
      case 'K': II->topk.end = atoi (argv[c+1]);
      break; case 'l': II->lb = atoi (argv[c+1]);
      break; case 'u': II->ub = atoi(argv[c+1]);
      break; case 'U': II->frq_ub = (WEIGHT)atof(argv[c+1]);
      break; case 'w': PP->weight_fname = argv[c+1];
      break; case '#': II->max_solutions = atoi(argv[c+1]);
      break; case ',': II->separator = argv[c+1][0];
      break; case 'Q': PP->outperm_fname = argv[c+1];
      break; default: goto NEXT;
    }
    c += 2;
    if ( argc<c+2 ){ PCE_error (); return; }
  }
*/
  NEXT:;
  PP->SG.fname = argv[1];
 // if ( II->topk.end==0 ) PP->dense = II->frq_lb = (WEIGHT)atof(argv[c+1]);
//  if ( argc>c+2 )//    PP->output_fname = argv[2];
}
/*  print graph by numbers  */void SGRAPH_print_odesGraph (SGRAPH *G){  VEC_ID i, j;  QUEUE_INT e;  printf ("%d\n" , SGRAPH_NODE_NUM(*G));  FLOOP (i, 0, SGRAPH_NODE_NUM(*G)){    printf ("%d\t", i);//    if ( G->edge.v && G->edge.v[i].t ){      for (j=0; j<G->edge.v[i].t ; j++){        e = G->edge.v[i].v[j];        printf ("%d ", e);      }      putchar ('\n');//    }  }}
/****************/
/* main routine */
/****************/
int PCE_main (int argc, char *argv[]){
  PROBLEM PP;
  SGRAPH *G = &PP.SG;

PROBLEM_init (&PP);
  PCE_read_param (argc, argv, &PP);
if ( ERROR_MES ) return (1);
  G->flag |= LOAD_PERM + LOAD_INCSORT + LOAD_RM_DUP + LOAD_EDGE;
  PP.II.flag |= ITEMSET_ADD;
  PROBLEM_load (&PP);  //SGRAPH_alloc()  //SGRAPH_load (G);  SGRAPH_print_odesGraph(G);


END_MAIN:;

  PROBLEM_end (&PP);

  return (ERROR_MES?1:0);
}

/*******************************************************************************/
#ifndef _NO_MAIN_
#define _NO_MAIN_
int main (int argc, char *argv[]){
  return (PCE_main (argc, argv));
}
#endif
/*******************************************************************************/

#endif
