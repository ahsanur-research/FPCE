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

//AuthorA new start
//#define PRINT_TREE 1
//#undef PRINT_TREE

#define  INF  1 << 30
//#define THREAD 100


//SGRAPH *subg;//Author-A:to use in getSubgraph function and
            //to generate children of current quasi-clique till PP->frq_lb = 1

int    degen        =  0; // degeneracy of the whole graph (calculated in PCE_main)
//int    degen_cutoff = -1; // degeneracy cutoff to calculate elements eligible roots array
//double LOCAL_DEGEN_CO_NUMER;


//AuthorA: used to print #of calls of PCE_iter before/after applying different bounds
unsigned long long int numCallPCE_Iter=0;


// AuthorK: For connected quasi-cliques
#define CONNECTED_CLIQUE 1
//#undef  CONNECTED_CLIQUE


int overAllDegree = 0;

// AuthorK: minimum value among three numbers

#define TWO_MAX(x,y) ((x)>(y)?(x):(y))
#define TWO_MIN(x,y) ((x)<(y)?(x):(y))
#define THREE_MIN(p,q,r) TWO_MIN(TWO_MIN(p,q),r)

// AuthorK: minumum value among three numbers


// AuthorK: K-core

int k_core   = -1;	   // maximum core of a graph
int min_core = -1; 	   // minimum core of a graph

// AuthorK: K-core


// AuthorK: core-no calculation of a vertex

int *degreeStatus; // will marked remove vertex ( size = #of nodes / 32 )
int *coreStatus;   // will contain core-no of a vertex

// AuthorK: core-no calculation of a vertex

#define STRING_SIZE 505
/*
FILE *filePointer;
char *commandlineString;
char *temporaryString;
char *finalString;
char *kcoreString;
*/
/***********************************************************************/
/*  PCE bounds  */
/***********************************************************************/

//AuthorA: define to apply a bound, undef to disable it
//#define APPLY_MAXDEG_BOUND 1;
//#undef APPLY_MAXDEG_BOUND

//#define APPLY_MINDEG_BOUND 1;
//#undef APPLY_MINDEG_BOUND

//#define APPLY_STRICT_MINDEG_BOUND 1;
//#undef APPLY_STRICT_MINDEG_BOUND

//#define APPLY_MAXSIZE_BOUND 1;
//#undef APPLY_MAXSIZE_BOUND


/*
// AuthorK: We don't need this extra #ifdef APPLY_DEL_BY_MU_BOUND condition
// 		   as we calculated upper bound in PCE-main

// AuthorK: Del by mu bound

//#define APPLY_DEL_BY_MU_BOUND 1;
//#undef  APPLY_DEL_BY_MU_BOUND

#ifdef APPLY_DEL_BY_MU_BOUND
    unsigned long long int numcalls_saved_by_del_by_mu_bound = 0;
#endif // APPLY_DEL_BY_MU_BOUND

// AuthorK: Del by mu bound
*/


// AuthorK: Edge bound
//#define APPLY_SIBLING_BOUND 1;
//#undef APPLY_SIBLING_BOUND;
#ifdef APPLY_SIBLING_BOUND
    unsigned long long int numcalls_saved_by_sibling_bound = 0;
#endif // APPLY_SIBLING_BOUND

//#define APPLY_EDGE_BOUND
#undef  APPLY_EDGE_BOUND

#ifdef APPLY_EDGE_BOUND
    int LHS; //AuthorA new
    unsigned long long int numcalls_saved_by_edge_bound = 0;
#endif // APPLY_EDGE_BOUND

//#define APPLY_ORDER_BOUND
#undef APPLY_ORDER_BOUND
#ifdef APPLY_ORDER_BOUND
unsigned long long int numcalls_saved_by_order_bound = 0;
#endif // APPLY_ORDER_BOUND


/*
// AuthorK: Local degeneracy bound

//#define APPLY_LOCAL_DEGEN_BOUND 1
//#undef  APPLY_LOCAL_DEGEN_BOUND

#ifdef APPLY_LOCAL_DEGEN_BOUND
	unsigned long long int numcalls_saved_by_local_degen_bound = 0;
#endif // APPLY_LOCAL_DEGEN_BOUND

// AuthorK: Local degeneracy bound
*/


// AuthorK: Turan filtering

int    R;
double storeTheta;

#define APPLY_TURAN_FILTERING 1
//#undef  APPLY_TURAN_FILTERING

// AuthorK: Turan filtering


// AuthorK: Neighbor bound

//#define APPLY_NEIGHBOR_BOUND 1;
//#undef  APPLY_NEIGHBOR_BOUND

#ifdef APPLY_NEIGHBOR_BOUND
    unsigned long long int numcalls_saved_by_neighbor_bound = 0;
#endif // APPLY_NEIGHBOR_BOUND

// AuthorK: Neighbor bound


/*
// AuthorK: We don't need this extra #ifdef APPLY_H_INDEX_BOUND condition
// 		   as we calculated h-index in PCE-main

// AuthorK: H-Index bound

//#define APPLY_H_INDEX_BOUND 1;
//#undef  APPLY_H_INDEX_BOUND

#ifdef APPLY_H_INDEX_BOUND
    unsigned long long int numcalls_saved_by_h_index_bound = 0;
#endif // APPLY_H_INDEX_BOUND

// AuthorK: H-Index bound
*/


#ifdef APPLY_MAXSIZE_BOUND
unsigned int maxSize = 0;//max. #of vertices in a pseudo-clique
unsigned long long int numcalls_saved_by_maxsize_bound=0;
#endif // APPLY_MAXSIZE_BOUND

#ifdef APPLY_MINDEG_BOUND
int bound;
unsigned long long int numcalls_saved_by_mindeg_bound=0;
#endif // APPLY_MINDEG_BOUND

#ifdef APPLY_MAXDEG_BOUND
unsigned long long int numcalls_saved_by_maxdeg_bound=0;
#endif // APPLY_MAXDEG_BOUND


// int PCE_connected=1;
void PCE_error (){
  ERROR_MES = "command explanation";
  print_err ("pce: CMqsS [options] input-filename threshhold [output-filename]\n\
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
  int c=1;
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

  NEXT:;
  PP->SG.fname = argv[c];
  if ( II->topk.end==0 ) PP->dense = II->frq_lb = (WEIGHT)atof(argv[c+1]);
  if ( argc>c+2 ) PP->output_fname = argv[c+2];
}


/******************************************************************/
/******************************************************************/
/******************************************************************/

/* add a vertex to clique and update the degrees of other vertices */
void PCE_add_vertex_to_clique (PROBLEM *PP, QUEUE_INT v, double *sum, double *all){
  QUEUE_INT *x;
  *all += PP->II.itemset.t;
  *sum += PP->occ.list[v];
  QUEUE_ins_ele_ (&PP->II.itemset, v);
  MQUE_FLOOP (PP->SG.edge.v[v], x)
      MALIST_mv (&PP->occ, PP->occ.list[*x] + 1, *x, 0);
}
/* remove a vertex to clique and update the degrees of other vertices */
void PCE_rm_vertex_from_clique (PROBLEM *PP, QUEUE_INT v, double *sum, double *all){
  QUEUE_INT *x;
  *sum -= PP->occ.list[v];
  *all -= PP->II.itemset.t;
  QUEUE_rm_ele_ (&PP->II.itemset, v);
  MQUE_FLOOP (PP->SG.edge.v[v], x)
      MALIST_mv (&PP->occ, PP->occ.list[*x] -1, *x, 0);
}
/* compute the minimum degree in the graph induced by PCE_clq */
QUEUE_INT PCE_min_degree_in_clique (PROBLEM *PP){
  QUEUE_INT *x, m=PP->II.itemset.t+1;
  MQUE_FLOOP (PP->II.itemset, x) ENMIN (m, PP->occ.list[*x]);
  return (m);
}


//AuthorA new

// version created by AuthorA for computing degeneracy
/* compute the minimum degree in the graph induced by PCE_clq */
QUEUE_INT PCE_min_degree_in_clique2 (PROBLEM *PP, QUEUE_INT *mindeg_vertex)
{
    QUEUE_INT *x, m=PP->II.itemset.t+1;
    MQUE_FLOOP (PP->II.itemset, x)
    {
    	//printf ( "node %d, degree = %d\n" , *x , PP -> occ.list[ *x ] );

        if (PP->occ.list[*x] < m)
        {
            m = PP->occ.list[*x];
            *mindeg_vertex = *x;
        }
    }
    //printf ( "\n" );

    return (m);
}

//AuthorA new end


/* collect minimum degree vertices in PCE_clq */
void PCE_cllect_degree_vertices (PROBLEM *PP, QUEUE_INT d){
  ITEMSET *II = &PP->II;
  QUEUE_INT *x;
  II->add.t = 0;
  MQUE_FLOOP (II->itemset, x)
    if ( PP->occ.list[*x] == d ) QUE_INS (II->add, *x);
}

//AuthorA new: generate children (K U {v}) of current quasi-clique K s.t.
//min_add_deg = \theta(K) <= deg_{K U {v}} (v) < min_deg (K)
void PCE_enum_type1_children (PROBLEM *PP, QUEUE_INT min_deg, QUEUE_INT min_add_deg){
  ALIST_ID i, j;
  FLOOP (i, min_add_deg, min_deg)
      MALIST_DO_FORWARD (PP->occ, i, j) QUE_INS (PP->itemcand, j);
}


//AuthorA new: generate children (K U {v}) of current quasi-clique K s.t.
//deg_{K U {v}} (v) == min_deg (K)
void PCE_enum_type2_children (PROBLEM *PP, QUEUE_INT min_deg, QUEUE_INT min_add_deg){
  SGRAPH *G = &PP->SG;
  ITEMSET *II = &PP->II;
  QUEUE_INT *x, u, k=0, kk=0;
  QUEUE_ID j;

  if ( min_deg < min_add_deg ) return;

  //AuthorA: probably iterate over vertices u in V[G] (in+out K) with degK = min_deg
  kk=0; MALIST_DO_FORWARD (PP->occ, min_deg, u){
    PP->itemary[u] = 0;
    PP->vecary[kk++] = u;
//    printf ("# %d\n", u);
  }

  //now PP->itemary[u] = 0 for all nodes (in/outside of K) u with degK = mindeg,
  //vecary contains a list of such nodes and kk is the #of such nodes

  qsort_VEC_ID (PP->vecary, kk, 1);//sort such nodes by their IDs

  MQUE_FLOOP (II->itemset, x){//foreach itemset/insider node *x
    PP->itemary[*x] -= G->edge.t;  // not to find vertices in II->itemset
  }

  PCE_cllect_degree_vertices (PP, min_deg); //AuthorA: populate II->add with nodes in K with degK = mindeg

  MQUE_FLOOP (II->add, x){//foreach insider (of K) node *x
//    if ( PP->occ.list[*x] == min_deg ){//if an insider node *x has degK = mindeg then
      //QUE_INS (II->add, *x);//AuthorA: adding here so that I don't have to do it later

      QUEUE_BE_LOOP_ (G->edge.v[*x], j, u){//foreach nbr u of *x taken in reverse order (of their IDs?)
        if ( u<*x ) break;
        PP->itemary[u]++;
//        printf ("#inc++ %d\n", u);
      }
      //now PP->itemary[u]>0 implies  PP->itemary[u] = #of nbrs of u outside K whose IDs > u
      k++;

      if(k >= min_deg) break;
//    }
  }

 //now S={u:PP->itemary[u]>0} contains nodes u outside K for which some potential competitors (nodes with degK=mindeg whose ID<u) exist in K
  //In S:  PP->itemary[u] is the number of such potential competitors; "potential" means not sure-competitors; since such competitor won't be a competitor
  //if it's a nbr of u because in that case deg_KU{u} of that node will be mindeg+1

  QUE_INS (II->add, G->edge.t);

  j=0; FLOOP (k, 0, kk){
    u = PP->vecary[k];
    while ( II->add.v[j]<u ) j++;
//    printf ("u=%d, j=%d,%d\n", u, j, II->add.v[j]);
    if ( PP->itemary[u] == j ) QUE_INS (PP->itemcand, u);
  }
}

//AuthorA new: generate children (K U {v}) of current quasi-clique K s.t.
//deg_{K U {v}} (v) == min_deg (K)+1
void PCE_enum_type3_children (PROBLEM *PP, QUEUE_INT min_deg, QUEUE_INT min_add_deg){
  SGRAPH *G = &PP->SG;
  ITEMSET *II = &PP->II;
  ALIST_ID u;
  int flag=1, vmin=0;
  QUEUE_ID jj=0, j;
  QUEUE_INT *x, *y, k=0, kk=0;

//  if ( min_deg+1 < min_add_deg ) return;
  MALIST_DO_FORWARD (PP->occ, min_deg+1, u){
      PP->itemary[u] = 0;
      PP->vecary[kk++] = u;
  }
  qsort_QUEUE_INT (PP->vecary, kk, 1);


  MQUE_FLOOP (II->itemset, x){
    PP->itemary[*x] -= G->edge.t;  // not to find vertices in PCE_clq
    if ( PP->occ.list[*x] == min_deg ){
      jj++;
      if ( flag ){ flag=0; vmin=*x; }
      MQUE_FLOOP (G->edge.v[*x], y){
        if ( *y >= vmin ) break;
        PP->itemary[*y]++;
//        printf ("inc++ %d\n", u);
      }
    }
  }

  PCE_cllect_degree_vertices (PP, min_deg+1);

  MQUE_FLOOP (II->add, x){
//    if ( PP->occ.list[*x] == min_deg+1 ){
      for ( y=G->edge.v[*x].v+G->edge.v[*x].t-1 ; y>=G->edge.v[*x].v ; y--){
        if ( *y < *x ) break;
        PP->itemary[*y]++;
//        printf ("inc %d\n", u);
      }
      k++;
      if (  k >= min_deg+1 ) break;
//    }
  }

  QUE_INS (II->add, G->edge.t);
  j=0; FLOOP (k, 0, kk){
    u = PP->vecary[k];
//    printf ("u=%d, j=%d, jj=%d, kk=%d, vmin=%d\n", u, j, jj, kk, vmin);
    if ( u>vmin ) break;
    while ( II->add.v[j]<u ) j++;
    if ( PP->itemary[u] == j+jj ) QUE_INS (PP->itemcand, u);
  }
}

QUEUE PCE_enum_children (PROBLEM *PP, QUEUE_INT min_deg, QUEUE_INT min_add_deg){
  PP->itemcand.t = 0;

  //AuthorA: vertices adjacent to K with degree < min_deg
  PCE_enum_type1_children (PP, min_deg, min_add_deg);

  //AuthorA: new subroutine to compute vertices adjacent to K with degree == min_deg or mindeg + 1
  //PCE_enum_type23_children (PP, min_deg, min_add_deg);

  //AuthorA: vertices adjacent to K with degree == min_deg
  PCE_enum_type2_children (PP, min_deg, min_add_deg);

  PCE_enum_type3_children (PP, min_deg, min_add_deg);

  return ( QUEUE_dup_ (&PP->itemcand) );
}

QUEUE PCE_enum_children_clq (PROBLEM *PP, //bitmap *bmp, unsigned lastNode){
                             QUEUE *prevCommonNbrs, unsigned lastNode){
  PP->itemcand.t = 0;
//AuthorA: vertices adjacent to K with degree_in_K == min_deg+1
  //return intersection of prevCommonNbrs and Nbr(lastNode)
//  printf("\nlastNode = %d, prevComNbrs: ", lastNode);
//  if(prevCommonNbrs != NULL)
//    QUEUE_print(prevCommonNbrs);


  unsigned int j, nbr;
  if(prevCommonNbrs == NULL){//initial case
    ITEMSET *II = &PP->II;
    if(II->itemset.t == 1){
      for(j=0; j < PP->subg->edge.v[lastNode].t; j++) {
          nbr = PP->subg->edge.v[lastNode].v[j];
          if(nbr < lastNode)
              QUE_INS (PP->itemcand, nbr);
      }
    }
    else if(II->itemset.t == 2){
      int u = TWO_MIN(II->itemset.v[0], II->itemset.v[1]), v = TWO_MAX(II->itemset.v[0], II->itemset.v[1]);
      for(j=0; j < PP->subg->edge.v[u].t; j++) {
        nbr = PP->subg->edge.v[u].v[j];
        if(nbr < u && SGRAPH_is_edge_bmp(nbr, v))//ensure that u,v, nbr forms a triangle since R >= 3. NOTE: II->itemset.t == 2 and prevCommonNbrs = NULL can only happen when R>=3, as per rules for forming Edge_list in PCE_main
          QUE_INS (PP->itemcand, nbr);
      }
    }
    else{
      print_err("Error: II->itemset.t > 2");
      EXIT;
    }
  } else
  {
//      bitmap bmp;
//      createBitMap(&bmp, PP->subg->edge.v[lastNode].v, PP->subg->edge.v[lastNode].t);
      QUEUE_ID i;
      QUEUE *Q = prevCommonNbrs;
      for ( i=Q->s; i!=Q->t ; ){
        if(Q->v[i] < lastNode && SGRAPH_is_edge_bmp(Q->v[i], lastNode))//isPresent(&bmp, Q->v[i]))//
            QUE_INS (PP->itemcand, Q->v[i]);
        QUEUE_INCREMENT(*Q,i);
      }
//      deallocateBitmap(&bmp);
  }

//  printf(", itemcand: ");
  //QUEUE_print(&PP->itemcand);
  //fflush(stdout);

  return ( QUEUE_dup_ (&PP->itemcand) );
}

int PCE_check_maximality (PROBLEM *PP, QUEUE_INT min_add_deg){
  ITEMSET *II = &PP->II;
  QUEUE_INT *x;
  QUEUE_ID i;
  ARY_FILL (PP->itemary, min_add_deg, II->itemset.t+1, 0);
  MQUE_FLOOP (II->itemset, x) PP->itemary[PP->occ.list[*x]]++;
  FLOOP (i, min_add_deg, II->itemset.t+1)
    if ( PP->occ.num[i] > PP->itemary[i] ) return (0);
  return (1);
}

//#define DEBUG
#undef DEBUG

/*************************************************************************/
/* pseudo clique enumeration, main iteration
   PP = K,
*/
/*************************************************************************/
#ifdef APPLY_SIBLING_BOUND
//returns 1 if there is no subtree under the treenode: K U {v}
unsigned PCE_iter ( PROBLEM *PP, QUEUE_INT v, double sum, double all, int local_degen){
    unsigned noSubtreeUnderCurTreenode = 1; //assume that there is no subtree under KU{v}. set it to false if execution goes to child-generation code and some children is generated
#else
void PCE_iter ( PROBLEM *PP, QUEUE_INT v, double sum, double all,
               int local_degen, QUEUE *prevCh){
#endif // APPLY_SIBLING_BOUND
  QUEUE_INT *x, min_deg, min_add_deg;
  QUEUE jump; //, prevChildren;
  ITEMSET *II = &PP->II;//II==PP->II now holds info of the old clique = K-{v}

//  II->iters++;

  //AuthorA:
  numCallPCE_Iter++;  // AuthorK: PCE_iter function call count


  /*
  #if defined (APPLY_EDGE_BOUND) || defined (APPLY_LOCAL_DEGEN_BOUND)
  int PCE_degen = local_degen;
  #endif // APPLY_EDGE_BOUND || APPLY_LOCAL_DEGEN_BOUND
  */

  //temporarily add vertex v to the old clique to get the clique K
  //after this call, II = PP->II = extended clique = K, sum =|E[K]|
  //all = clq(K) = |K|(|K|-1)/2
  PCE_add_vertex_to_clique (PP, v, &sum, &all);  // add vertex v to the current clique

//printf ("\nv=%d, sum=%0.0f, all = %0.0f, II->itemset.t=%d, Itemset: ", v, sum, all, II->itemset.t);
//QUEUE_print (&(II->itemset));

  /*
  #if defined (APPLY_EDGE_BOUND) || defined (APPLY_LOCAL_DEGEN_BOUND)
  local_degen = coreStatus [ v ];
  if ( local_degen > PCE_degen )
  	local_degen = PCE_degen;
  #endif // APPLY_EDGE_BOUND || APPLY_LOCAL_DEGEN_BOUND
  */


  // AuthorK: Turan Filtering

  #ifdef APPLY_TURAN_FILTERING

  if ( II -> itemset.t < R )
  {
  	II -> frq_lb = 1;
  }

  else
  {
  	II -> frq_lb = storeTheta;
  }

  #endif // APPLY_TURAN_FILTERING

  // AuthorK: Turan Filtering


#ifdef DEBUG
  printf ("v=%d, sum=%0.0f, all = %0.0f, II->itemset.t=%d, Itemset: ", v, sum, all, II->itemset.t);
  QUEUE_print (&(II->itemset));
#endif

  //  printf ("start:  %d:  ", v);
  //  QUEUE_print (&PCE_clq);

  //if |K| < 2 then density(K) = 1, otherwise compute density(K)
  if ( II->itemset.t < 2 ) II->frq = 1.0;
  else II->frq = sum/all;//sum = |E[K]|, all = clq(K) = |K|(|K|-1)/2
  //so II->frq = density(K)

  //min_deg = PCE_min_degree_in_clique (PP);//min_deg = min. deg in II = K
  min_deg = PP->occ.list[v];

  //min_add_deg = theta(K) because sum = |E[K]|, II->frq_lb = theta,
  //all = clq(K) = |K|(|K|-1)/2, II->itemset.t = |K|
  min_add_deg = (QUEUE_INT)ceil (II->frq_lb*(all+II->itemset.t) - sum);
  if ( min_add_deg < 0 ) min_add_deg = 0;

  // AuthorK: For connected quasi-cliques
  #ifdef CONNECTED_CLIQUE
  	  if ( min_add_deg == 0 ) min_add_deg = 1;
  #endif // CONNECTED_CLIQUE

  if(II->itemset.t < R) goto GEN_CHILD;//skip unnecessary codes for pruning-attempts or for checking maximality

  #ifdef APPLY_EDGE_BOUND
  //#if defined (APPLY_EDGE_BOUND) || defined (APPLY_LOCAL_DEGEN_BOUND)

    int PCE_degen = local_degen;
    local_degen   = coreStatus [ v ];

    if ( local_degen > PCE_degen )
  	  local_degen = PCE_degen;

  //#endif // APPLY_EDGE_BOUND || APPLY_LOCAL_DEGEN_BOUND
  #endif // APPLY_EDGE_BOUND


  /*
  // AuthorK: Local degeneracy bound

  #ifdef APPLY_LOCAL_DEGEN_BOUND

    if ( ( II -> itemset.t ) >= ( II -> lb ) )
  	  goto END_LOCAL_DEGEN;

    if ( min_deg > local_d )
  	  local_d = min_deg;

    int dc 			       = TWO_MAX ( local_degen , local_d );
    int local_degen_co = ceil ( ( LOCAL_DEGEN_CO_NUMER - sum ) / ( 1.0 * ( ( II -> lb ) - ( II -> itemset.t ) ) ) );

    if ( dc < local_degen_co )
    {
  	  numcalls_saved_by_local_degen_bound ++;

  	  #ifdef PRINT_TREE
      printf ( "Local degen bound\n" );
	  #endif // PRINT_TREE

  	  goto END;
    }

  #endif // APPLY_LOCAL_DEGEN_BOUND

  END_LOCAL_DEGEN:;

  // AuthorK: Local degeneracy bound
  */


  //AuthorA new start
  #ifdef PRINT_TREE
  //QUEUE_print (&(II->itemset));
  //printf ("density=%f, v=%d, E[K]=%0.0f, KC2=%0.0f, n=%d\n", II->frq, v, sum, all, II->itemset.t);
      int i;

      FLOOP(i,0,II->itemset.t-2)printf("|");
      if(II->itemset.t > 1)printf("%c",192);
  //QUEUE_print (&(II->itemset));
      printf ("%d => (%0.2f, %d, %d): ", v, II->frq, II->itemset.t, min_deg);
      QUEUE_print (&(II->itemset));
  #endif // PRINT_TREE
  //AuthorA new end


  //printf ( "Queue Print ==> " );
  //QUEUE_print (&(II->itemset));

  //if ( II -> frq >= II -> frq_lb && II -> itemset.t >= II -> lb ) QUEUE_print (&(II->itemset));    // AuthorK: Print pseudo clique

  /*
  //Generate Tree

  //QUEUE_print (&(II->itemset));
  //printf ("density=%f, v=%d, E[K]=%0.0f, KC2=%0.0f, n=%d\n", II->frq, v, sum, all, II->itemset.t);
  int i;

  FLOOP(i,0,II->itemset.t-2)printf("|");
  if(II->itemset.t > 1)printf( "%c" , 192 );
  //QUEUE_print (&(II->itemset));
  //printf ("%d => (%0.2f, %d, %d, %.2lf, %d): ", v, II->frq, II->itemset.t, min_deg, sum ,PP->SG.edge.v[v].t);
  printf ("%d => (%0.2f, %d, %d ): ", v, II->frq, II->itemset.t, min_deg );
  if ( II -> frq >= II -> frq_lb && II -> itemset.t >= II -> lb ) printf ( " Pseudo clique ==> " );
  QUEUE_print (&(II->itemset));

// |E|K|| #of vertices in K           = sum
// Current degree of a vertex         = PP->SG.edge.v[v].t
// Minimum degree in K 				        = min_deg
// theta ( K ) 						            = min_add_deg
// Output cliques with size at least  = II -> lb
// n 								                  = II -> itemset.t
// #of nodes 						              = G->edge.t
// clq(K) = |K|(|K|-1)/2 			        = all

  //Generate Tree
  */



  //  printf ("%d %d\n", (PCE_problem&1),PCE_check_maximality( min_add_deg ));
  if ( ( II->itemset.t >= II->lb ) && ( (PP->problem&1) || PCE_check_maximality(PP, min_add_deg ) ) ){
      ITEMSET_output_itemset (II, NULL, 0);
      //output K if it can't be extended further (line 1 of EnumPseudoClique algorithm in paper)
  //AuthorA new
  #ifdef PRINT_TREE
	  printf("*");
  #endif // PRINT_TREE
  }
  #ifdef PRINT_TREE
	  printf("\n");
  #endif // PRINT_TREE
  //AuthorA new end
#ifdef APPLY_ORDER_BOUND
  if ( II->itemset.t == II->ub ) {
//printf(": pruned by order bound\n");
numcalls_saved_by_order_bound++;
    goto END;
  }
  #endif

  #ifdef APPLY_MINDEG_BOUND
  #ifdef APPLY_STRICT_MINDEG_BOUND
      //if min_deg = deg_K(v*) = deg_K(v) is equal to deg(v) in total graph G
      //then total degree of v is contributed in K
      //so the next vertex to be added (w) can't be a neighbor of v, so next value of
      //min_deg = deg_KU{w}(w) = deg_KU{w}(v*(KU{w})) <= deg_K(v) = current min_deg
      //otherwise, next value of min_deg can be at most min_deg+1
      if(min_deg == PP->SG.edge.v[v].t)
          bound = min_deg;
      else
          bound = min_deg+1;
  #else
      bound = min_deg+1;
  #endif // APPLY_STRICT_MINDEG_BOUND
      if(bound<min_add_deg){
       numcalls_saved_by_mindeg_bound++;
      goto END;
      }

    printf ( "%d   " , PP->SG.edge.v[v].t );       //Current degree of a vertex.
  #endif // APPLY_MINDEG_BOUND


  /*
  // AuthorK: Del by mu bound permanently disable

  // AuthorK: Del by mu bound

  #ifdef APPLY_DEL_BY_MU_BOUND

      if ( ( II -> itemset.t + 1 ) > ( ( sum + min_deg + 1 ) / ( storeTheta ) ) + 1 )
      {
          numcalls_saved_by_del_by_mu_bound ++;
          goto END;
      }

  #endif // APPLY_DEL_BY_MU_BOUND

  // AuthorK: Del by mu bound
  */


  // AuthorK: Edge bound

  #ifdef APPLY_EDGE_BOUND

  //	printf ( "vertex = %d , local degeneracy  = %d\n" , v , local_degen );
  //	int    S_critical;
  //	int    S1;
  	  int    S          = II -> lb;
      int    P          = II -> itemset.t;
      int    SminusP      = S - P;
    																	//double mu     = II -> frq_lb;                                     //storeTheta = II -> frq_lb;
    																	// compareLeft  = ( ( mu * S * ( S - 1 ) ) / 2 );
    																	// compareRight = sum + ( min_deg * l ) + ( ( l * ( l + 1 ) ) / 2 );

	//AuthorA new: updated if condition which uses degeneracy now: this will, hopefully, make the RHS lower than before
	//and as such will make the inequality tighter
	int g = local_degen - min_deg - 1; //no. of nodes with highest degree that can be added sequentially with current
	//quasi-clique s.t. the min_deg of the resulting subgraph remains < degen

    if ( SminusP > 0 )
    {
        if(g >= 0){ //Apply lemma 9 of our paper since it is tighter than lemma 8
            int max_valueOf_mindeg_S = TWO_MIN(local_degen, min_deg+SminusP);
            int RHS = floor( sum + max_valueOf_mindeg_S*(SminusP + min_deg + 0.5) - 0.5*(min_deg+1)*(2*min_deg+1) ); //( min_deg + 0.5*(g+1) )*g +
            if ( LHS > RHS )
            {
                numcalls_saved_by_edge_bound ++;
#ifdef DEBUG
printf("Pruned by EDGE_BOUND1: LHS = %d, RHS = %d, local_degen = %d, SminusP = %d, g = %d\n", LHS, RHS, local_degen, SminusP, g);
#endif
                goto END;
            }
        }//g>=0
        else{// i.e. if local_degen == min_deg [since c(P) >= delta(P)]
            int RHS = ( int ) floor( sum + SminusP*( min_deg) );
            if( LHS > RHS){
                numcalls_saved_by_edge_bound ++;
#ifdef DEBUG
printf("Pruned by EDGE_BOUND2: LHS = %d, RHS = %d, local_degen = %d, SminusP = %d, g = %d\n", LHS, RHS, local_degen, SminusP, g);
#endif
                goto END;
            }
        }

        //Apply lemma 8 of our paper when lemma 9 is not applicable
        int RHS = ( int ) floor( sum + SminusP*( min_deg + ( SminusP + 1 )*0.5 ) );
    	if ( LHS > RHS )
    	{
        	numcalls_saved_by_edge_bound ++;
#ifdef DEBUG
printf("Pruned by EDGE_BOUND3: LHS = %d, RHS = %d, local_degen = %d, SminusP = %d, g = %d\n", LHS, RHS, local_degen, SminusP, g);
#endif
        	goto END;
    	}
    }
/*
//AuthorA: commenting since this is equivalent to the later checking: min_add_deg > min_deg+1
	else
	{
		int S1 = P + 1;
        if ( ( int ) ceil( ( ( storeTheta * S1 * ( S1 - 1 ) ) / 2 ) - sum ) > ( min_deg + 1 ) )
        {
            numcalls_saved_by_edge_bound ++;
printf("Pruned by EDGE_BOUND4: sum = %d, storeTheta = %f, g = %d\n", sum, storeTheta, g);
            goto END;
        }
	}//else
*/
    //AuthorA new end

    //AuthorA new: commenting out else part as it is useless: the curve of LHS-RHS, now, is a convex upward parabola
    //as such there is no finite point l_critical at which LHS-RHS is maximum
    //it's better to take two extreme values of l and check at which point LHS>RHS, as done by if/else-if above
  /*
    else    //AuthorA
    {
       	double l_critical    = ( storeTheta * P - min_deg - 1 ) / ( 1 - storeTheta );//find global maxima
       	//find integer points just above and below global maxima
       	int l_critical_lb = ceil(l_critical);
       	int l_critical_ub = floor(l_critical+1);

        if( ( l_critical_lb < SminusP ) && ( l_critical_lb > 1 ) ){
            S_critical        = P + l_critical_lb;
            int compareLeft   = ( ( storeTheta * S_critical * ( S_critical - 1 ) ) / 2 );
            int compareRight  = sum + ( min_deg * l_critical_lb ) + ( ( l_critical_lb * ( l_critical_lb + 1 ) ) / 2 );

            if ( compareLeft > compareRight )
            {
                numcalls_saved_by_edge_bound ++;
                goto END;
            }
        }//if

        if( (l_critical_lb != l_critical_ub) && ( l_critical_ub < SminusP ) && ( l_critical_ub > 1 ) ){
            S_critical        = P + l_critical_ub;
            int compareLeft   = ( ( storeTheta * S_critical * ( S_critical - 1 ) ) / 2 );
            int compareRight  = sum + ( min_deg * l_critical_ub ) + ( ( l_critical_ub * ( l_critical_ub + 1 ) ) / 2 );

            if ( compareLeft > compareRight )
            {
                numcalls_saved_by_edge_bound ++;
                goto END;
            }
        }//if
    }//else
  */
  #endif // APPLY_EDGE_BOUND

  // AuthorK: Edge bound


  // AuthorK: Neighbor bound

  #ifdef APPLY_NEIGHBOR_BOUND

    overAllDegree += PP->SG.edge.v[v].t;

    // compareLeft = ( theta * ( ( k_min * ( k_min - 1 ) ) / 2 ) ) - sum - ( ( k_min - k) * ( k_min - k - 1 ) ) / 2 ;
    // compareRight   = overAllDegree;

    if ( ( overAllDegree - sum ) < ( ( ( storeTheta ) * ( ( ( II -> lb ) * ( ( II -> lb ) - 1 ) ) / 2 ) ) - ( ( ( II -> lb ) - ( II -> itemset.t ) ) * ( ( II -> lb ) - ( II -> itemset.t ) - 1 ) ) / 2 ) )
    {
        numcalls_saved_by_neighbor_bound ++;
        goto END;
    }

  #endif // APPLY_NEIGHBOR_BOUND

  // AuthorK: Neighbor bound


  #ifdef APPLY_MAXSIZE_BOUND
    int n = II->itemset.t+1, m = sum+min_deg+1;//max #edges in KU{w} is: |E[K]|+deg_K(v*)+1
    double mu = II->frq_lb;
    int maxSize = (mu+2+sqrt((mu+2)*(mu+2)+8*(m-n)*mu))/(2*mu);
    //printf("current K=%d, n=%d,m=%d,mu=%lf, maxSize=%d\n", II->itemset.t, n,m,mu,maxSize);

    if(II->itemset.t>=maxSize){
            numcalls_saved_by_maxsize_bound++;
            goto END;
    }
  #endif // APPLY_MAXSIZE_BOUND
/*
printf("Itemset: ");fflush(stdout);
QUEUE_print(&II->itemset);fflush(stdout);

printf("lastNode = %d, mindeg = %d, minadddeg = %d\n", v, min_deg, min_add_deg);fflush(stdout);
*/
GEN_CHILD:
  // jump = set of all vertices w that are adjacent to (and therefore can be added to) current pseudo-clique K such that KU{w} is a pseudo-clique
  // i.e. that satisfy the properties of lemma 3 or lemma 4 in the paper
  if ( II -> itemset.t < R )
  {
      jump = PCE_enum_children_clq (PP, prevCh, v);//prevCh, v);
  }//if R-clique couldn't be  found yet
  else {

    if ( min_add_deg > min_deg+1 ) //AUTHOR-A: DON'T CALL IT EDGE-BOUND AS IT IS ALSO DONE BY PCE
    {
      goto END;  //AuthorA: no children possible in this case
    }

    jump = PCE_enum_children (PP, min_deg, min_add_deg);
  }
#ifdef DEBUG
qsort_QUEUE_INT(jump.v, jump.t, 1);
printf("Children: ");fflush(stdout);
QUEUE_print(&jump);fflush(stdout);
#endif
  //the for loop in line 3 of the Algorithm EnumPseudoClique (G = (V, E), K) in paper
  //  int it = 1;

#ifdef APPLY_SIBLING_BOUND
    //printf("\nIn sibling bound\n");
    unsigned noSubTreeUnderPrevSibling = 0;//assume false first so that the if inside loop isn't entered the first time
    int prevSibling = -1;
  MQUE_BLOOP (jump, x){//go backward direction so that the nodes outside P having higher deg_P comes first in the loop
    noSubtreeUnderCurTreenode = 0;//there exist some children of the current Treenode
    //if the a vertex u exist such that no subtree exist under PU{u} then no subtree will exist under PU{v} for all v with deg_P(u) > deg_P(v)
    //since we are going backward in this loop, deg_P(u) > deg_P(v) for all v that comes next in the loop from the jump QUEUE
    //so we need not make the calls for PU{v} for all such v if PU{v} need not be reported (i.e. |PU{v}| < II->lb)
    if(noSubTreeUnderPrevSibling && II->itemset.t < (II->lb - 1) && PP->occ.list[prevSibling] > PP->occ.list[*x]){
        numcalls_saved_by_sibling_bound++;
        //printf("\nPruned by SIBLING bound: "); QUEUE_print (&(II->itemset));
        break;
    }//if
    noSubTreeUnderPrevSibling = PCE_iter (PP, *x, sum, all, local_degen); // recursive call; analogous to EnumPseudoClique(G, K U {w}) where w is a potential vertex which can be added with K
    prevSibling = *x;
#else
  MQUE_FLOOP (jump, x){
    PCE_iter (PP, *x, sum, all, local_degen, &jump); // recursive call; analogous to EnumPseudoClique(G, K U {w}) where w is a potential vertex which can be added with K
#endif // SIBLING_BOUND
  //printf("iteration#%d end\n", it);
  //it++;
  }//end of loop
  QUEUE_end (&jump);

  END:;
  PCE_rm_vertex_from_clique (PP, v, &sum, &all);  // remove vertex v from the current clique K

  // AuthorK: Neighbor bound

  #ifdef APPLY_NEIGHBOR_BOUND
    overAllDegree -= PP->SG.edge.v[v].t;
  #endif // APPLY_NEIGHBOR_BOUND

  // AuthorK: Neighbor bound


  #ifdef APPLY_SIBLING_BOUND
    return noSubtreeUnderCurTreenode;
  #endif // SIBLING_BOUND
}


/****************/
/* main routine */
/****************/
int PCE_main (int argc, char *argv[]){
  QUEUE_INT i;
  PROBLEM PP;
  SGRAPH *G = &PP.SG;

// QUEUE PP->itemcand;   /* queue for candidates */
// PP->itemary:  #vertices of the current clique adjacent to each vertex */
// PP->tmp: used as a temporary array, for itemset
// PP->itemcand: used to keep the list of candidates, temporary (solid one is stored in local variable "jump"

PROBLEM_init (&PP);
  PCE_read_param (argc, argv, &PP);
if ( ERROR_MES ) return (1);
  G->flag |= LOAD_PERM + LOAD_INCSORT + LOAD_RM_DUP + LOAD_EDGE;
  PP.II.flag |= ITEMSET_ADD;
  PROBLEM_load (&PP);

  PROBLEM_alloc (&PP, G->edge.t, G->edge.t, G->edge.eles, NULL, PROBLEM_ITEMARY +PROBLEM_VECARY +PROBLEM_ITEMCAND);
  MALIST_alloc (&PP.occ, G->edge.t, G->edge.t+1); // No need to call MALIST_end to clear this memory; as it's called in PROBLEM_end(...)
  FLOOP (i, 0, G->edge.t) MALIST_ins_tail (&PP.occ, 0, i, 0);

  //SGRAPH_init_hash(G);
  SGRAPH_init_bmp(G);

  int *R_1coreNodes, R_1coreNodes_endIndex = -1;


  //AuthorA new end

  if ( !ERROR_MES ){
#ifdef APPLY_STRICT_MINDEG_BOUND
    PP.II.itemset.prevMinDeg = 0;//AuthorA: init prevMinDeg
#endif // APPLY_STRICT_MINDEG_BOUND

        ITEMSET *II  = &PP.II;

        int v, i, j, min_deg;
        double sum = 0.0, all = 0.0;
        QUEUE_INT mindeg_vertex;

//		degen_cutoff = ( int ) ceil ( ( ( II -> frq_lb ) / 2 ) * ( ( II -> lb ) - 0.5 ) );

		storeTheta   = II -> frq_lb;
		//R			       = ( int ) floor ( ( 1 / ( 1 - ( ( storeTheta * ( ( II -> lb ) - 1 ) ) / ( II -> lb ) ) ) ) + 1 );    // AuthorK: Previous turan filtering equation
		R			       = ( int ) ceil  ( 1 / ( 1 - ( ( storeTheta * ( ( II -> lb ) - 1 ) ) / ( II -> lb ) ) ) );

    printf ( "R                   = %d\n" , R );fflush(stdout);

		#ifdef APPLY_TURAN_FILTERING

//		if ( ( R - 1 ) > degen_cutoff )
//			degen_cutoff = R - 1;

		#endif // APPLY_TURAN_FILTERING

//		printf ( "degen_cutoff        = %d\n" , degen_cutoff );
		//AuthorA start: coreStatus calculation via 3rd party and reading from the generated file
    malloc2 ( coreStatus , G -> edge.t , EXIT );
    unsigned maxdeg = getcores(G, coreStatus);
   //int arrayIndex = 0;
    //int coreValue;

    int minValue   = INF;
    int maxValue   = -1;

    for ( int i = 0 ; i < G -> edge.t ; i += 1 )
    {
        minValue = TWO_MIN ( coreStatus [ i ] , minValue );
        maxValue = TWO_MAX ( coreStatus [ i ] , maxValue );
    }

    //for ( int i = 0 ; i < G -> edge.t ; i += 1 )
    //  printf ( "%d\n" , coreStatus [ i ] );

    int min_degree;

    min_degree = min_core = minValue;
    degen      = k_core   = maxValue;

//    printf ( "minimum degree      = %d\n" , min_degree );
    printf ( "minimum core        = %d\n" , min_core );
//    printf ( "K-core              = %d\n" , k_core );
    printf ( "degeneracy = max. core = %d\n" , degen );

    //AuthorA end: coreStatus calculation+reading


/*

///AuthorA: commenting out previous code -- which uses eligible_roots; we aren't using that anymore since
//using R_1corenodes (instead of eligible_roots) suffices our purpose and never the condition: degen_cutoff > R-1
//comes true, in our experience

		for ( i = 0 ; i < G -> edge.t ; i += 1 ){
			if ( coreStatus [ i ] >= degen_cutoff )
				eligible_roots [ ++eligible_roots_endIndex ] = i;
		}

    if(degen_cutoff > R-1){///todo:find a scenerio when this happens (so far couldn't find)
        malloc2(R_1coreNodes, G->edge.t,EXIT);

        for ( i = 0 ; i < G -> edge.t ; i += 1 ){
            if ( coreStatus [ i ] >= R-1 )
                R_1coreNodes [ ++R_1coreNodes_endIndex ] = i;
        }
    }
    else{
        R_1coreNodes = eligible_roots;
        R_1coreNodes_endIndex = eligible_roots_endIndex;
    }
*/

    malloc2(R_1coreNodes, G->edge.t,EXIT);

    for ( i = 0 ; i < G -> edge.t ; i += 1 ){
        if ( coreStatus [ i ] >= R-1 )
            R_1coreNodes [ ++R_1coreNodes_endIndex ] = i;
    }

		//printf ( "eligible roots = " );
		//for ( i = 0 ; i <= eligible_roots_endIndex ; i += 1 )
		//	printf ( "%d " , eligible_roots [ i ] );
		//printf ( "\n" );

        //calculate upper-bound
        int ub_from_degen1 = II->ub = G->edge.t;//set default value of ub
        double numerator = (degen+1)*II->frq_lb;
        double denominator = numerator - degen;
        if(denominator > 0) ///todo: check why denominator is sometimes < 0. What does it mean?
            ub_from_degen1 = floor(numerator/denominator);//calculated using Turan's theorem
        int ub_from_degen2 = floor(((2.0*degen)/II->frq_lb) + 0.5);//found in lemma 4(iii) of (Komusiewicz, 2015)
        int new_ub = (ub_from_degen1 < ub_from_degen2) ? ub_from_degen1 : ub_from_degen2;
        //if(II->ub > new_ub)
        //    II->ub = new_ub;//ub = max of a quasi-clique in this graph

        II -> ub = TWO_MIN ( II -> ub , new_ub );

        printf("Upper-bound, II->ub = %d\n", II->ub);

        if(new_ub < II->lb){
            printf("No quasi-clique of size at least %d can exist in this graph", II->lb);
            goto END_MAIN;
        }
        printf ( "\n" );


/*
// AuthorK: Local degeneracy bound

#ifdef APPLY_LOCAL_DEGEN_BOUND

	//int S 			     = II->lb;
	LOCAL_DEGEN_CO_NUMER = ( II->frq_lb * ( II->lb ) * ( ( II->lb ) - 1 ) ) / 2.0;

#endif // APPLY_LOCAL_DEGEN_BOUND

// AuthorK: Local degeneracy bound
*/


#ifdef APPLY_EDGE_BOUND
    //int S = II->lb;
	LHS   = ceil( ( II->frq_lb * ( II->lb ) * ( ( II->lb ) - 1 ) ) / 2.0 ); //we want to avoid recalculating LHS in different conditions
	//and in different calls of PCE_iter; so calculating it here once and for all
#endif //APPLY_EDGE_BOUND
//AuthorA new end


	// AuthorK: Turan filtering

    //ITEMSET *II = &PP.II;
    //storeTheta  = II -> frq_lb;

    //R			= ( int ) ceil  ( ( 1 / ( 1 - ( ( storeTheta * ( ( II -> lb ) - 1 ) ) / ( II -> lb ) ) ) ) );

   	//R			= ( int ) floor ( ( 1 / ( 1 - ( ( storeTheta * ( ( II -> lb ) - 1 ) ) / ( II -> lb ) ) ) ) + 1 );  		// AuthorK: Previous turan filtering equation

	  //printf ( "R = %d\n" , R );

    // AuthorK: Turan filtering


    //AuthorA new start

    //updated code so that PCE follows reverse-degeneracy order of vertices (in the whole graph)
    //while exploring the reverse-search tree (at least in the level-0 of the reverse-search tree).
    //This is just a heuristic that is helpful to find quasi-cliques early during execution
    //but it seems to have no effect on the execution time itself. If we can derive some inequality
    //to prune the remaining search tree after we find all the quasi-cliques, then this heuristic
    //may become useful for saving time. But I couldn't figure out any idea to know beforehand
    //that all the quasi-cliques are found and we can now stop.

    ///todo: find an upper-bound that may help us to achieve this goal

    //FLOOP (i, mindegVertices_startIndex, G->edge.t) //PCE_iter (&PP, i, 0, 0);//AuthorA new : commented out here

    //printf ( "eligible roots endIndex = %d\n", eligible_roots_endIndex );

//    if(degen_cutoff > R-1)
        getSubgraph(G, R_1coreNodes, R_1coreNodes_endIndex+1,maxdeg, &PP.subg);
//    else
//        getSubgraph(G, eligible_roots, eligible_roots_endIndex+1,maxdeg, &PP.subg);

    for(i = R_1coreNodes_endIndex; i >= 0; i--)
        PCE_iter (&PP, R_1coreNodes[i], 0, 0, coreStatus [ R_1coreNodes[i] ], NULL);
//    printf ( "here\n");fflush(stdout);
    //FLOOP (i, 0, G->edge.t) PCE_iter (&PP, i, 0, 0, coreStatus [ i ]);

    //AuthorA new end
    ITEMSET_last_output (&PP.II);
/*
  typedef struct{
    int u,v;
  }Edge;

  int numNodes = R_1coreNodes_endIndex+1;
  unsigned int numEdges = 0;
  Edge *edgelist;
  malloc2(edgelist, PP.subg->edge.eles/2, EXIT);
  for(i = R_1coreNodes_endIndex; i >= 0; i--){
    int u = R_1coreNodes[i];
    for(j=0; j < PP.subg->edge.v[u].t; j++) {
      int v = PP.subg->edge.v[u].v[j];
      Edge e;
      if(u > v){
        e.u = u;//parent
        e.v = v;//child in the reverse-search tree
        edgelist[numEdges++] = e;
      }
    }
  }

  printf("#nodes = %d, #edges = %u\n", numNodes, numEdges);

  int k, *z;

  if(R >= 3){
    QUEUE Q[THREAD];
    for(k=0;k<THREAD; k++){
      QUEUE_alloc (&Q[k], G->edge.t+2);
    }

    for (k = 0; k < numEdges; k++){
      double sum = 0, all = 0;
      unsigned int u = edgelist[k].u, v = edgelist[k].v;
//      printf("(%d,%d)\n", i, j);
      PCE_add_vertex_to_clique (&PP_list [ k ], u, &sum, &all);  // add vertex to the current clique
      PCE_add_vertex_to_clique (&PP_list [ k ], v, &sum, &all);  // add vertex to the current clique
      PP_list [ k ].II.frq = 1.0;//sum/all;

      Q[k].t = 0;
      for(j=0; j < PP_list [ k ].subg->edge.v[u].t; j++) {
        int nbr = PP_list [ k ].subg->edge.v[u].v[j];
        if(nbr < v && SGRAPH_is_edge_bmp(nbr, v))//ensure that {u,v,nbr} forms a triangle since R >= 3. NOTE: II->itemset.t == 2 and prevCommonNbrs = NULL can only happen when R>=3, as per rules for forming Edge_list in PCE_main
          QUE_INS (Q[k], nbr);
      }
      MQUE_FLOOP(Q[k],z)
        PCE_iter (&PP_list [ k ], *z, sum, all, TWO_MIN(coreStatus [u], coreStatus[v]), &Q[k]);

      PCE_rm_vertex_from_clique (&PP_list [ k ], v, &sum, &all);  // remove vertex from the current clique
      PCE_rm_vertex_from_clique (&PP_list [ k ], u, &sum, &all);  // remove vertex from the current clique
    }

    for(k=0;k<THREAD; k++){
      QUEUE_end (&Q[k]);
    }
  }//if R-1 >= 2
  else if(R == 2){
    for (k = 0; k < numEdges; k++){
	   	double sum = 0, all = 0;
	   	unsigned int u = edgelist[k].u, v = edgelist[k].v;
      //printf("(%d,%d)\n", i, j);
  		PCE_add_vertex_to_clique (&PP_list [ k ], u, &sum, &all);  // add vertex v to the current clique
	    PCE_iter (&PP_list [ k ], v, sum, all, coreStatus[u], NULL);//Author1: In this case, using prevCh=NULL is OK because ( II -> itemset.t < R ) condition checking in this PCE_iter call will fail and as such we shall not call PCE_enum_children_clq from there
      PCE_rm_vertex_from_clique (&PP_list [ k ], u, &sum, &all);  // remove vertex v from the current clique
    }
  }
  else{
    print_err("Error: cannot enumerate dense subgraphs when R < 2.");
    EXIT;
  }

  //Author-A: the following function simply outputs some stat about itemsets in console
  //so no problem if we don't call it
  for ( int i = 0 ; i < THREAD ; i += 1 )
    ITEMSET_last_output (&PP_list [ i ].II);

  free2(edgelist);
*/
}//if not ERR MSG

  // AuthorK: Connected Clique
  #ifdef CONNECTED_CLIQUE
  printf ( "CONNECTED_CLIQUE\n" );
  #endif // CONNECTED_CLIQUE
  // AuthorK: Connected Clique

  // AuthorK: Turan Filtering
  #ifdef APPLY_TURAN_FILTERING
  printf ( "APPLY_TURAN_FILTERING\n" );
  #endif // APPLY_TURAN_FILTERING
  // AuthorK: Turan Filtering


  // AuthorK: #calls of PCE iteration save

  #ifdef APPLY_SIBLING_BOUND
  printf("#calls of PCE_iter saved by SIBLING bound = %llu\n" , numcalls_saved_by_sibling_bound );
  #endif // APPLY_SIBLING_BOUND

  #ifdef APPLY_EDGE_BOUND
  printf("#calls of PCE_iter saved by EDGE bound = %llu\n" , numcalls_saved_by_edge_bound );
  #endif // APPLY_EDGE_BOUND

 #ifdef APPLY_ORDER_BOUND
  printf("#calls of PCE_iter saved by ORDER bound = %llu\n" , numcalls_saved_by_order_bound );
#endif // APPLY_ORDER_BOUND

  /*
  #ifdef APPLY_LOCAL_DEGEN_BOUND
  printf ( "#calls of PCE_iter saved by LOCAL_DEGEN_BOUND bound = %llu\n" , numcalls_saved_by_local_degen_bound );
  #endif // APPLY_LOCAL_DEGEN_BOUND
  */

  #ifdef APPLY_NEIGHBOR_BOUND
  printf("#calls of PCE_iter saved by NEIGHBOR_BOUND bound = %llu\n" , numcalls_saved_by_neighbor_bound );
  #endif // APPLY_NEIGHBOR_BOUND

  /*
  // AuthorK: H-Index bound permanently disable
  #ifdef APPLY_H_INDEX_BOUND
  printf("#calls of PCE_iter saved by H_INDEX_BOUND bound = %llu\n" , numcalls_saved_by_h_index_bound );
  #endif // APPLY_H_INDEX_BOUND
  */

  #ifdef APPLY_MAXSIZE_BOUND
  printf("#calls of PCE_iter saved by MAXSIZE bound = %llu\n", numcalls_saved_by_maxsize_bound);
  #endif // APPLY_MAXSIZE_BOUND

  #ifdef APPLY_MINDEG_BOUND
  printf("#calls of PCE_iter saved by MINDEG bound = %llu\n", numcalls_saved_by_mindeg_bound);
  #endif // APPLY_MINDEG_BOUND

  #ifdef APPLY_MAXDEG_BOUND
  printf("#calls of PCE_iter saved by MAXDEG bound = %llu\n", numcalls_saved_by_maxdeg_bound);
  #endif // APPLY_MAXDEG_BOUND

  // AuthorK: #calls of PCE iteration save


//  printf("\nTime taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
  printf("#calls of PCE_iter = %llu\n", numCallPCE_Iter); //AuthorA

    freeSubgraph(&PP.subg, R_1coreNodes, R_1coreNodes_endIndex+1);
//AuthorA new cleanup
END_MAIN:;
    //free2(mindegVals);
    //free2(mindegVertices);
//    if(degen_cutoff > R-1)


//    if(R_1coreNodes != eligible_roots)
        free2(R_1coreNodes);

//    free2(eligible_roots);

    //free2(vertex2index);
//AuthorA new cleanup ends

	free2 ( coreStatus );

  PROBLEM_end (&PP);
  //SGRAPH_free_hash();//AuthorA
  SGRAPH_free_bmp(G);
  printf("freed everything ...");fflush(stdout);
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
