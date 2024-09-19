#include "sgraph.c"
#include "stdlib2.h"
#include "vec.h"

#define uli unsigned long int
typedef struct{
uli *bits;
uli sz;
unsigned maxID;
}bitmap;

bitmap bmpSubgNodes;   //, **bmpSubgEdges;
//unsigned *nodeID2Index;
//unsigned *subgNodes;
//unsigned numSubgNodes;

unsigned int bitsInWord = sizeof(uli)*8;
SGRAPH INIT_SGRAPH2 = {TYPE_SGRAPH,NULL,0,INIT_SETFAMILY_,INIT_SETFAMILY_,INIT_SETFAMILY_,0,NULL,NULL,NULL};

//#define DEBUG 1

void insertBit(bitmap *bmp, uli id){
    if(id > bmp->maxID){
        printf("Error: insertion of id>maxID impossible");
        exit(1);
    }
    uli k = 1;
    uli bmp_index = id/bitsInWord;
    uli bmp_val = k << ( id%bitsInWord );

    bmp->bits[bmp_index] = bmp->bits[bmp_index] | bmp_val;
}

unsigned isPresent(bitmap *bmp, uli id){
    if(id > bmp->maxID) return 0;

    uli k = 1;
    uli bmp_index = id / bitsInWord;
    uli bmp_val = k << ( id % bitsInWord );

    return ( (bmp->bits[bmp_index] & bmp_val) != 0 );
}
/*
//AuthorA: incomplete/buggy: DON'T USE
unsigned isPresentInRange(bitmap *bmp, uli id, uli idmax){
    if(id >= idmax) return 0;

    uli k = 1, x = -1;
    uli mask = x >> (bitsInWord - ( idmax % bitsInWord ) + 1);//keep only leftmost (idmax-1) bits
    uli bmp_index = id / bitsInWord;
    uli bmp_val = k << ( id % bitsInWord );

    return ( (bmp->bits[bmp_index] & mask & bmp_val) != 0 );
}
*/

void createBitMap(bitmap* bmp, const int nodes[], unsigned int n)
{
    if(bmp == NULL){
        printf("Error - failed to allocate...");
        exit(1);
    }
    unsigned int i;
    bmp->maxID = 0;
    for (i = 0; i < n; i++) {
        if(nodes[i] > bmp->maxID)
            bmp->maxID = nodes[i];
    }

    double x = bmp->maxID*1.0/bitsInWord;
    bmp->sz = (x - (uli)x > 0)? (uli)x+1 : (uli)x;//set bmp->sz = ceil(maxID*1.0/bitsInWord)
    malloc2 ( bmp->bits, bmp->sz, EXIT );
#ifdef DEBUG
    printf("%d, %d, %lf %lu\n", n, bitsInWord, x, bmp->sz);
#endif // DEBUG
    for (i = 0; i < bmp->sz; i++)
        bmp->bits[i] = 0;

#ifdef DEBUG
    for (i = 0; i < bmp->sz; i++)
        printf("%lu ", bmp->bits[i]);
#endif // DEBUG
//printf("\nInserting:\n");fflush(stdout);

    for (i = 0; i < n; i++) {
//printf("%d, ", nodes[i]);
        insertBit(bmp, nodes[i]);
    }
//printf("\n");
#ifdef DEBUG
    printf("\n");
    for (i = 0; i < bmp->sz; i++)
        printf("%lu ", bmp->bits[i]);
    printf("\n");fflush(stdout);
#endif // DEBUG
}

void deallocateBitmap(bitmap *bmp){
//    if(bmp != NULL){
        if(bmp->bits != NULL){
            free2 ( bmp->bits );
            bmp->bits = NULL;
        }
//    }
}

void getSubgraph(const SGRAPH *g, const int nodes[], unsigned int n,
                 unsigned int md, SGRAPH **sub)
{
    unsigned int tn = g->edge.t;

    uli m =  md*n/2;//max. no. of edges in a n-node graph with max. deg = md
///todo: use a tighter bound on max. edges later
    unsigned i, j;

    createBitMap(&bmpSubgNodes, nodes, n);
//printf("n=%u, md=%u, tn=%u\n", n, md, tn);fflush(stdout);
    *sub = &INIT_SGRAPH2;
//
//
    //SGRAPH_alloc (sub, n, m, 0);
    SETFAMILY_alloc (&((*sub)->edge), tn, NULL, tn, m);
    (*sub)->flag = g->flag;
    (*sub)->edge.t = n;
    //SETFAMILY_alloc (&G->edge, nodes, NULL, nodes, edge_num);

    for (i = 0; i < n; i++) {
        malloc2 ( (*sub)->edge.v[nodes[i]].v, g->edge.v[nodes[i]].t, EXIT );
    }

//    printf ( "#edges(G) = %ld, #edges(*sub) = %ld\n" , g->edge.eles/2, (*sub)->edge.eles/2 );
//    fflush(stdout);
    for (i = 0; i < n; i++) {
        //if(isPresent(nodes[i])){
            //malloc2 ( (*sub) ->edge.v[ nodes[i] ].v, md, EXIT );
//            printf("%u: ", nodes[i]);fflush(stdout);
            for (int j=0; j < g->edge.v[nodes[i]].t ; j++){
                unsigned int nbr = g->edge.v[nodes[i]].v[j];
                if( nodes[i]<nbr && isPresent(&bmpSubgNodes, nbr) ){//if nbr is also present in nodes array
//                    printf("%u, ", nbr);fflush(stdout);
                    SGRAPH_edge_mk( (*sub), nodes[i], nbr, 0);//create edge in subgraph *sub
                }//if
            }
        //}//if

 //       printf("\n");fflush(stdout);
    }

    printf ( "#edges = %ld, #edges(*sub) = %ld\n" , g->edge.eles/2, (*sub)->edge.eles/2 );
    fflush(stdout);

//    printf ( "&edges = %x, %x\n" , g->edge.v, (*sub)->edge.v);
//    fflush(stdout);

}

void freeSubgraph(SGRAPH **sub, const int nodes[], unsigned int n)
{
    unsigned i;
    for (i = 0; i < n; i++) {
        free2 ( (*sub)->edge.v[nodes[i]].v );
    }

    SETFAMILY_end (&((*sub)->edge));
    deallocateBitmap(&bmpSubgNodes);
}
