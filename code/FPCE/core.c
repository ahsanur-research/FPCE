#include "sgraph.h"
#include "stdlib2.h"

unsigned getcores(SGRAPH *g, unsigned *core) {
/*    FILE *fp    = fopen   ( kcoreString , "w" );
    if(fp == NULL){
        printf("error opening file: %s\n", kcoreString);
        return;
    }

*/
    unsigned int n = g->edge.t, d, md = 0, i, start, num;
    unsigned int v, u, w, du, pu, pw;
    //unsigned int vert[n], pos[n], deg[n]; //, core[n];
    unsigned int *vert , *pos , *deg;
    malloc2 ( vert , g -> edge.t , EXIT );
    malloc2 ( pos  , g -> edge.t , EXIT );
    malloc2 ( deg  , g -> edge.t , EXIT );
    //printf("hello\n");fflush(stdout);
    for (i = 0; i < n; i++) {
        deg[i] = g->edge.v[i].t;
        //printf("deg[%u] = %u, neighbors: ", i, deg[i]);fflush(stdout);

/*        printf ( "%d => " , i );
        for (int j=0; j<g->edge.v[i].t ; j++){
            int e = g->edge.v[i].v[j];
            printf ("%d ", e);
        }
        printf ( "\n" );
*/
        if (deg[i] > md)
            md = deg[i];
    }
    //printf("max deg = %u\n", md);fflush(stdout);
    unsigned int bin[md+1];
    //now md = max. deg. in the input graph
    for (d = 0; d <= md; d++) bin[d] = 0;
    for (v = 0; v < n; v++)
        bin[deg[v]] ++;
    //now bin[d] is the number of vertices with degree d
    start = 0;
    for (d = 0; d <= md; d++) {
        /* printf ( "bin [ %d ] = %d\n" , d , bin [ d ] ); */
        num = bin[d];
        bin[d] = start;
        start += num;
        /* printf ( "bin [ %d ] = %d\n\n" , d , bin [ d ] ); */
    }
    //now bin[d] = sum_{i=0 to d-1} (number of vertices with degree i)
    for (v = 0; v < n; v++) {
        pos[v] = bin[deg[v]];
        vert[pos[v]] = v;
        bin[deg[v]]++;
    }

/*     printf ( "\n" );

    for (v = 0; v < n; v++)
        printf ( "vert [ %d ] = %d\n" , v , vert [ v ] );

    printf ( "\n" );

    for (d = 0; d <= md; d++)
        printf ( "bin [ %d ] = %d\n" , d , bin [ d ] ); */

    //now bin[d] = sum_{i=0 to d} (number of vertices with degree i),
    //vert array contains the vertex ids in the increasing order of their degrees,
    //and pos[v] is the index of vertex v in the vert array
    for (d = md; d >= 1; d--)
        bin[d] = bin[d-1];
    //now bin[d] = sum_{i=0 to d-1} (number of vertices with degree i), again
    bin[0] = 0;

    /* printf ( "\n\n" ); */

/*     for (d = 0; d <= md; d++)
        printf ( "bin [ %d ] = %d\n" , d , bin [ d ] ); */

    /* printf ( "\n\n" ); */

    for (i = 0; i < n; i++) {
/*         for ( int kk = 0 ; kk < n ; kk += 1 )
            printf ( "deg  [ %d ] = %d; " , kk , deg [kk ] );
        printf ( "\n" );

        for ( int kk = 0 ; kk < n ; kk += 1 )
            printf ( "vert [ %d ] = %d; " , kk , vert [kk ] );
        printf ( "\n" ); */

        v = vert[i];
        core[v] = deg[v];
        /* printf ( "core [ %d ] = %d\n" , v , core [ v ] ); */
        //find those neighbors of v which have degree > deg[v] in order to update (decrement) their degrees
        for (unsigned j=0; j < g->edge.v[v].t ; j++){//G->edge.v[i] is the data structure of i-th node
            unsigned u = g->edge.v[v].v[j];//G->edge.v[i].v is the adj. list of i-th node
            if (deg[u] > deg[v]) {
                //in this case, we need to decrement deg[u]
                //to do that we find the first vertex w in u's bin i.e.
                //w is the first vertex in vert that has degree = deg[u]
                du = deg[u];
                pu = pos[u];
                pw = bin[du];
                w = vert[pw];
                //if u!=w then we swap u and w in vert array
                //as well as swap their indices in pos array

                /* printf ( "u = %d; w = %d\n" , u , w ); */

                if (u != w) {
                    pos[u] = pw;
                    vert[pu] = w;
                    pos[w] = pu;
                    vert[pw] = u;
                }

                bin[du]++;//update bin array accordingly
                deg[u]--;
            }
        }//foreach neighbor u
        /* printf ( "\n" ); */
    }//foreach node v

/*    for (v = 0; v < n; v++)
        fprintf(fp, "%u\n", core[v]);
    fclose(fp);
*/
    free2 ( vert );
    free2 ( pos  );
    free2 ( deg  );

//    printf ( "#edges = %ld\n" , g->edge.eles/2 );

    return md;

}
