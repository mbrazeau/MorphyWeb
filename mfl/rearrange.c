/*
 *  rearrange.c
 *  Morphy
 *
 *  Functions useful to tree rearrangements
 *
 */

#include "morphy.h"

extern int ntax;
extern int numnodes;

void bswap(node *p, node *q)
{
    /* Exchanges the position of two nodes */
    
    node *pbase, *qbase;
    
    pbase = p->outedge;
    qbase = q->outedge;
    
    p->outedge = qbase;
    q->outedge = pbase;
    pbase->outedge = q;
    qbase->outedge = p;
    
}

