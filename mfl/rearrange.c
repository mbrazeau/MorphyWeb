/*
 *  rearrange.c
 *  Morphy
 *
 *  Functions useful to tree rearrangements
 *
 */

#include "morphy.h"

void mfl_bswap(node *p, node *q)
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

struct node * mfl_remove_branch(node *n)
{
    node *p, *q, *nb;
    
    nb = n->outedge;
    
    p = nb->next->outedge;
    q = nb->next->next->outedge;
    nb->next->outedge = NULL;
    nb->next->next->outedge = NULL;
    joinNodes(p, q);
    
    nb = nb->next;
    return nb;
}