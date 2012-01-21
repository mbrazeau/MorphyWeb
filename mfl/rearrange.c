/*
 *  rearrange.c
 *  Morphy
 *
 *  Functions useful to tree rearrangements
 *
 */

#include "morphy.h"
#define TREELIMIT 500

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

void mfl_insert_branch(node *br, node *target)
{
    // Inserts a branch with a ring base into another branch
    
    node *p, *bout, *tdesc;
    
    tdesc = target->outedge;
    
    // Find an available node in the ring 
    p = br->next;
    while (p != br) {
        if (!p->outedge) {
            bout = p;
            p = br;
        }
        else {
            p = p->next;
        }
    }
    
    joinNodes(br, target);
    joinNodes(bout, tdesc);
}

void mfl_nni_traversal(node *n, tree *swapingon, tree **treeset, int ntax, int numnodes, long int *current, bool *undertreelimit)
{
    node *p;
    
    if (n->tip || !(*undertreelimit)) {
        return;
    }
    
    if (!n->outedge->tip)
    {
        if ((*current < TREELIMIT) && ((*current + 2) < TREELIMIT)) {
            mfl_bswap(n->next->outedge, n->outedge->next->outedge);
            printNewick(swapingon->trnodes[0]);
            printf("\n");
            treeset[*current] = copytree(swapingon, ntax, numnodes);
            treeset[*current]->index = *current;
            *current = *current + 1;
            mfl_bswap(n->next->outedge, n->outedge->next->outedge);
            printf("current: %li\n", *current);
            
            mfl_bswap(n->next->next->outedge, n->outedge->next->outedge);
            printNewick(swapingon->trnodes[0]);
            printf("\n");
            treeset[*current] = copytree(swapingon, ntax, numnodes);
            treeset[*current]->index = *current;
            *current = *current + 1;
            mfl_bswap(n->next->next->outedge, n->outedge->next->outedge);
            printf("current: %li\n", *current);
        }
        else {
            *undertreelimit = false;
            return;
        }

        //return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_nni_traversal(p->outedge, swapingon, treeset, ntax, numnodes, current, undertreelimit);
        p = p->next;
    }
}

void test_nni(int ntax, int numnodes)
{
    long int i, travlimit = 1;
    long int branch_count = 0;   // The branch iteration of the NNI traversal
    long int *bc_pointer = &branch_count;
    long int numreps = TREELIMIT;
    tree **treeset;
    bool undertreelimit = true, *undertreelimit_p = &undertreelimit;
    
    treeset = (tree**) malloc(TREELIMIT * sizeof(tree*));
    
    treeset[0] = randunrooted(ntax, numnodes);
    treeset[0]->index = 0;

    for (i = 0; i < 10 /*numreps && undertreelimit*/; ++i) {
        mfl_nni_traversal(treeset[i]->trnodes[0]->outedge, treeset[i], treeset, ntax, numnodes, bc_pointer, undertreelimit_p);
    }
    
    for (i = 0; i < *bc_pointer; ++i) {
        printf("Tree %i:\n", i + 1);
        printNewick(treeset[i]->trnodes[0]);
        printf(";\n");
    }
    
    freetree(treeset[0], numnodes);
    free(treeset);
}