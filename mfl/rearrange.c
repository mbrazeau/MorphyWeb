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

void mfl_nni_traversal(node *n, tree *swapingon, tree **treeset, int ntax, int numnodes, int *current)
{
    node *p;
    
    if (n->tip) {
        return;
    }
    
    if (!n->outedge->tip)
    {
        mfl_bswap(n->next->outedge, n->outedge->next->outedge);
        //printNewick(swapingon->trnodes[0]);
        //printf("\n");
        treeset[*current] = copytree(swapingon, ntax, numnodes);
        treeset[*current]->index = *current;
        *current = *current + 1;
        mfl_bswap(n->next->outedge, n->outedge->next->outedge);


        mfl_bswap(n->next->next->outedge, n->outedge->next->outedge);
        //printNewick(swapingon->trnodes[0]);
        //printf("\n");
        treeset[*current] = copytree(swapingon, ntax, numnodes);
        treeset[*current]->index = *current;
        *current = *current + 1;
        mfl_bswap(n->next->next->outedge, n->outedge->next->outedge);
        //return;
    }
    
    p = n->next;
    while (p != n) {
        mfl_nni_traversal(p->outedge, swapingon, treeset, ntax, numnodes, current);
        p = p->next;
    }
}

void test_nni(int ntax, int numnodes)
{
    
    int counter = 1;
    int *cptr = &counter;
    tree **treeset;
    
    treeset = (tree**) malloc(TREELIMIT * sizeof(tree*));
    
    treeset[0] = randunrooted(ntax, numnodes);
    treeset[0]->index = 0;
    //dump_tree(treeset[0], ntax, numnodes);
    printNewick(treeset[0]->trnodes[0]);
    printf("\nin test_nni\n");
    mfl_nni_traversal(treeset[0]->trnodes[0]->outedge, treeset[0], treeset, ntax, numnodes, cptr);
    printNewick(treeset[0]->trnodes[0]);
    printf(";\n");
    //dump_tree(treeset[0], ntax, numnodes);
    printNewick(treeset[1]->trnodes[0]);
    printf(";\n");
    //dump_tree(treeset[1], ntax, numnodes);
    printNewick(treeset[2]->trnodes[0]);
    printf(";\n");
    //dump_tree(treeset[2], ntax, numnodes);
    printNewick(treeset[3]->trnodes[0]);
    printf(";\n");
    //dump_tree(treeset[3], ntax, numnodes);
    printNewick(treeset[4]->trnodes[0]);
    printf(";\n");
    //dump_tree(treeset[4], ntax, numnodes);

    freetree(treeset[0], numnodes);
    free(treeset);
}