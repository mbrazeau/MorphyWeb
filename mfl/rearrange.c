/*
 *  rearrange.c
 *  Morphy
 *
 *  Functions useful to tree rearrangements
 *
 */

#include "morphy.h"
#define TREELIMIT 450

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

void mfl_nni_traversal(node *n, tree *swapingon, tree **treeset, int ntax, 
                       int nchar, int numnodes, long int *current, 
                       charstate *tipdata, bool *undertreelimit, 
                       long int *currentbesttree, bool *foundbettertree)
{
    int trlength = 0;
    node *p;
    
    if (n->tip || !(*undertreelimit)) {
        return;
    }
    
    if (!n->outedge->tip)
    {
        /* Nearly identical steps are conducted twice because all nearest-
         * neighbor interchanges produce two distinct tree topologies */
        
        if ((*current + 1 <= TREELIMIT) && ((*current + 2) <= TREELIMIT)) {
            mfl_bswap(n->next->outedge, n->outedge->next->outedge);
            mfl_root_tree(swapingon, 1, ntax);
            trlength = mfl_get_treelen(swapingon, tipdata, ntax, nchar);
            unroot(ntax, swapingon);
            if (trlength <= *currentbesttree) {
                *foundbettertree = true;
                *currentbesttree = trlength;
                treeset[*current] = copytree(swapingon, ntax, numnodes);
                treeset[*current]->index = *current;
                treeset[*current]->length = trlength;
                trlength = 0;
                *current = *current + 1;
            }
            mfl_bswap(n->next->outedge, n->outedge->next->outedge);
            //printf("current: %li\n", *current);
            
            mfl_bswap(n->next->next->outedge, n->outedge->next->outedge);
            mfl_root_tree(swapingon, 1, ntax);
            trlength = mfl_get_treelen(swapingon, tipdata, ntax, nchar);
            unroot(ntax, swapingon);
            if (trlength <= *currentbesttree) {
                *foundbettertree = true;
                *currentbesttree = trlength;
                treeset[*current] = copytree(swapingon, ntax, numnodes);
                treeset[*current]->index = *current;
                treeset[*current]->length = trlength;
                trlength = 0;
                *current = *current + 1;
            }
            mfl_bswap(n->next->next->outedge, n->outedge->next->outedge);
            //printf("current: %li\n", *current);
        }
        else {
            printf("Hit tree limit\n");
            *undertreelimit = false;
            return;
        }
    }
    
    p = n->next;
    while (p != n) {
        mfl_nni_traversal(p->outedge, swapingon, treeset, ntax, nchar, numnodes, 
                          current, tipdata, undertreelimit, currentbesttree,
                          foundbettertree);
        p = p->next;
    }
}


void mfl_nni_search(int ntax, int nchar, int numnodes, charstate *tipdata, 
                    tree **treeset, int starttreelen)
{
    long int i, currentbest;
    long int branch_count = 0;   // The branch iteration of the NNI traversal
    long int *bc_pointer = &branch_count, *currentbest_p = &currentbest;
    long int numreps = 1000;
    bool undertreelimit = true, *undertreelimit_p = &undertreelimit;
    bool foundbettertree = false, *foundbettertree_p = &foundbettertree;
    
    for (i = 0; i < numreps && *undertreelimit_p; ++i) {
        *foundbettertree_p = false;
        mfl_nni_traversal(treeset[i]->trnodes[0]->outedge, treeset[i], treeset, 
                          ntax, nchar, numnodes, bc_pointer, tipdata, 
                          undertreelimit_p, currentbest_p, foundbettertree_p);
        if (!foundbettertree_p) {
            freetree(treeset[i], numnodes);
            treeset[i] = randunrooted(ntax, numnodes);
        }
    }
    
    printf("bc_pointer: %li\n", *bc_pointer);
    
    for (i = 0; i < *bc_pointer; ++i) {
        printf("Tree %li:\n", i + 1);
        mfl_root_tree(treeset[i], treeset[i]->trnodes[6]->outedge->index, ntax);
        printNewick(treeset[i]->root);
        printf(";\n");
        printf("Length: %i\n", treeset[i]->length);
    }
    
    freetree(treeset[0], numnodes);
    free(treeset);
}