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

void mfl_nni_traversal(node *n, tree *swapingon, tree **treeset, int ntax, 
                       int nchar, int numnodes, long int *current, 
                       charstate *tipdata, bool *undertreelimit, 
                       long int *currentbesttree, bool *foundbettertree)
{
    int trlength = 0;
    node *p;
    
    if (n->start) {
        mfl_nni_traversal(n->outedge, swapingon, treeset, ntax, nchar, numnodes, 
                          current, tipdata, undertreelimit, currentbesttree,
                          foundbettertree);
        return;
    }
    
    if (n->tip || !(*undertreelimit)) {
        return;
    }
    
    if (!n->outedge->tip)
    {
        /* Nearly identical steps are conducted twice because all nearest-
         * neighbor interchanges produce two distinct tree topologies */
        
        if (*current + 1 <= TREELIMIT) {
            /* ^ This double condition might be preventing the last replicate */
            
            mfl_bswap(n->next->outedge, n->outedge->next->outedge);
            mfl_root_tree(swapingon, 1, ntax);
            trlength = mfl_get_treelen(swapingon, tipdata, ntax, nchar);
            unroot(ntax, swapingon);
            if (trlength < *currentbesttree) {
                *foundbettertree = true;
                *currentbesttree = trlength;
                swapingon->length = trlength;
                mfl_reinit_treebuffer(treeset, swapingon, current, numnodes);
                trlength = 0;
                *current = *current + 1;
                return;
            }
            if (trlength == *currentbesttree) {
                treeset[*current] = copytree(swapingon, ntax, numnodes);
                treeset[*current]->index = *current;
                treeset[*current]->length = trlength;
                trlength = 0;
                *current = *current + 1;
            }
            mfl_bswap(n->next->outedge, n->outedge->next->outedge);
            
            if (*current + 1 <= TREELIMIT){ 
                mfl_bswap(n->next->next->outedge, n->outedge->next->outedge);
                mfl_root_tree(swapingon, 1, ntax);
                trlength = mfl_get_treelen(swapingon, tipdata, ntax, nchar);
                unroot(ntax, swapingon);
                if (trlength < *currentbesttree) {
                    *foundbettertree = true;
                    *currentbesttree = trlength;
                    swapingon->length = trlength;
                    mfl_reinit_treebuffer(treeset, swapingon, current, numnodes);
                    trlength = 0;
                    *current = *current + 1;
                    return;
                }
                if (trlength == *currentbesttree) {
                    treeset[*current] = copytree(swapingon, ntax, numnodes);
                    treeset[*current]->index = *current;
                    treeset[*current]->length = trlength;
                    trlength = 0;
                    *current = *current + 1;
                }
                mfl_bswap(n->next->next->outedge, n->outedge->next->outedge);
            }
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
    long int i, j, currentbest;
    long int trbufpos = 0;   // The branch iteration of the NNI traversal
    long int *nxtintrbuf = &trbufpos, *currentbest_p = &currentbest;
    long int numreps = 1000;
    bool undertreelimit = true, *undertreelimit_p = &undertreelimit;
    bool foundbettertree = false, *foundbettertree_p = &foundbettertree;
    
    do {
        for (j = 0; j == 0 || j < *nxtintrbuf; ++j) {
            mfl_nni_traversal(treeset[j]->trnodes[0], treeset[j], treeset, 
                              ntax, nchar, numnodes, nxtintrbuf, tipdata, 
                              undertreelimit_p, currentbest_p, foundbettertree_p);
        }
        ++i;
    } while (i < numreps);
    
    printf("Next in trbuf: %li\n", *nxtintrbuf);
    
    printf("\nThe optimal tree(s) found by nearest neighbor interchange:\n");
    for (i = 0; i < *nxtintrbuf; ++i) {
        printf("Tree %li:\n", i + 1);
        mfl_root_tree(treeset[i], treeset[i]->trnodes[6]->outedge->index, ntax);
        printNewick(treeset[i]->root);
        printf(";\n");
        printf("Length: %i\n", treeset[i]->length);
    }
    printf("\n\n");
    mfl_clear_treebuffer(treeset, nxtintrbuf, numnodes);
    
    free(treeset);
}